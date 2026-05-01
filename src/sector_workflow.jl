# Public API (exported from GIModel.jl):
#   solve_sector, compute_sector, compare, write_residual_report
#   RadialChannelKey, ChannelRadialSolution, SectorComputation

"""
    RadialChannelKey(m1_GeV, m2_GeV, L_label)

Dict key for one **spin-independent radial channel**: constituent quark masses (GeV) and
orbital label (`S`, `P`, …), same convention as `ReferenceState.L`. Masses are rounded
to 12 significant figures so distinct finite-difference solves that would share the same
physics collapse to one cache entry.
"""
struct RadialChannelKey
    m1_GeV::Float64
    m2_GeV::Float64
    L_label::String
    function RadialChannelKey(m1::Real, m2::Real, L_label::AbstractString)
        new(
            round(Float64(m1); sigdigits = 12),
            round(Float64(m2); sigdigits = 12),
            String(L_label),
        )
    end
end

"""
    ChannelRadialSolution(eigenvalues_GeV, eigenvectors, r)

Output of one `channel_solution` call, stored in `SectorComputation.channel_cache`:

  - `eigenvalues_GeV`: lowest radial eigenvalues (GeV) of the central Hamiltonian on the mesh.
  - `eigenvectors`: columns are reduced radial functions ``u_n(r)`` for each level.
  - `r`: uniform interior radial grid (same spacing as in the FD builder).
"""
struct ChannelRadialSolution
    eigenvalues_GeV::Vector{Float64}
    eigenvectors::Matrix{Float64}
    r::Vector{Float64}
end

"""
    SectorComputation(params, m_fallback, channel_cache)

Heavy lifting from `compute_sector`: precomputed radial FD solves per distinct channel.

# Fields

  - `params`: `GIParameters` used to build each central Hamiltonian.
  - `m_fallback`: mass (GeV) used when `parse_quark_masses` fails for a row
    (typically `params.masses[flavor]` for the `flavor` passed to `compute_sector`).
  - `channel_cache`: map `RadialChannelKey` → `ChannelRadialSolution`.

`channel_cache` **deduplicates** work: every reference state with the same rounded
`(m₁, m₂, L)` shares one eigenproblem. `compare` then picks radial level `n`,
applies contact hyperfine and fine-structure corrections, and forms residuals — that part
depends on `J`, multiplicity, etc., and is not stored here.
"""
struct SectorComputation
    params::GIParameters
    m_fallback::Float64
    channel_cache::Dict{RadialChannelKey,ChannelRadialSolution}
end

function _reference_masses(params::GIParameters, m_fallback::Float64, s::ReferenceState)
    try
        return parse_quark_masses(params, String(s.sector), String(s.quark_content))
    catch
        return m_fallback, m_fallback
    end
end

function solve_sector(
    params::GIParameters,
    flavor::String;
    maxn::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    m = params.masses[flavor]
    results = Dict{Tuple{Int,String},Float64}()
    for (symbol, L) in L_SYMBOLS
        levels = solve_channel(
            params,
            m,
            m,
            L;
            nlevels = maxn,
            ngrid = ngrid,
            rmax = rmax,
            kinetic = kinetic,
            eigensolver = eigensolver,
        )
        for n = 1:length(levels)
            results[(n, symbol)] = levels[n]
        end
    end
    results
end

function compute_sector(
    params::GIParameters,
    reference::Vector{ReferenceState},
    flavor::String;
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    m_fallback = params.masses[flavor]
    channel_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    for state in reference
        m1, m2 = _reference_masses(params, m_fallback, state)
        key = RadialChannelKey(m1, m2, state.L)
        if !haskey(channel_cache, key)
            Lval = L_SYMBOLS[state.L]
            ev, vecs, r = channel_solution(
                params,
                m1,
                m2,
                Lval;
                nlevels = 6,
                ngrid = ngrid,
                rmax = rmax,
                kinetic = kinetic,
                eigensolver = eigensolver,
            )
            channel_cache[key] = ChannelRadialSolution(ev, vecs, r)
        end
    end
    SectorComputation(params, m_fallback, channel_cache)
end

function compare(
    computed::SectorComputation,
    reference::Vector{ReferenceState};
    contact_hyperfine::Bool = true,
    use_fine_structure::Bool = true,
)
    params = computed.params
    m_fallback = computed.m_fallback
    channel_cache = computed.channel_cache
    rows = NamedTuple[]
    for state in reference
        m1, m2 = _reference_masses(params, m_fallback, state)
        key = RadialChannelKey(m1, m2, state.L)
        haskey(channel_cache, key) || continue
        sol = channel_cache[key]
        ev = sol.eigenvalues_GeV
        vectors = sol.eigenvectors
        r = sol.r
        state.n <= length(ev) || continue
        h = r[2] - r[1]
        central = ev[state.n]
        contact_shift = 0.0
        spin_orbit_vector_shift = 0.0
        spin_orbit_thomas_shift = 0.0
        spin_orbit_shift = 0.0
        tensor_shift = 0.0
        fine_structure_shift = 0.0
        fine_structure_mass_convention = "disabled"
        if contact_hyperfine
            contact_shift = contact_hyperfine_shift_active(
                params,
                m1,
                m2,
                state.L,
                state.multiplicity,
                vectors[:, state.n],
                r,
            )
        end
        if use_fine_structure && params.fine_structure
            comp = fine_structure_components(
                params,
                m1,
                m2,
                state.L,
                state.multiplicity,
                state.J,
                collect(vectors[:, state.n]),
                collect(r),
                h;
                enabled = true,
                k_spin_orbit = params.k_spin_orbit,
                k_tensor = params.k_tensor,
            )
            spin_orbit_vector_shift = comp.spin_orbit_vector
            spin_orbit_thomas_shift = comp.spin_orbit_thomas
            spin_orbit_shift = comp.spin_orbit
            tensor_shift = comp.tensor
            fine_structure_shift = comp.total
            fine_structure_mass_convention =
                isapprox(m1, m2; rtol = 0.0, atol = 0.0) ? "equal_mass" :
                "unequal_mass_equal_share_LdotS"
        end
        predicted = central + contact_shift + fine_structure_shift
        push!(
            rows,
            (
                sector = state.sector,
                state = state.composition,
                L = state.L,
                n = state.n,
                J = state.J,
                multiplicity = state.multiplicity,
                m1_GeV = m1,
                m2_GeV = m2,
                fine_structure_mass_convention = fine_structure_mass_convention,
                reference_GeV = state.mass_GeV,
                central_GeV = central,
                contact_shift_GeV = contact_shift,
                spin_orbit_vector_shift_GeV = spin_orbit_vector_shift,
                spin_orbit_thomas_shift_GeV = spin_orbit_thomas_shift,
                spin_orbit_shift_GeV = spin_orbit_shift,
                tensor_shift_GeV = tensor_shift,
                fine_structure_shift_GeV = fine_structure_shift,
                predicted_GeV = predicted,
                residual_MeV = 1000 * (predicted - state.mass_GeV),
                confidence = state.confidence,
            ),
        )
    end
    rows
end

function write_residual_report(
    path::AbstractString,
    title::AbstractString,
    rows;
    kinetic::Symbol = :relativistic,
    contact_hyperfine::Bool = true,
    appendix_a_smearing::Bool = false,
    appendix_a_derivative_g::Bool = false,
    appendix_a_closed_form::Bool = false,
    appendix_a_momentum_sandwich::Bool = false,
    contact_momentum_sandwich::Bool = false,
    fine_structure_momentum_sandwich::Bool = false,
    fine_structure_smeared_kernels::Bool = false,
    coulomb_1d_smear::Bool = false,
    use_fine_structure::Bool = true,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io)
        hyperfine_note = if contact_hyperfine && contact_momentum_sandwich
            "with GI momentum-sandwiched smeared S-wave contact hyperfine"
        elseif contact_hyperfine
            "with smeared S-wave contact hyperfine"
        else
            "without S-wave contact hyperfine"
        end
        fs_note =
            if use_fine_structure &&
               fine_structure_momentum_sandwich &&
               fine_structure_smeared_kernels
                " first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; "
            elseif use_fine_structure && fine_structure_momentum_sandwich
                " first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches; "
            elseif use_fine_structure && fine_structure_smeared_kernels
                " first-order L·S (vector+Thomas) and OGE-tensor with smeared-G/S derivative kernels; "
            elseif use_fine_structure
                " first-order L·S (vector+Thomas) and OGE-tensor; "
            else
                " no first-order L·S/tensor; "
            end
        central_note = if appendix_a_momentum_sandwich
            "closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, "
        elseif appendix_a_closed_form
            "closed-form Gaussian-smeared GI G̃(r) and S̃(r), without the central Coulomb momentum sandwich, "
        elseif appendix_a_derivative_g
            "Appendix-A derivative proxy for G(r), `G + ∇²G/(4σ²)`, with pointwise S(r); not the full audited (A12)–(A13) expansion, "
        elseif appendix_a_smearing
            "experimental (A7)–(A8)-style 3D isotropic smearing of pointwise Coulomb G and confinement S (Table II σ₀, s), "
        elseif coulomb_1d_smear
            "1D Gaussian renormalization of G(r) only (pointwise S); same σ as contact (A9); not the full (A12)–(A13) expansion, "
        else
            "pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), "
        end
        println(
            io,
            "Model: finite-difference + `$kinetic` kinetic, $hyperfine_note,$fs_note",
            "GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`;",
            " ",
            central_note,
            "see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).",
        )
        println(io)
        println(io, "| state | reference GeV | baseline GeV | residual MeV | confidence |")
        println(io, "|---|---:|---:|---:|---|")
        for row in rows
            label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
            println(
                io,
                @sprintf(
                    "| `%s` | %.3f | %.3f | %+7.1f | %s |",
                    label,
                    row.reference_GeV,
                    row.predicted_GeV,
                    row.residual_MeV,
                    row.confidence,
                ),
            )
        end

        if !isempty(rows) && hasproperty(rows[1], :central_GeV)
            println(io)
            println(io, "## Contribution Breakdown (diagnostic)")
            println(io)
            println(
                io,
                "All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).",
            )
            println(io)
            println(
                io,
                "| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |",
            )
            println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for row in rows
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                total_shift = row.contact_shift_GeV + row.fine_structure_shift_GeV
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.3f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %.3f |",
                        label,
                        row.central_GeV,
                        1000 * row.contact_shift_GeV,
                        1000 * row.spin_orbit_vector_shift_GeV,
                        1000 * row.spin_orbit_thomas_shift_GeV,
                        1000 * row.spin_orbit_shift_GeV,
                        1000 * row.tensor_shift_GeV,
                        1000 * total_shift,
                        row.predicted_GeV,
                    ),
                )
            end
        end

        if use_fine_structure &&
           !isempty(rows) &&
           hasproperty(rows[1], :fine_structure_mass_convention) &&
           any(!isapprox(row.m1_GeV, row.m2_GeV; rtol = 0.0, atol = 0.0) for row in rows)
            println(io)
            println(io, "## Fine-Structure Mass Convention (audit note)")
            println(io)
            println(
                io,
                "Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).",
            )
            println(io)
            println(io, "| state | m1 GeV | m2 GeV | convention |")
            println(io, "|---|---:|---:|---|")
            for row in rows
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.6f | %.6f | `%s` |",
                        label,
                        row.m1_GeV,
                        row.m2_GeV,
                        row.fine_structure_mass_convention
                    )
                )
            end
        end

        residuals = [abs(row.residual_MeV) for row in rows]
        if !isempty(residuals)
            println(io)
            println(
                io,
                @sprintf(
                    "Mean absolute residual: %.1f MeV.",
                    sum(residuals) / length(residuals)
                )
            )
            println(io, @sprintf("Max absolute residual: %.1f MeV.", maximum(residuals)))
        end
        groups = Dict{Tuple{Int,String},Vector{eltype(rows)}}()
        for row in rows
            push!(get!(groups, (row.n, row.L), eltype(rows)[]), row)
        end
        println(io)
        println(io, "## Spin-Averaged Diagnostics")
        println(io)
        println(io, "Weighted by `2J+1` within each available `(n, L)` group.")
        println(io)
        println(io, "| multiplet | states | reference GeV | baseline GeV | residual MeV |")
        println(io, "|---|---:|---:|---:|---:|")
        for key in sort(collect(keys(groups)); by = x -> (x[2], x[1]))
            group = groups[key]
            weights = [2 * row.J + 1 for row in group]
            weight_sum = sum(weights)
            ref =
                sum(w * row.reference_GeV for (w, row) in zip(weights, group)) / weight_sum
            pred =
                sum(w * row.predicted_GeV for (w, row) in zip(weights, group)) / weight_sum
            println(
                io,
                @sprintf(
                    "| `%d%s` | %d | %.3f | %.3f | %+7.1f |",
                    key[1],
                    key[2],
                    length(group),
                    ref,
                    pred,
                    1000 * (pred - ref)
                )
            )
        end
        println(io)
        if appendix_a_momentum_sandwich
            println(
                io,
                "The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.",
            )
        else
            println(
                io,
                "This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.",
            )
        end
    end
end
