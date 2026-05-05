# Cached radial solves (`compute_sector`), reference-row comparison (`compare`),
# markdown reports (`write_residual_report`).
#
# Public API (exported from GIModel.jl): compute_sector, compare,
#   mixing_prone_state, nonmixing_deviation_summary, write_residual_report

"""
    compute_sector(params, annotated::AbstractVector; …)

Run [`channel_solution`](@ref) once per distinct [`RadialChannelKey`](@ref) needed by `annotated`,
cache eigenpairs in [`SectorComputation`](@ref), return it for [`compare`](@ref).

`annotated` is normally `Vector{ReferenceStateWithMasses}` from [`attach_constituent_masses`](@ref); each
element must have `.constituent_masses` and `.state` (with `.L`, `.n`, …) like [`ReferenceStateWithMasses`](@ref).
"""
function compute_sector(
    params::GIParameters,
    annotated::AbstractVector;
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    channel_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    for row in annotated
        masses = row.constituent_masses
        state = row.state
        key = RadialChannelKey(masses, state.L)
        if !haskey(channel_cache, key)
            Lval = L_SYMBOLS[state.L]
            ev, vecs, r = channel_solution(
                params,
                masses,
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
    return SectorComputation(params, channel_cache)
end

function compare(
    computed::SectorComputation,
    annotated::AbstractVector;
    contact_hyperfine::Bool = true,
    use_fine_structure::Bool = true,
)
    params = computed.params
    channel_cache = computed.channel_cache
    rows = NamedTuple[]
    for row in annotated
        state = row.state
        masses = row.constituent_masses
        key = RadialChannelKey(masses, state.L)
        haskey(channel_cache, key) || continue
        sol = channel_cache[key]
        ev = sol.eigenvalues_GeV
        state.n <= length(ev) || continue
        central = ev[state.n]
        contact_shift = 0.0
        spin_orbit_vector_shift = 0.0
        spin_orbit_thomas_shift = 0.0
        spin_orbit_shift = 0.0
        tensor_shift = 0.0
        fine_structure_shift = 0.0
        fine_structure_mass_convention = "disabled"
        multiplet = FineStructureMultiplet(state)
        wave = RadialWaveOnUniformMesh(sol, state.n)
        if contact_hyperfine
            contact_shift = contact_hyperfine_shift_active(params, masses, multiplet, wave)
        end
        if use_fine_structure && params.fine_structure
            comp = fine_structure_components(
                params,
                masses,
                multiplet,
                wave;
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
                isapprox(masses.m1_GeV, masses.m2_GeV; rtol = 0.0, atol = 0.0) ?
                "equal_mass" : "unequal_mass_equal_share_LdotS"
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
                m1_GeV = masses.m1_GeV,
                m2_GeV = masses.m2_GeV,
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
    return rows
end

function _row_sector(row)
    return String(getproperty(row, :sector))
end

function _row_L(row)
    return String(getproperty(row, :L))
end

function _same_j_singlet_triplet_prone(row)
    Ls = _row_L(row)
    haskey(L_SYMBOLS, Ls) || return false
    Lval = L_SYMBOLS[Ls]
    return Lval > 0 && getproperty(row, :J) == Lval && getproperty(row, :multiplicity) in (1, 3)
end

function _tensor_sd_prone(row)
    return getproperty(row, :multiplicity) == 3 &&
           getproperty(row, :J) == 1 &&
           _row_L(row) in ("S", "D") &&
           (getproperty(row, :n) > 1 || _row_L(row) == "D")
end

function _open_flavor_like(row)
    if hasproperty(row, :m1_GeV) && hasproperty(row, :m2_GeV)
        return !isapprox(
            getproperty(row, :m1_GeV),
            getproperty(row, :m2_GeV);
            rtol = 0.0,
            atol = 1.0e-12,
        )
    end
    return _row_sector(row) in (
        "strange",
        "charmed",
        "charmed_strange",
        "bottom_light",
        "bottom_strange",
        "bottom_charm",
        "b_flavored",
    )
end

"""
    mixing_prone_state(row) -> Bool

Heuristic guardrail for residual scorecards that should not depend on explicit
mixing machinery. It excludes isoscalar flavor-mixing rows, same-`J`
`^1L_J`/`^3L_J` candidates, and the common triplet `S`/`D`, `J=1` tensor/radial
mixing candidates.
"""
function mixing_prone_state(row)
    sector = _row_sector(row)
    sector == "isoscalar" && return true
    if _open_flavor_like(row)
        _same_j_singlet_triplet_prone(row) && return true
    end
    return _tensor_sd_prone(row)
end

function nonmixing_deviation_summary(rows)
    sectors = sort(unique(_row_sector(row) for row in rows))
    out = NamedTuple[]
    for sector in sectors
        sector_rows = [row for row in rows if _row_sector(row) == sector]
        kept = [row for row in sector_rows if !mixing_prone_state(row)]
        absres = [abs(getproperty(row, :residual_MeV)) for row in kept]
        push!(
            out,
            (
                sector = sector,
                n_included = length(kept),
                n_excluded_mixing_prone = length(sector_rows) - length(kept),
                mean_abs_deviation_MeV = isempty(absres) ? NaN : sum(absres) / length(absres),
                max_abs_deviation_MeV = isempty(absres) ? NaN : maximum(absres),
            ),
        )
    end
    return out
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
