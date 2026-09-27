# Residual markdown reports and mixing-prone scorecard helpers, ported verbatim
# from the pre-split GIModel sector_comparison.jl. Rows are the NamedTuples
# produced by compare_reference.
#
# Public API (exported from GIPaper.jl): mixing_prone_state,
#   nonmixing_deviation_summary, write_residual_report


function _row_sector(row)
    return String(getproperty(row, :sector))
end

function _row_L(row)
    return String(getproperty(row, :L))
end

function _same_j_singlet_triplet_prone(row)
    Ls = _row_L(row)
    Lval = try
        orbital_angular_momentum(Ls)
    catch error
        error isa ArgumentError || rethrow()
        return false
    end
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
mixing machinery. It excludes isoscalar flavor-mixing rows, unassigned same-`J`
`^1L_J`/`^3L_J` candidates, and the common triplet `S`/`D`, `J=1` tensor/radial
mixing candidates. Open-flavor same-`J` rows become scoreable once the
antisymmetric spin-orbit block has been applied.
"""
function mixing_prone_state(row)
    sector = _row_sector(row)
    if sector == "isoscalar"
        scheme =
            hasproperty(row, :isoscalar_annihilation_scheme) ?
            String(getproperty(row, :isoscalar_annihilation_scheme)) : "none"
        # Ideal-mixing and general-Eq.(16) rows are genuine predictions; the
        # calibrated pseudoscalar control is a fit and unassigned rows still
        # lack their flavor partner, so both stay excluded.
        scheme in ("ideal", "general_eq16") || return true
        return _tensor_sd_prone(row)
    end
    if hasproperty(row, :same_j_mixing_scheme) &&
       getproperty(row, :same_j_mixing_scheme) == "antisymmetric_spin_orbit"
        return _tensor_sd_prone(row)
    end
    if hasproperty(row, :tensor_mixing_scheme) &&
       getproperty(row, :tensor_mixing_scheme) == "tensor_mixing"
        return false
    end
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

# Report-header sentence for each central-potential construction.
_central_note(::AppendixAMomentumSandwich) =
    "closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, "
_central_note(::AppendixAClosedForm) =
    "closed-form Gaussian-smeared GI G̃(r) and S̃(r), without the central Coulomb momentum sandwich, "
_central_note(::AppendixADerivativeG) =
    "Appendix-A derivative proxy for G(r), `G + ∇²G/(4σ²)`, with pointwise S(r); older comparator below the closed-form (A12)–(A14) modes, "
_central_note(::AppendixASmearing3D) =
    "experimental (A7)–(A8)-style 3D isotropic smearing of pointwise Coulomb G and confinement S (Table II σ₀, s), "
_central_note(::Coulomb1DSmearing) =
    "1D Gaussian renormalization of G(r) only (pointwise S); same σ as contact (A9); diagnostic comparator below the closed-form (A12)–(A14) modes, "
_central_note(::PointwiseCentral) =
    "pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), "

function write_residual_report(
    path::AbstractString,
    title::AbstractString,
    rows;
    solver::RadialSolver = FiniteDifferenceSolver(),
    contact_hyperfine::Bool = true,
    central::CentralPotentialMethod = PointwiseCentral(),
    contact_momentum_sandwich::Bool = false,
    fine_structure_momentum_sandwich::Bool = false,
    fine_structure_smeared_kernels::Bool = false,
    use_fine_structure::Bool = true,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io)
        hyperfine_note = if contact_hyperfine && contact_momentum_sandwich
            "with GI momentum-sandwiched smeared contact hyperfine (every L)"
        elseif contact_hyperfine
            "with smeared contact hyperfine (every L)"
        else
            "without contact hyperfine"
        end
        fs_note =
            if use_fine_structure &&
               fine_structure_momentum_sandwich &&
               fine_structure_smeared_kernels
                " fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; "
            elseif use_fine_structure && fine_structure_momentum_sandwich
                " fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches; "
            elseif use_fine_structure && fine_structure_smeared_kernels
                " fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with smeared-G/S derivative kernels; "
            elseif use_fine_structure
                " fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization; "
            else
                " no L·S/tensor; "
            end
        central_note = _central_note(central)
        println(
            io,
            "Model: $hyperfine_note,$fs_note",
            "GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`;",
            " ",
            central_note,
            "see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.",
        )
        println(io)
        println(io, numerics_provenance(solver))
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
                "The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.",
            )
            println(io)
            println(
                io,
                "| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |",
            )
            println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
            for row in rows
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                annihilation_shift =
                    hasproperty(row, :annihilation_shift_GeV) ? row.annihilation_shift_GeV : 0.0
                total_shift =
                    row.contact_shift_GeV + row.fine_structure_shift_GeV + annihilation_shift
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.3f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %.3f |",
                        label,
                        row.central_GeV,
                        1000 * row.contact_shift_GeV,
                        1000 * row.spin_orbit_vector_shift_GeV,
                        1000 * row.spin_orbit_thomas_shift_GeV,
                        1000 * row.spin_orbit_shift_GeV,
                        1000 * row.tensor_shift_GeV,
                        1000 * annihilation_shift,
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
                "Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.",
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

        if !isempty(rows) &&
           hasproperty(rows[1], :same_j_mixing_scheme) &&
           any(row.same_j_mixing_scheme == "antisymmetric_spin_orbit" for row in rows)
            println(io)
            println(io, "## Same-J Antisymmetric Spin-Orbit Mixing")
            println(io)
            println(
                io,
                "Rows below use the mixed eigenvalues from the `(^1L_J, ^3L_J)` mass block. Components are ordered as singlet/triplet in the unmixed basis.",
            )
            println(io)
            println(
                io,
                "| state | unmixed GeV | mixed GeV | offdiag MeV | theta deg | singlet component | triplet component |",
            )
            println(io, "|---|---:|---:|---:|---:|---:|---:|")
            for row in rows
                row.same_j_mixing_scheme == "antisymmetric_spin_orbit" || continue
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.3f | %.3f | %+7.1f | %+7.2f | %+7.3f | %+7.3f |",
                        label,
                        row.same_j_unmixed_GeV,
                        row.predicted_GeV,
                        1000 * row.same_j_offdiag_GeV,
                        row.same_j_mixing_angle_deg,
                        row.same_j_component_singlet,
                        row.same_j_component_triplet,
                    ),
                )
            end
        end

        if !isempty(rows) &&
           hasproperty(rows[1], :tensor_mixing_scheme) &&
           any(row.tensor_mixing_scheme == "tensor_mixing" for row in rows)
            println(io)
            println(io, "## Same-J Tensor Mixing")
            println(io)
            println(
                io,
                "Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.",
            )
            println(io)
            println(
                io,
                "| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |",
            )
            println(io, "|---|---:|---:|---:|---:|---:|")
            for row in rows
                row.tensor_mixing_scheme == "tensor_mixing" || continue
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.3f | %.3f | %+7.1f | %+7.3f | %+7.3f |",
                        label,
                        row.tensor_unmixed_GeV,
                        row.predicted_GeV,
                        1000 * row.tensor_offdiag_GeV,
                        row.tensor_component_lowL,
                        row.tensor_component_highL,
                    ),
                )
            end
        end

        if !isempty(rows) &&
           hasproperty(rows[1], :isoscalar_annihilation_scheme) &&
           any(row.isoscalar_annihilation_scheme != "none" for row in rows)
            println(io)
            println(io, "## Isoscalar Annihilation Mixing")
            println(io)
            println(
                io,
                "Rows below use the selected isoscalar annihilation scheme. `calibrated_p1` is the Fig. 5/Table III pseudoscalar control; `p1`/`p2` are the paper Eq. (18a,b) formula modes; `general_eq16` is the literal Eq. (16) block (`^3S_1`, `^3P_2`); `ideal` assigns the heavier row to the `s sbar` channel prediction per the Table III ideal-mixing prescription.",
            )
            println(io)
            println(io, "| state | scheme | unmixed GeV | mixed GeV | shift MeV |")
            println(io, "|---|---|---:|---:|---:|")
            for row in rows
                row.isoscalar_annihilation_scheme == "none" && continue
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                println(
                    io,
                    @sprintf(
                        "| `%s` | `%s` | %.3f | %.3f | %+7.1f |",
                        label,
                        row.isoscalar_annihilation_scheme,
                        row.isoscalar_annihilation_unmixed_GeV,
                        row.predicted_GeV,
                        1000 * row.annihilation_shift_GeV,
                    ),
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
        if central isa AppendixAMomentumSandwich
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
