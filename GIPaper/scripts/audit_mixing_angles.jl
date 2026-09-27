#!/usr/bin/env julia
# Audit of the 13 same-J mixing angles quoted in the text of Godfrey & Isgur,
# Phys. Rev. D 32, 189 (1985): strange (Fig. 4 caption), charmed (Fig. 7
# caption), and b-flavored (Fig. 9 caption) sectors.
#
# Paper convention (identical matrix printed in all three captions):
#
#   [ Q_low  ]   [  cos(theta_nL)  sin(theta_nL) ] [ n ^1L_L ]
#   [ Q_high ] = [ -sin(theta_nL)  cos(theta_nL) ] [ n ^3L_L ]
#
# so Q_low = cos(theta) |n ^1L_L> + sin(theta) |n ^3L_L>. The production
# spectrum diagonalizes all requested radial levels in one shared block. For a
# quoted nL pair this audit projects the lower of its two assigned physical
# states onto the matching singlet/triplet basis rows and evaluates
# theta = atan2(v_triplet, v_singlet). It also reports the norm outside that
# two-row projection instead of pretending the full block was 2x2.
#
# Angles are read from the public spectrum API only: compute_spectrum ->
# spectrum_state(...).mixings (mechanism "antisymmetric_spin_orbit").

using Pkg
Pkg.activate(@__DIR__)

using Dates
using CSV
using Printf

root = dirname(@__DIR__)
using GIModel
using GIPaper: quark_for   # paper flavor label -> quark object

params_path = default_parameters_path()
params, mq = load_parameters_and_quark_masses(params_path)

const L_OF = Dict("P" => 1, "D" => 2, "F" => 3, "G" => 4)

# Preserve the complete GI85 caption benchmark and add the independently
# sourced GK91 comparison without changing model inputs or fitting angles.
const targets = collect(CSV.File(joinpath(root, "data", "mixing_angle_targets.csv")))
original_targets = filter(r -> r.source_id == "GI1985", targets)
sectors = map(unique([(String(r.f1), String(r.f2)) for r in original_targets])) do (f1, f2)
    rows = filter(r -> r.f1 == f1 && r.f2 == f2, original_targets)
    (label = "$f1 $(f2)bar", f1 = Symbol(f1), f2 = Symbol(f2),
     rows = [(r.n, String(r.L), Float64(r.angle_deg)) for r in rows],
     src = String(first(rows).source_location))
end

# Angle quoted if the low/high labels were exchanged (theta defined mod 180;
# a label swap shifts by 90 within (-90, 90]).
complement_deg(theta) = theta > 0 ? theta - 90.0 : theta + 90.0

same_j_mix(state) =
    only(filter(m -> m.mechanism == "antisymmetric_spin_orbit", state.mixings))

results = NamedTuple[]
for sec in sectors
    levels = BasisState[]
    for (n, L, _) in sec.rows
        push!(levels, BasisState(n, L, 1, L_OF[L]))
        push!(levels, BasisState(n, L, 3, L_OF[L]))
    end
    meson = Meson(quark_for(mq, sec.f1), quark_for(mq, sec.f2))
    println("computing spectrum for ", flavor_label(meson), " ...")
    # Production solver settings (defaults): ngrid = 450, rmax = 24.0.
    spec = compute_spectrum(params, meson; levels = levels)
    for (n, L, paper) in sec.rows
        s1 = spectrum_state(spec, n, L, 1, L_OF[L])   # n ^1L_L
        s3 = spectrum_state(spec, n, L, 3, L_OF[L])   # n ^3L_L
        m1 = same_j_mix(s1)
        m3 = same_j_mix(s3)
        m1.result === m3.result || error("$n$L pair does not share one mixing result")
        basis = m1.result.block.basis
        is = findfirst(b -> b.label == s1.label, basis)
        it = findfirst(b -> b.label == s3.label, basis)
        (isnothing(is) || isnothing(it)) && error("$n$L pair absent from mixing basis")
        physical = s1.mass_GeV <= s3.mass_GeV ? m1 : m3
        cs, ct = physical.components[is], physical.components[it]
        phase = cs < 0 ? -1.0 : 1.0
        theta = atand(phase * ct, phase * cs)
        projection_norm2 = cs^2 + ct^2
        operator = spin_orbit_mixing_components(
            params,
            meson.constituent_masses,
            L,
            radial_wave(spec, s1.corrected),
            radial_wave(spec, s3.corrected),
        )
        operator_scale = abs(operator.vector) + abs(operator.thomas)
        cancellation = iszero(operator_scale) ? 0.0 : abs(operator.total) / operator_scale
        push!(results, (
            sector = sec.label,
            nL = string(n, L),
            paper = paper,
            theta = theta,
            delta = theta - paper,
            comp = complement_deg(theta),
            comp_delta = complement_deg(theta) - paper,
            diag_s = m1.unmixed_GeV,           # ^1L_L corrected (pre-mixing) mass
            diag_t = m3.unmixed_GeV,           # ^3L_L corrected (pre-mixing) mass
            split_MeV = 1000 * (m1.unmixed_GeV - m3.unmixed_GeV),
            offdiag_MeV = 1000 * m1.result.block.matrix[is, it],
            vector_MeV = 1000 * operator.vector,
            thomas_MeV = 1000 * operator.thomas,
            cancellation_percent = 100 * cancellation,
            low = min(s1.mass_GeV, s3.mass_GeV),
            high = max(s1.mass_GeV, s3.mass_GeV),
            outside_percent = max(0.0, 100 * (1 - projection_norm2)),
            src = sec.src,
        ))
    end
end

# One row per source/state: printed values are retained alongside convention-
# mapped targets and convention-independent projected singlet probabilities.
comparisons = map(targets) do target
    result = only(filter(r -> r.sector == "$(target.f1) $(target.f2)bar" &&
                             r.nL == "$(target.n)$(target.L)", results))
    (source_id = target.source_id, sector = target.sector, n = target.n, L = target.L,
     printed_angle_deg = target.printed_angle_deg, target_angle_deg = target.angle_deg,
     computed_angle_deg = result.theta,
     angular_distance_deg = abs(mod(result.theta - target.angle_deg + 90, 180) - 90),
     target_singlet_probability = cosd(target.angle_deg)^2,
     computed_singlet_probability = cosd(result.theta)^2,
     outside_percent = result.outside_percent,
     source_location = target.source_location, conversion = target.conversion)
end
CSV.write(joinpath(root, "docs", "residual_reports", "mixing_reference_comparison.csv"), comparisons)

# --- report ------------------------------------------------------------------

outpath = joinpath(root, "docs", "residual_reports", "mixing_angles.md")
open(outpath, "w") do io
    println(io, "# Same-J Mixing-Angle Audit (13 GI85 captions and 6 GK91 targets)")
    println(io)
    println(io, "Generated by `scripts/audit_mixing_angles.jl` on ",
        Dates.format(Dates.now(), "yyyy-mm-dd"), ".")
    println(io, numerics_provenance(FiniteDifferenceSolver()))
    println(io)
    println(io, "Computed-vs-paper comparison of every same-J `(^1L_L, ^3L_L)` mixing angle")
    println(io, "quoted in the text of Godfrey & Isgur, Phys. Rev. D 32, 189 (1985):")
    println(io, "strange (Fig. 4 caption), charmed (Fig. 7 caption), and b-flavored")
    println(io, "(Fig. 9 caption) sectors. Model angles come from the public spectrum API")
    println(io, "(`compute_spectrum` at production solver settings `ngrid = 450`,")
    println(io, "`rmax = 24.0`; `StateMixing` records with mechanism")
    println(io, "`antisymmetric_spin_orbit`), quark masses `u/d 0.220`, `s 0.419`,")
    println(io, "`c 1.628`, `b 4.977` GeV from `data/parameters.provisional.toml`.")
    println(io)
    println(io, "**Closed as an implementation investigation.** The missing all-L contact term is fixed. ",
        "The large GI85 caption-angle mismatch is a documented historical reference discrepancy: ",
        "the corrected calculation agrees with the later GI-model results. The origin of the ",
        "1985 values is unknown; no parameters are adjusted to match them.")
    println(io, "Targets and provenance: [ledger](../../data/mixing_angle_targets.md), ",
        "[machine-readable comparison](mixing_reference_comparison.csv). ",
        "See the [resolution](../investigations/mixing_composition_investigation.md).")
    println(io)
    println(io, "## Convention mapping (derived once, applied uniformly)")
    println(io)
    println(io, "The paper prints the same rotation matrix in all three captions:")
    println(io)
    println(io, "```")
    println(io, "[ Q_low  ]   [  cos(theta_nL)  sin(theta_nL) ] [ n ^1L_L ]")
    println(io, "[ Q_high ] = [ -sin(theta_nL)  cos(theta_nL) ] [ n ^3L_L ]")
    println(io, "```")
    println(io)
    println(io, "so the *lower-mass* eigenstate is")
    println(io, "`Q_low = cos(theta) |n ^1L_L> + sin(theta) |n ^3L_L>` with `cos(theta) >= 0`")
    println(io, "implied by the quoted range.")
    println(io)
    println(io, "The production spectrum diagonalizes one complete block containing every")
    println(io, "requested radial singlet and triplet with this `L,J`. For each quoted `nL`,")
    println(io, "we take the lower of the two physical states assigned to its singlet/triplet")
    println(io, "precursors, project its shared eigenvector onto those two exact basis rows,")
    println(io, "fix the singlet component positive, and report")
    println(io, "`theta = atan2(v_triplet,nL, v_singlet,nL)`. Thus the angular convention is")
    println(io, "identical to the paper while radial-state mixing remains present, so")
    println(io)
    println(io, "> **theta_paper = theta_model (identity angular mapping; the model value is a projection of the complete radial block).**")
    println(io)
    println(io, "Two residual sign conventions enter only through the off-diagonal element")
    println(io, "`c` and are pinned once, globally:")
    println(io)
    println(io, "1. **Constituent ordering.** `Meson(mq, f1, f2)` is the `f1 f2bar` meson;")
    println(io, "   the antisymmetric spin-orbit operator carries `(1/m1^2 - 1/m2^2)/2` and")
    println(io, "   flips sign under quark <-> antiquark exchange. Flavors are ordered")
    println(io, "   exactly as the paper's sector labels: `u sbar`, `c ubar`, `c sbar`,")
    println(io, "   `b ubar`, `b sbar`, `b cbar`.")
    println(io, "2. **Basis phase.** The angular matrix element is")
    println(io, "   `<^1L_L| L.(S1-S2) |^3L_L> = +sqrt(L(L+1))` with `S1` the quark (`f1`)")
    println(io, "   spin (`spin_orbit_mixing_components`, `src/spin_fine_structure.jl`).")
    println(io)
    println(io, "Consistency check (an outcome, not a fit): with this fixed mapping the")
    println(io, "model reproduces the paper's sector-wide sign pattern from the mass")
    println(io, "ordering alone — strange (`u` quark, heavy antiquark) angles positive,")
    println(io, "charmed (heavy quark, light antiquark) angles negative.")
    println(io)
    println(io, "## Computed vs paper")
    println(io)
    println(io, "Both columns use the paper convention above. `Delta = computed - paper`.")
    println(io)
    println(io, "| Sector | nL | Computed theta (deg) | Paper theta (deg) | Delta (deg) |")
    println(io, "|---|---|---:|---:|---:|")
    for r in results
        @printf(io, "| %s | %s | %+.1f | %+.0f | %+.1f |\n",
            r.sector, r.nL, r.theta, r.paper, r.delta)
    end
    println(io)
    println(io, "## Later GI-model reference: Godfrey–Kokoski (1991), Table I")
    println(io)
    println(io, "Same model parameters and lower-state convention. The printed K angle is -5 degrees ",
        "for s ubar; constituent exchange gives +5 degrees for our u sbar. ",
        "P_s is the singlet probability within the same-n projection. ",
        "The 1-degree printed step is precision, not an assigned uncertainty or a pass threshold.")
    println(io)
    println(io, "| Sector | GI85 angle | GK91 angle | Computed angle | Computed − GK91 | GK91 P_s | Computed P_s |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|")
    for c in filter(r -> r.source_id == "GK1991", comparisons)
        original = only(filter(r -> r.source_id == "GI1985" && r.sector == c.sector &&
                                    r.n == c.n && r.L == c.L, comparisons))
        @printf(io, "| %s | %+.0f | %+.0f | %+.1f | %+.1f | %.3f | %.3f |\n",
            c.sector, original.target_angle_deg, c.target_angle_deg, c.computed_angle_deg,
            c.computed_angle_deg - c.target_angle_deg, c.target_singlet_probability,
            c.computed_singlet_probability)
    end
    println(io)
    println(io, "No later target is substituted for a GI85 row. The seven caption rows without ",
        "a GK91 counterpart remain visible; missing later targets are not counted as agreement.")
    println(io)
    println(io, "## Reading the deviations")
    println(io)
    worst = sort(results; by = r -> abs(r.delta), rev = true)[1:min(5, length(results))]
    println(io, "The five largest direct-angle differences are:")
    for r in worst
        @printf(io, "- `%s %s`: model %+.1f deg, paper %+.0f deg (Delta %+.1f deg); ",
            r.sector, r.nL, r.theta, r.paper, r.delta)
        @printf(io, "same-n projection leaves %.2f%% norm in other radial rows.\n",
            r.outside_percent)
    end
    println(io)
    println(io, "Near-degenerate singlet/triplet diagonals remain sensitive to small matrix")
    println(io, "changes; when the physical low/high labeling is reversed, the complementary")
    println(io, "angle in the diagnostics is the relevant comparison. The explicit outside-")
    println(io, "projection percentage distinguishes that convention issue from genuine")
    println(io, "cross-radial composition in the complete block.")
    println(io)
    println(io, "The antisymmetric element is itself cancellation-sensitive: its vector")
    println(io, "color-magnetic and scalar/Thomas pieces usually have opposite signs.")
    println(io, "The survival percentage below is `|vector + Thomas| / (|vector| + |Thomas|)`.")
    println(io, "A small value indicates sensitivity to the relative strength of the two")
    println(io, "Appendix-A spin-orbit kernels; it does not establish an implementation error.")
    println(io)
    println(io, "## Diagnostics (raw block data)")
    println(io)
    println(io, "Corrected (pre-mixing) diagonal masses `E_s = E(^1L_L)`, `E_t = E(^3L_L)`,")
    println(io, "same-n off-diagonal `c = c_vector + c_Thomas`, the two assigned physical masses, projected angle,")
    println(io, "complementary angle (low/high labels exchanged), and norm carried by other")
    println(io, "radial rows. All angles already use the paper")
    println(io, "convention (identity mapping).")
    println(io)
    println(io, "| Sector | nL | E_s (GeV) | E_t (GeV) | E_s - E_t (MeV) | vector c (MeV) | Thomas c (MeV) | total c (MeV) | survives (%) | low (GeV) | high (GeV) | theta (deg) | complement (deg) | other radial norm (%) |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in results
        @printf(io,
            "| %s | %s | %.4f | %.4f | %+.2f | %+.2f | %+.2f | %+.2f | %.1f | %.4f | %.4f | %+.2f | %+.2f | %.2f |\n",
            r.sector, r.nL, r.diag_s, r.diag_t, r.split_MeV,
            r.vector_MeV, r.thomas_MeV, r.offdiag_MeV,
            r.cancellation_percent, r.low, r.high, r.theta, r.comp,
            r.outside_percent)
    end
    println(io)
    println(io, "## Paper sources")
    println(io)
    for sec in sectors
        println(io, "- `", sec.label, "`: ", sec.src, ".")
    end
end

println("wrote ", outpath)
println()
@printf("%-8s %-3s %10s %8s %8s %12s\n",
    "sector", "nL", "computed", "paper", "delta", "complement")
for r in results
    @printf("%-8s %-3s %+10.1f %+8.0f %+8.1f %+12.1f\n",
        r.sector, r.nL, r.theta, r.paper, r.delta, r.comp)
end
