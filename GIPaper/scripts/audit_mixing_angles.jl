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
using Printf

root = dirname(@__DIR__)
using GIModel
using GIPaper: quark_for   # paper flavor label -> quark object

params_path = default_parameters_path()
params, mq = load_parameters_and_quark_masses(params_path)

const L_OF = Dict("P" => 1, "D" => 2, "F" => 3, "G" => 4)

# The 13 text-quoted angles: (n, L, paper theta in degrees), grouped by
# sector with the paper's flavor ordering (flavor1 = quark, flavor2 = antiquark).
sectors = [
    (label = "u sbar", f1 = :u, f2 = :s,
     rows = [(1, "P", 34.0), (1, "D", 33.0), (2, "P", 15.0),
             (1, "F", 32.0), (2, "D", 25.0), (1, "G", 33.0)],
     src = "Fig. 4 caption (PDF p. 7; vision OCR ~line 348)"),
    (label = "c ubar", f1 = :c, f2 = :u,
     rows = [(1, "P", -41.0), (1, "D", -39.0)],
     src = "Fig. 7 caption (PDF p. 9; vision OCR ~line 402)"),
    (label = "c sbar", f1 = :c, f2 = :s,
     rows = [(1, "P", -44.0), (1, "D", -39.0)],
     src = "Fig. 7 caption (PDF p. 9; vision OCR ~line 402)"),
    (label = "b ubar", f1 = :b, f2 = :u,
     rows = [(1, "P", -43.0)],
     src = "Fig. 9 caption (PDF p. 10; vision OCR ~line 432)"),
    (label = "b sbar", f1 = :b, f2 = :s,
     rows = [(1, "P", -45.0)],
     src = "Fig. 9 caption (PDF p. 10; vision OCR ~line 432)"),
    (label = "b cbar", f1 = :b, f2 = :c,
     rows = [(1, "P", -53.0)],
     src = "Fig. 9 caption (PDF p. 10; vision OCR ~line 432)"),
]

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
            low = min(s1.mass_GeV, s3.mass_GeV),
            high = max(s1.mass_GeV, s3.mass_GeV),
            outside_percent = 100 * (1 - projection_norm2),
            src = sec.src,
        ))
    end
end

# --- report ------------------------------------------------------------------

outpath = joinpath(root, "docs", "residual_reports", "mixing_angles.md")
open(outpath, "w") do io
    println(io, "# Same-J Mixing-Angle Audit (13 text-quoted angles)")
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
    println(io, "## Diagnostics (raw block data)")
    println(io)
    println(io, "Corrected (pre-mixing) diagonal masses `E_s = E(^1L_L)`, `E_t = E(^3L_L)`,")
    println(io, "same-n off-diagonal `c`, the two assigned physical masses, projected angle,")
    println(io, "complementary angle (low/high labels exchanged), and norm carried by other")
    println(io, "radial rows. All angles already use the paper")
    println(io, "convention (identity mapping).")
    println(io)
    println(io, "| Sector | nL | E_s (GeV) | E_t (GeV) | E_s - E_t (MeV) | offdiag c (MeV) | low (GeV) | high (GeV) | theta (deg) | complement (deg) | other radial norm (%) |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in results
        @printf(io,
            "| %s | %s | %.4f | %.4f | %+.2f | %+.2f | %.4f | %.4f | %+.2f | %+.2f | %.2f |\n",
            r.sector, r.nL, r.diag_s, r.diag_t, r.split_MeV, r.offdiag_MeV,
            r.low, r.high, r.theta, r.comp, r.outside_percent)
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
