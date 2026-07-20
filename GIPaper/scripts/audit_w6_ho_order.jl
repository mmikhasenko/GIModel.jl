#!/usr/bin/env julia
# =============================================================================
# W6 audit — paper-order (HO first-order) spin-distorted waves vs Table VII.
# =============================================================================
# Resolves the residual structure left by the central-wave Table VII gluonic
# audit (eta_c 0.87 vs psi 1.11, chi_0c 0.78 vs chi_2c 1.14, Upsilon tower
# ~1.1x). Three findings, each scored below:
#
#  1. NOT basis fidelity: the paper-style HO diagonalization of the CENTRAL
#     Hamiltonian reproduces the FD smeared wavefunction-at-origin S_L to
#     ≤0.5% (charm) / ≤3% (bottom) — far too small, and in the wrong
#     direction, to explain the 15-20% row residuals.
#
#  2. The residuals are SPIN-DEPENDENT WAVEFUNCTION DISTORTION: the smeared
#     contact term (attractive for ¹S₀, repulsive for ³S₁) and the spin-orbit/
#     tensor terms (attractive for ³P₀, repulsive for ³P₂) shift the paper's
#     wavefunction at the origin away from the central wave, with exactly the
#     observed signs on all five splitting patterns.
#
#  3. The paper's treatment is FIRST ORDER in the HO eigenbasis: resumming the
#     distortion nonperturbatively (FD) overshoots (the wave collapses into
#     the smeared attraction; e.g. eta_c 0.87 -> 1.11, chi_0b 0.81 -> 1.32),
#     while first-order perturbation theory in the HO central eigenbasis with
#     the calibrated k_spin_orbit/k_tensor scales lands all 16 gluonic rows in
#     [0.92, 1.09] (median ~1.02).
#
# Wave sources: contact_hyperfine_operator (S-waves, mult 1/3),
# fine_structure_grid_operator (3P_J), ho_first_order_distorted_states
# (paper-order treatment), contact_hyperfine_nonperturbative_states /
# H_central + V (nonperturbative FD reference).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using Printf
using LinearAlgebra
using GIModel

const ROOT = dirname(@__DIR__)                                  # GIPaper/
const PARAMS_PATH = joinpath(dirname(ROOT), "data", "parameters.provisional.toml")
const REPORT = joinpath(ROOT, "docs", "residual_reports", "w6_ho_order_validation.md")

const NGRID = 1200
const RMAX = 24.0
const NPTS = 900
const NB = 24

params, mq = load_parameters_and_quark_masses(PARAMS_PATH)
params_ho = with_basis(params, HarmonicOscillatorBasis)

# Phase convention: outermost antinode positive (shared with audit_table_vii.jl).
function fix_outer_antinode_positive!(u)
    peak = maximum(abs, u)
    i = findlast(x -> abs(x) > 0.2 * peak, u)
    (i !== nothing && u[i] < 0) && (u .*= -1)
    return u
end

smeared_S(u, r, mQ, L) = wavefunction_origin_smearing(
    RadialWaveOnUniformMesh(fix_outer_antinode_positive!(copy(u)), r), mQ; L = L, npoints = NPTS)

# --- spin-dependent grid operator for a sector ------------------------------
# L=0: smeared contact (multiplicity 1 or 3); L=1: triplet 3P_J LS+tensor.
function spin_operator(masses, L, mult, J, r, h)
    L == 0 && return GIModel.contact_hyperfine_operator(params, masses, "S", mult, r)
    return fine_structure_grid_operator(params, masses, J, r, h; L = L)
end

# --- the three treatments of one (flavor, L, mult/J) sector -----------------
function sector_waves(fl::Symbol, L::Int, mult::Int, J::Int; nlevels = 4)
    mQ = mq[String(fl)]
    masses = ConstituentMasses(mQ, mQ)
    r, h = GIModel.radial_grid(NGRID, RMAX)
    V = spin_operator(masses, L, mult, J, r, h)
    # central (spin-independent) FD wave — the W4 audit's choice
    cvals, cvecs, cr = channel_solution(params, masses, L;
        nlevels = max(nlevels, 4), ngrid = NGRID, rmax = RMAX)
    # paper-order: first-order PT in the HO central eigenbasis
    pvals, pwaves, pr = ho_first_order_distorted_states(params_ho, masses, L, V;
        nlevels = max(nlevels, 4), ngrid = NGRID, rmax = RMAX, nbasis = NB)
    # nonperturbative FD resummation (overshoots — kept as the reference)
    H, hr = GIModel.relativistic_hamiltonian(params, masses, L; ngrid = NGRID, rmax = RMAX)
    nvals, nvecs = GIModel.lowest_eigenpairs(Symmetric(Matrix(H) + Matrix(V)), max(nlevels, 4))
    return (central = (cvals, cvecs, cr), paper = (pvals, pwaves, pr),
            nonpert = (nvals, Matrix(nvecs), collect(Float64, hr)))
end

# --- Table VII gluonic rows --------------------------------------------------
# (label, flavor, L, n, mult, J, channel, paper MeV^1/2)
const ROWS = [
    ("eta_c -> 2g",      :c, 0, 1, 1, 0, :S0_2g, +4.700),
    ("psi -> 3g",        :c, 0, 1, 3, 1, :S1_3g, +0.420),
    ("eta'_c -> 2g",     :c, 0, 2, 1, 0, :S0_2g, -2.700),
    ("psi' -> 3g",       :c, 0, 2, 3, 1, :S1_3g, -0.280),
    ("eta_b -> 2g",      :b, 0, 1, 1, 0, :S0_2g, +2.500),
    ("Upsilon -> 3g",    :b, 0, 1, 3, 1, :S1_3g, +0.210),
    ("eta'_b -> 2g",     :b, 0, 2, 1, 0, :S0_2g, -1.700),
    ("Upsilon' -> 3g",   :b, 0, 2, 3, 1, :S1_3g, -0.150),
    ("Upsilon'' -> 3g",  :b, 0, 3, 3, 1, :S1_3g, +0.130),
    ("Upsilon''' -> 3g", :b, 0, 4, 3, 1, :S1_3g, -0.110),
    ("chi_2c -> 2g",     :c, 1, 1, 3, 2, :P2_2g, +0.880),
    ("chi_0c -> 2g",     :c, 1, 1, 3, 0, :P0_2g, +2.500),
    ("chi_2b -> 2g",     :b, 1, 1, 3, 2, :P2_2g, +0.350),
    ("chi_0b -> 2g",     :b, 1, 1, 3, 0, :P0_2g, +0.820),
    ("chi'_2b -> 2g",    :b, 1, 2, 3, 2, :P2_2g, -0.370),
    ("chi'_0b -> 2g",    :b, 1, 2, 3, 0, :P0_2g, -0.820),
]

function gluonic_amp(ch, S, M, mQ)
    gluonic_annihilation_amplitude(ch, S, GIModel.alpha_s_q(M), mQ) * sqrt(1000)
end

median_of(xs) = isempty(xs) ? NaN : sort(xs)[cld(length(xs), 2)]

function run_audit()
    # Part 1: central-wave S_L, FD vs HO (basis-fidelity control)
    part1 = NamedTuple[]
    for (fl, L, nmax) in ((:c, 0, 2), (:c, 1, 1), (:b, 0, 4), (:b, 1, 2))
        mQ = mq[String(fl)]
        masses = ConstituentMasses(mQ, mQ)
        _, fdv, fdr = channel_solution(params, masses, L;
            nlevels = max(nmax, 4), ngrid = NGRID, rmax = RMAX)
        _, hov, hor = channel_solution(params_ho, masses, L;
            nlevels = max(nmax, 4), ngrid = NGRID, rmax = RMAX)
        for n = 1:nmax
            Sfd = smeared_S(fdv[:, n], fdr, mQ, L)
            Sho = smeared_S(hov[:, n], hor, mQ, L)
            push!(part1, (fl = fl, L = L, n = n, Sfd = Sfd, Sho = Sho, ratio = Sho / Sfd))
        end
    end

    # Part 2: gluonic rows under the three treatments
    cache = Dict{Tuple{Symbol,Int,Int,Int},Any}()
    part2 = NamedTuple[]
    for (label, fl, L, n, mult, J, ch, paper) in ROWS
        mQ = mq[String(fl)]
        sec = get!(cache, (fl, L, mult, J)) do
            sector_waves(fl, L, mult, J; nlevels = max(n, 4))
        end
        row = Dict{Symbol,Any}(:label => label, :paper => paper)
        for (key, (vals, vecs, r)) in
            ((:central, sec.central), (:paper_order, sec.paper), (:nonpert, sec.nonpert))
            S = smeared_S(vecs[:, n], r, mQ, L)
            model = gluonic_amp(ch, S, vals[n], mQ)
            row[key] = (M = vals[n], S = S, model = model, ratio = abs(model) / abs(paper),
                        sign_ok = sign(model) == sign(paper))
        end
        push!(part2, (; (k => row[k] for k in (:label, :paper, :central, :paper_order, :nonpert))...))
    end
    return part1, part2
end

function write_report(part1, part2)
    mkpath(dirname(REPORT))
    open(REPORT, "w") do io
        println(io, "# W6 — paper-order (HO first-order) validation of the spin-distorted waves")
        println(io)
        println(io, "Scores the Table VII gluonic subtable under three treatments of the")
        println(io, "spin-dependent operators (smeared contact for S-waves, calibrated")
        println(io, "spin-orbit + tensor for ³P_J), against the central-wave W4 baseline")
        println(io, "(`table_vii_annihilation_em.md`). Zero new parameters: the operators are")
        println(io, "the spectrum-calibrated blocks (`contact_hyperfine_operator`,")
        println(io, "`fine_structure_grid_operator` with k_spin_orbit/k_tensor).")
        println(io)
        println(io, "## 1. Basis fidelity control: central-wave S_L, HO vs FD")
        println(io)
        println(io, "The paper-style HO diagonalization (24 states, paper β convention)")
        println(io, "reproduces the FD smeared wavefunction-at-origin to ≤3%. The 15-20% row")
        println(io, "residuals of the central-wave audit are NOT a basis-fidelity artifact —")
        println(io, "and truncation lowers S_L, the wrong direction to explain them.")
        println(io)
        println(io, "| flavor | L | n | S_FD | S_HO | S_HO/S_FD |")
        println(io, "|:-:|:-:|:-:|---:|---:|:-:|")
        for p in part1
            @printf(io, "| %s | %d | %d | %+.4f | %+.4f | %.4f |\n",
                p.fl, p.L, p.n, p.Sfd, p.Sho, p.ratio)
        end
        println(io)
        println(io, "## 2. Gluonic subtable under the three treatments")
        println(io)
        println(io, "`ratio` = |model|/|paper| per treatment. `central` = spin-independent")
        println(io, "wave (the W4 audit); `paper-order` = first-order PT in the HO central")
        println(io, "eigenbasis (`ho_first_order_distorted_states`); `nonpert` = the operator")
        println(io, "resummed exactly on the FD grid.")
        println(io)
        println(io, "| decay | paper | central | paper-order | nonpert FD | M_paper-order |")
        println(io, "|---|---:|:-:|:-:|:-:|---:|")
        for p in part2
            @printf(io, "| `%s` | %+.3f | %.2f | **%.2f** | %.2f | %.3f |\n",
                p.label, p.paper, p.central.ratio, p.paper_order.ratio,
                p.nonpert.ratio, p.paper_order.M)
        end
        println(io)
        for (key, name) in ((:central, "central"), (:paper_order, "paper-order"), (:nonpert, "nonpert FD"))
            rs = [getfield(p, key).ratio for p in part2]
            @printf(io, "- %s: median %.3f, spread [%.2f, %.2f]\n",
                name, median_of(rs), minimum(rs), maximum(rs))
        end
        println(io)
        println(io, "## 3. Splitting-pattern collapse")
        println(io)
        println(io, "Ratio-of-ratios inside a multiplet sharing one central wave; 1.00 means")
        println(io, "the treatment carries the paper's full spin-splitting of the origin.")
        println(io)
        println(io, "| pattern | central | paper-order |")
        println(io, "|---|:-:|:-:|")
        pairs = [("eta_c/psi", 1, 2), ("eta'_c/psi'", 3, 4), ("eta_b/Upsilon", 5, 6),
                 ("eta'_b/Upsilon'", 7, 8), ("chi_0c/chi_2c", 12, 11),
                 ("chi_0b/chi_2b", 14, 13), ("chi'_0b/chi'_2b", 16, 15)]
        for (name, i, j) in pairs
            @printf(io, "| `%s` | %.2f | %.2f |\n", name,
                part2[i].central.ratio / part2[j].central.ratio,
                part2[i].paper_order.ratio / part2[j].paper_order.ratio)
        end
        println(io)
        println(io, "## Conclusion")
        println(io)
        println(io, "The W4 residual structure (singlet low / triplet high, ³P₀ low / ³P₂")
        println(io, "high) is the paper's spin-dependent wavefunction distortion, treated at")
        println(io, "first order in the HO eigenbasis. The nonperturbative resummation")
        println(io, "overshoots — evidence that the paper's Table VII waves carry the spin")
        println(io, "blocks at paper order, consistent with the calibrated k_spin_orbit/")
        println(io, "k_tensor bridge used for the mass spectrum. The former attribution of")
        println(io, "the charm rows to FD-vs-HO wavefunction-at-origin infidelity is refuted")
        println(io, "by part 1.")
    end
end

part1, part2 = run_audit()
write_report(part1, part2)

println("Part 1 — central S_L HO/FD: worst |1-ratio| = ",
    @sprintf("%.4f", maximum(abs(1 - p.ratio) for p in part1)))
for (key, name) in ((:central, "central"), (:paper_order, "paper-order"), (:nonpert, "nonpert"))
    rs = [getfield(p, key).ratio for p in part2]
    @printf("Part 2 — %-11s median %.3f  spread [%.2f, %.2f]\n",
        name, median_of(rs), minimum(rs), maximum(rs))
end
all(p -> p.paper_order.sign_ok, part2) || error("paper-order sign mismatch")
println("report: ", REPORT)
