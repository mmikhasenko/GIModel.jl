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
#  3. The paper's treatment is FULL DIAGONALIZATION in the FINITE HO basis, not
#     first-order PT. The light ¹S₀ (pion) mass discriminates (§1b): first-order
#     PT over-raises it to ~0.28 GeV, while the finite-HO full diagonalization
#     keeps it light (~0.10 GeV) like the fine-grid FD resummation. Full
#     diagonalization in the finite paper-β basis lands all 16 gluonic rows in
#     [0.92, 1.14] (median 1.02); the fine-grid FD *over*-resums (eta_b 1.21 vs
#     the finite-basis 1.14) — a grid-resolution effect, not perturbation order.
#     First-order PT coincides only in the heavy-quark (small-V) limit.
#
# Wave sources: contact_hyperfine_operator (S-waves, mult 1/3),
# fine_structure_grid_operator (3P_J), ho_full_distorted_states (paper-order
# treatment), ho_first_order_distorted_states (heavy-quark PT proxy),
# H_central + V on the fine grid (nonperturbative FD reference).

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
solver_ho = OscillatorSolver(nbasis = NB, ngrid = NGRID, rmax = RMAX)

# Phase convention: outermost antinode positive (shared with audit_table_vii.jl).
function fix_outer_antinode_positive!(u)
    peak = maximum(abs, u)
    i = findlast(x -> abs(x) > 0.2 * peak, u)
    (i !== nothing && u[i] < 0) && (u .*= -1)
    return u
end

function smeared_S(wave::RadialWave, mQ, L)
    r, _ = GIModel.radial_grid(NGRID, RMAX)
    sampled = wave isa MeshWave ? wave : sample_wave(wave, r)
    fixed = MeshWave(fix_outer_antinode_positive!(copy(sampled.u)), sampled.r)
    return wavefunction_origin_smearing(fixed, mQ; L = L, npoints = NPTS)
end

# --- spin-dependent grid operator for a sector ------------------------------
# L=0: smeared contact (multiplicity 1 or 3); L=1: triplet 3P_J LS+tensor.
function spin_operator(masses, L, mult, J, r, h)
    L == 0 && return GIModel.contact_hyperfine_operator(params, masses, "S", mult, r)
    return fine_structure_grid_operator(params, masses, J, r, h; L = L)
end

# --- the four treatments of one (flavor, L, mult/J) sector ------------------
function sector_waves(fl::Symbol, L::Int, mult::Int, J::Int; nlevels = 4)
    mQ = mq[String(fl)]
    masses = ConstituentMasses(mQ, mQ)
    r, h = GIModel.radial_grid(NGRID, RMAX)
    V = spin_operator(masses, L, mult, J, r, h)
    # central (spin-independent) FD wave — the pre-harmonization gluonic choice
    central = channel_solution(
        params, masses, L;
        nlevels = max(nlevels, 4),
        solver = FiniteDifferenceSolver(ngrid = NGRID, rmax = RMAX),
    )
    # first-order PT in the HO central eigenbasis (heavy-quark proxy)
    first_order = ho_first_order_distorted_states(params, masses, L, V;
        solver = solver_ho, nlevels = max(nlevels, 4))
    # paper-order: FULL diagonalization of H_central + V in the finite HO basis
    paper = ho_full_distorted_states(params, masses, L, V;
        solver = solver_ho, nlevels = max(nlevels, 4))
    # nonperturbative FD resummation on the fine grid (over-resums)
    H, hr = GIModel.relativistic_hamiltonian(params, masses, L; ngrid = NGRID, rmax = RMAX)
    nvals, nvecs = GIModel.lowest_eigenpairs(Symmetric(Matrix(H) + Matrix(V)), max(nlevels, 4))
    nonpert = ChannelRadialSolution(nvals, Matrix(nvecs), collect(Float64, hr))
    return (; central, first_order, paper, nonpert)
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
        fd_solution = channel_solution(
            params, masses, L;
            nlevels = max(nmax, 4),
            solver = FiniteDifferenceSolver(ngrid = NGRID, rmax = RMAX),
        )
        ho_solution = channel_solution(params, masses, L;
            solver = solver_ho, nlevels = max(nmax, 4))
        for n = 1:nmax
            Sfd = smeared_S(radial_wave(fd_solution, n), mQ, L)
            Sho = smeared_S(radial_wave(ho_solution, n), mQ, L)
            push!(part1, (fl = fl, L = L, n = n, Sfd = Sfd, Sho = Sho, ratio = Sho / Sfd))
        end
    end

    # Part 1b: light ¹S₀ (pion) mass discriminator — first-order PT vs the
    # finite-HO full diagonalization vs the fine-grid FD resummation. This is
    # what rules out first-order PT as the paper's actual mechanism.
    mn = ConstituentMasses(mq["q"], mq["q"])
    r0, _ = GIModel.radial_grid(NGRID, RMAX)
    Vpi = GIModel.contact_hyperfine_operator(params, mn, "S", 1, r0)
    fo = ho_first_order_distorted_states(params, mn, 0, Vpi;
        solver = solver_ho, nlevels = 4)
    fu = ho_full_distorted_states(params, mn, 0, Vpi;
        solver = solver_ho, nlevels = 4)
    Hn, _ = GIModel.relativistic_hamiltonian(params, mn, 0; ngrid = NGRID, rmax = RMAX)
    nv, _ = GIModel.lowest_eigenpairs(Symmetric(Matrix(Hn) + Matrix(Vpi)), 4)
    pion = (first_order = fo.eigenvalues_GeV[1], full = fu.eigenvalues_GeV[1], nonpert = nv[1])

    # Part 2: gluonic rows under the four treatments
    cache = Dict{Tuple{Symbol,Int,Int,Int},Any}()
    part2 = NamedTuple[]
    for (label, fl, L, n, mult, J, ch, paper) in ROWS
        mQ = mq[String(fl)]
        sec = get!(cache, (fl, L, mult, J)) do
            sector_waves(fl, L, mult, J; nlevels = max(n, 4))
        end
        row = Dict{Symbol,Any}(:label => label, :paper => paper)
        for (key, solution) in
            ((:central, sec.central), (:first_order, sec.first_order),
             (:full, sec.paper), (:nonpert, sec.nonpert))
            M = solution.eigenvalues_GeV[n]
            S = smeared_S(radial_wave(solution, n), mQ, L)
            model = gluonic_amp(ch, S, M, mQ)
            row[key] = (M = M, S = S, model = model, ratio = abs(model) / abs(paper),
                        sign_ok = sign(model) == sign(paper))
        end
        push!(part2, (; (k => row[k] for k in
            (:label, :paper, :central, :first_order, :full, :nonpert))...))
    end
    return part1, pion, part2
end

function write_report(part1, pion, part2)
    mkpath(dirname(REPORT))
    open(REPORT, "w") do io
        println(io, "# W6 — paper-order (finite HO-basis) validation of the spin-distorted waves")
        println(io)
        println(io, "Scores the Table VII gluonic subtable under four treatments of the")
        println(io, "spin-dependent operators (smeared contact for S-waves, calibrated")
        println(io, "spin-orbit + tensor for ³P_J), against the spin-independent central-wave")
        println(io, "baseline. Zero new parameters: the operators are the spectrum-calibrated")
        println(io, "blocks (`contact_hyperfine_operator`, `fine_structure_grid_operator` with")
        println(io, "k_spin_orbit/k_tensor). The **paper-order** treatment — full diagonalization")
        println(io, "of `H_central + V_spin` in the finite paper-β HO basis")
        println(io, "(`ho_full_distorted_states`) — is what the harmonized Table VII audit uses")
        println(io, "for every subtable.")
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
        println(io, "## 1b. What order is \"paper order\"? The light ¹S₀ mass discriminator")
        println(io)
        println(io, "The residuals are spin-dependent wavefunction distortion (§2). But *how*")
        println(io, "the paper carries the distortion is settled by the light `¹S₀` nonstrange")
        println(io, "mass, where the contact term is enormous. First-order PT and full")
        println(io, "diagonalization agree for heavy spin splittings but diverge sharply here:")
        println(io)
        @printf(io, "| treatment | pion ¹S₀ mass (GeV) |\n")
        println(io, "|---|---:|")
        @printf(io, "| first-order PT (`ho_first_order_distorted_states`) | %.4f |\n", pion.first_order)
        @printf(io, "| **finite-HO full diag (`ho_full_distorted_states`)** | **%.4f** |\n", pion.full)
        @printf(io, "| nonperturbative FD (fine grid) | %.4f |\n", pion.nonpert)
        println(io)
        println(io, "First-order PT over-raises the pion to ≈0.28 GeV; the paper keeps it light")
        println(io, "(≈0.10 GeV), which BOTH the finite-HO full diagonalization and the fine-grid")
        println(io, "FD resummation reproduce. So the paper's mechanism is **full diagonalization**")
        println(io, "(the contact is resummed, not truncated at first order); first-order PT is")
        println(io, "only a heavy-quark proxy. The finite HO basis — not a perturbation order —")
        println(io, "is why the heavy gluonic ratios sit *below* the fine-grid FD overshoot (§2):")
        println(io, "a finite basis resums the contact less aggressively than the fine FD grid.")
        println(io)
        println(io, "## 2. Gluonic subtable under the four treatments")
        println(io)
        println(io, "`ratio` = |model|/|paper| per treatment. `central` = spin-independent wave;")
        println(io, "`1st-PT` = first-order PT in the HO central eigenbasis; **`paper` = full")
        println(io, "diagonalization in the finite HO basis** (the harmonized-audit treatment);")
        println(io, "`nonpert` = the operator resummed on the fine FD grid.")
        println(io)
        println(io, "| decay | paper amp | central | 1st-PT | **paper** | nonpert FD | M_paper |")
        println(io, "|---|---:|:-:|:-:|:-:|:-:|---:|")
        for p in part2
            @printf(io, "| `%s` | %+.3f | %.2f | %.2f | **%.2f** | %.2f | %.3f |\n",
                p.label, p.paper, p.central.ratio, p.first_order.ratio,
                p.full.ratio, p.nonpert.ratio, p.full.M)
        end
        println(io)
        for (key, name) in ((:central, "central"), (:first_order, "1st-order PT"),
                            (:full, "paper (full diag)"), (:nonpert, "nonpert FD"))
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
        println(io, "| pattern | central | paper (full diag) |")
        println(io, "|---|:-:|:-:|")
        pairs = [("eta_c/psi", 1, 2), ("eta'_c/psi'", 3, 4), ("eta_b/Upsilon", 5, 6),
                 ("eta'_b/Upsilon'", 7, 8), ("chi_0c/chi_2c", 12, 11),
                 ("chi_0b/chi_2b", 14, 13), ("chi'_0b/chi'_2b", 16, 15)]
        for (name, i, j) in pairs
            @printf(io, "| `%s` | %.2f | %.2f |\n", name,
                part2[i].central.ratio / part2[j].central.ratio,
                part2[i].full.ratio / part2[j].full.ratio)
        end
        println(io)
        println(io, "## Conclusion")
        println(io)
        println(io, "The W4 residual structure (singlet low / triplet high, ³P₀ low / ³P₂ high)")
        println(io, "is the paper's spin-dependent wavefunction distortion, carried by **full")
        println(io, "diagonalization of `H_central + V_spin` in the finite paper-β HO basis** —")
        println(io, "the paper's literal method, now used across the whole harmonized Table VII")
        println(io, "audit. First-order PT reproduces the heavy gluonic ratios (its perturbative")
        println(io, "limit) but is refuted as the mechanism by the light `¹S₀` mass (§1b): it")
        println(io, "over-raises the pion, whereas full diagonalization keeps it light like the")
        println(io, "FD resummation. The fine-grid FD *over*-resums (η_b 1.21 vs the finite-basis")
        println(io, "1.14) — a grid-resolution effect, not perturbation order. The former")
        println(io, "attribution of the charm rows to FD-vs-HO wavefunction-at-origin infidelity")
        println(io, "is refuted by §1; the former attribution to a first-order treatment is")
        println(io, "refuted by §1b.")
    end
end

part1, pion, part2 = run_audit()
write_report(part1, pion, part2)

println("Part 1 — central S_L HO/FD: worst |1-ratio| = ",
    @sprintf("%.4f", maximum(abs(1 - p.ratio) for p in part1)))
@printf("Part 1b — pion ¹S₀ mass: 1st-PT %.4f | full-diag %.4f | nonpert FD %.4f\n",
    pion.first_order, pion.full, pion.nonpert)
for (key, name) in ((:central, "central"), (:first_order, "1st-PT"),
                    (:full, "paper full"), (:nonpert, "nonpert"))
    rs = [getfield(p, key).ratio for p in part2]
    @printf("Part 2 — %-11s median %.3f  spread [%.2f, %.2f]\n",
        name, median_of(rs), minimum(rs), maximum(rs))
end
all(p -> p.full.sign_ok, part2) || error("paper full-diag sign mismatch")
println("report: ", REPORT)
