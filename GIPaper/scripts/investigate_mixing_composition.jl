#!/usr/bin/env julia
# Why did GIModel not reproduce the GI (1985) mixing compositions?
#
#   julia --project=GIPaper/scripts GIPaper/scripts/investigate_mixing_composition.jl --label after
#   julia --project=GIPaper/scripts GIPaper/scripts/investigate_mixing_composition.jl --report
#
# `--label L` computes every target/diagnostic with the GIModel it is loaded
# against and writes <out>/*_L.csv. The "before" ledger is this same script run
# in a checkout of commit 88e9a99 (S-wave-only contact); `--out DIR` redirects
# the CSVs. `--report` reads both labels and writes the results note + figure.
# Narrative, hypotheses and conclusion: docs/investigations/mixing_composition_investigation.md.
using GIModel, GIPaper, CSV, DataFrames, LinearAlgebra, Printf, Statistics
using GIPaper: quark_for
const G = GIModel

root = dirname(@__DIR__)
arg(name, default) = (i = findfirst(==(name), ARGS); isnothing(i) ? default : ARGS[i+1])
out = arg("--out", joinpath(root, "docs", "investigations", "mixing_composition"))
mkpath(out)
params, mq = load_parameters_and_quark_masses(default_parameters_path())
const LOF = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)

# Source-verified targets (page images, see the investigation note).
# (f1, f2, n, L, θ_GI, printed lower/upper J=L masses in MeV or missing)
const ANGLES = [
    (:u,:s,1,"P", 34.0, 1340, 1380), (:u,:s,1,"D", 33.0, 1780, 1810), (:u,:s,2,"P", 15.0, 1900, 1930),
    (:u,:s,1,"F", 32.0, 2120, 2150), (:u,:s,2,"D", 25.0, 2230, 2260), (:u,:s,1,"G", 33.0, 2410, 2440),
    (:c,:u,1,"P",-41.0, 2440, 2490), (:c,:u,1,"D",-39.0, missing, missing),
    (:c,:s,1,"P",-44.0, 2530, 2570), (:c,:s,1,"D",-39.0, missing, missing),
    (:b,:u,1,"P",-43.0, missing, missing), (:b,:s,1,"P",-45.0, missing, missing),
    (:b,:c,1,"P",-53.0, missing, missing)]
# GI amplitudes; dominant component first. Fig. 3 is the light isovector.
const TENSOR = [
    ("Fig. 3 rho(1.45)", :u, :u, (2, "S"), [(2,"S") => 1.00, (1,"D") => 0.04]),
    ("Fig. 4 K*(1.58)",  :u, :s, (2, "S"), [(2,"S") => 1.00, (1,"D") => 0.04]),
    ("Fig. 6 psi(3.82)", :c, :c, (1, "D"), [(1,"D") => 1.00, (1,"S") => 0.01, (2,"S") => -0.03, (3,"S") => -0.01])]

meson(f1, f2) = Meson(quark_for(mq, f1), quark_for(mq, f2))
same_j(s) = only(filter(m -> m.mechanism == "antisymmetric_spin_orbit", s.mixings))
pcontact(m, L, w) = (32π / (9m.m1_GeV * m.m2_GeV)) * G.radial_expect_momentum_sandwich(
    params, m, L, w, params.factors.epsilon_c, (r, _) -> G.smeared_contact_kernel(params, m, r))

# Lower physical state of an nL pair projected on its same-n singlet/triplet rows.
function pair_angle(spec, n, L)
    J = LOF[L]
    s1, s3 = spectrum_state(spec, n, L, 1, J), spectrum_state(spec, n, L, 3, J)
    m1, m3 = same_j(s1), same_j(s3)
    basis = m1.result.block.basis
    is = findfirst(b -> b.label == s1.label, basis); it = findfirst(b -> b.label == s3.label, basis)
    phys = s1.mass_GeV <= s3.mass_GeV ? m1 : m3
    cs, ct = phys.components[is], phys.components[it]
    ph = cs < 0 ? -1 : 1
    H = m1.result.block.matrix
    (theta = atand(ph * ct, ph * cs), singlet = cs^2 / (cs^2 + ct^2), outside = max(0.0, 1 - cs^2 - ct^2),
     gap = 1000 * (H[is, is] - H[it, it]), V = 1000 * H[is, it],
     low = 1000min(s1.mass_GeV, s3.mass_GeV), high = 1000max(s1.mass_GeV, s3.mass_GeV), s1, s3)
end

function angle_rows(solver)
    rows = NamedTuple[]
    for (f1, f2) in unique((a[1], a[2]) for a in ANGLES)
        mine = filter(a -> (a[1], a[2]) == (f1, f2), ANGLES)
        levels = vcat([[BasisState(a[3], a[4], 1, LOF[a[4]]), BasisState(a[3], a[4], 3, LOF[a[4]])] for a in mine]...)
        me = meson(f1, f2); m = me.constituent_masses
        spec = compute_spectrum(params, me; levels, solver)
        for (_, _, n, L, gi, lo, hi) in mine
            p = pair_angle(spec, n, L)
            ws, wt = radial_wave(spec, p.s1.corrected), radial_wave(spec, p.s3.corrected)
            op = spin_orbit_mixing_components(params, m, L, ws, wt)
            push!(rows, (sector = "$f1$f2", nL = "$n$L", theta = p.theta, theta_GI = gi,
                singlet = p.singlet, singlet_GI = cosd(gi)^2, outside = p.outside,
                distance = abs(mod(p.theta - gi + 90, 180) - 90),
                gap = p.gap, V = p.V, V_vector = 1000op.vector, V_thomas = 1000op.thomas,
                contact_S_singlet = 1000pcontact(m, LOF[L], ws), contact_S_triplet = 1000pcontact(m, LOF[L], wt),
                low = p.low, high = p.high, low_GI = lo, high_GI = hi,
                # conditional two-state inversion of GI's printed masses and angle
                gap_GI = ismissing(lo) ? missing : -(hi - lo) * cosd(2gi),
                V_GI = ismissing(lo) ? missing : -(hi - lo) * sind(2gi) / 2))
        end
    end
    rows
end

function tensor_rows(; nS = 5, nD = 3)
    rows = NamedTuple[]
    for (tag, f1, f2, (n0, L0), target) in TENSOR
        levels = vcat([BasisState(n, "S", 3, 1) for n in 1:nS], [BasisState(n, "D", 3, 1) for n in 1:nD])
        spec = compute_spectrum(params, meson(f1, f2); levels)
        st = spectrum_state(spec, n0, L0, 3, 1)
        mix = only(filter(x -> x.mechanism == "tensor_mixing", st.mixings))
        basis = mix.result.block.basis; c = collect(mix.components)
        dom = findfirst(b -> (b.n, b.L_label) == (n0, L0), basis)
        c ./= c[dom]
        for (b, x) in zip(basis, c)
            gi = something(findfirst(t -> t.first == (b.n, b.L_label), target), 0)
            push!(rows, (state = tag, mass = 1000st.mass_GeV, component = b.label,
                amplitude_outer = x,                        # GIModel phase: outer lobe positive
                amplitude_origin = x * (-1)^(b.n - 1) * (-1)^(n0 - 1),  # wave positive at origin
                GI = gi == 0 ? missing : target[gi].second))
        end
    end
    rows
end

# Spectrum-wide test independent of mixing: every catalog, split by L and 2S+1.
function global_rows()
    rows = NamedTuple[]
    dd = GIPaper.paper_data_dir()
    for fn in readdir(dd)
        startswith(fn, "reference_spectrum_") || continue
        base = fn[length("reference_spectrum_")+1:end-4]
        ref = GIPaper.load_reference_spectrum(joinpath(dd, fn))
        for r in compare_reference(params, mq, ref; solver = FiniteDifferenceSolver(), contact_hyperfine = true,
                use_fine_structure = true, isoscalar_pseudoscalar_annihilation = base == "isoscalar" ? :table_iii : :none,
                strange_mass_GeV = base == "isoscalar" ? mq["s"] : nothing)
            push!(rows, (catalog = base, state = r.state, n = r.n, L = r.L, multiplicity = r.multiplicity,
                J = r.J, reference = 1000r.reference_GeV, predicted = 1000r.predicted_GeV, residual = r.residual_MeV))
        end
    end
    rows
end

# Stage 2 in the limit: full coupled-channel grid problem [H(1L_L) W; W H(3L_L)]
# versus second-stage blocks truncated to N radial eigenstates per sector.
function antisym_grid(m, L, r)
    h = r[2] - r[1]; m1, m2 = m.m1_GeV, m.m2_GeV; f = params.factors
    p2 = eigen(G.p2_operator(params, m1, L, r, h))
    sw(pair, eps, k) = (B = G.momentum_relativization_matrix(pair.m1_GeV, pair.m2_GeV,
        G.gi_spin_dependent_side_exponent(eps), p2); B * Diagonal(k.(r)) * B)
    p11, p22 = G._mass_pair(m1, m1), G._mass_pair(m2, m2)
    term(p, mi) = (sw(p, f.epsilon_so_vector, x -> G._vector_so_kernel(params, p, x)) -
                   sw(p, f.epsilon_so_scalar, x -> G._scalar_so_kernel(params, p, x))) / (4mi^2)
    sqrt(L * (L + 1.0)) * (term(p11, m1) - term(p22, m2))
end

function stage2_rows(solver = FiniteDifferenceSolver())
    rows = NamedTuple[]
    for (f1, f2, n, Ls, gi, _, _) in ANGLES
        L = LOF[Ls]; m = meson(f1, f2).constituent_masses
        mats(mult) = G._fixed_channel_matrices(solver, params, m, FineStructureMultiplet(Ls, mult, L), nothing, SpinTerms(), 1)
        Ms, Mt = mats(1), mats(3); r = Ms.r
        es, Us = eigen(Symmetric(Matrix(Ms.total))); et, Ut = eigen(Symmetric(Matrix(Mt.total)))
        # production phase: outer lobe positive
        for U in (Us, Ut), k in axes(U, 2)
            i = findlast(x -> abs(x) > 0.2maximum(abs, U[:, k]), U[:, k]); U[i, k] < 0 && (U[:, k] .*= -1)
        end
        W = antisym_grid(m, L, r); Wn = Us' * W * Ut
        low(cs, ct) = (ph = cs < 0 ? -1 : 1; atand(ph * ct, ph * cs))
        for N in (n, n + 1, n + 3, 10, 40)
            E, U = eigen(Symmetric([Diagonal(es[1:N]) Wn[1:N, 1:N]; Wn[1:N, 1:N]' Diagonal(et[1:N])]))
            k = minimum(sort(sortperm([U[n, k]^2 + U[N+n, k]^2 for k in axes(U, 2)], rev = true)[1:2]))
            push!(rows, (sector = "$f1$f2", nL = "$n$Ls", basis = "N=$N", theta = low(U[n, k], U[N+n, k]),
                outside = 1 - U[n, k]^2 - U[N+n, k]^2))
        end
        ng = length(r)
        Ef, Uf = eigen(Symmetric([Matrix(Ms.total) W; W Matrix(Mt.total)]))
        proj(k) = (dot(Us[:, n], Uf[1:ng, k]), dot(Ut[:, n], Uf[ng+1:end, k]))
        k = minimum(sort(sortperm([sum(abs2, proj(k)) for k in axes(Uf, 2)], rev = true)[1:2]))
        cs, ct = proj(k)
        push!(rows, (sector = "$f1$f2", nL = "$n$Ls", basis = "coupled grid", theta = low(cs, ct), outside = 1 - cs^2 - ct^2))
    end
    rows
end

# Controlled substitutions in the antisymmetric element only (diagnostics, not fixes).
function variant_rows()
    rows = NamedTuple[]; f = params.factors
    for (f1, f2, n, Ls, gi, lo, hi) in filter(a -> !ismissing(a[6]), ANGLES)
        L = LOF[Ls]; me = meson(f1, f2); m = me.constituent_masses
        spec = fixed_spectrum(params, me; levels = [BasisState(n, Ls, 1, L), BasisState(n, Ls, 3, L)])
        wl = radial_wave(spec, spectrum_state(spec, n, Ls, 1, L)); wr = radial_wave(spec, spectrum_state(spec, n, Ls, 3, L))
        function element(; smear = :ii, factor = :ii, vf = true, tf = true)
            total = 0.0
            for (i, sgn, mi) in ((1, 1.0, m.m1_GeV), (2, -1.0, m.m2_GeV))
                pii = G._mass_pair(mi, mi)
                sp = smear == :ii ? pii : m; fp = factor == :ii ? pii : m
                x(eps, k, on) = on ? G.radial_cross_expect_momentum_sandwich(params, fp, L, wl, L, wr, eps, (r, _) -> k(r)) :
                                     G.radial_overlap(wl, wr, k)
                total += sgn * sqrt(L * (L + 1.0)) * (x(f.epsilon_so_vector, r -> G._vector_so_kernel(params, sp, r), vf) -
                                                     x(f.epsilon_so_scalar, r -> G._scalar_so_kernel(params, sp, r), tf)) / (4mi^2)
            end
            1000total
        end
        VGI = -(hi - lo) * sind(2gi) / 2
        for (name, kw) in (("production", (;)), ("sigma_12 smearing", (smear = :pair,)), ("f_12 momentum factor", (factor = :pair,)),
                           ("no momentum factor", (vf = false, tf = false)), ("no factor on Thomas", (tf = false,)),
                           ("no factor on vector", (vf = false,)))
            v = element(; kw...)
            push!(rows, (sector = "$f1$f2", nL = "$n$Ls", variant = name, V = v, V_GI = VGI, ratio = VGI / v))
        end
    end
    rows
end

label = arg("--label", nothing)
if !isnothing(label)
    t0 = time()
    CSV.write(joinpath(out, "angles_$label.csv"), angle_rows(FiniteDifferenceSolver()))
    CSV.write(joinpath(out, "tensor_$label.csv"), tensor_rows())
    CSV.write(joinpath(out, "global_$label.csv"), global_rows())
    if label == "after"
        # Stage 1: independent grids/domain and the native HO solver.
        for (name, solver) in (("fd900", FiniteDifferenceSolver(ngrid = 900)), ("fd600_rmax32", FiniteDifferenceSolver(ngrid = 600, rmax = 32.0)),
                               ("ho", OscillatorSolver()), ("ho_nbasis6_unconverged", OscillatorSolver(nbasis = 6, converge = false)))
            CSV.write(joinpath(out, "angles_after_$name.csv"), angle_rows(solver))
        end
        CSV.write(joinpath(out, "stage2_after.csv"), stage2_rows())
        CSV.write(joinpath(out, "variants_after.csv"), variant_rows())
    end
    @printf("wrote %s CSVs to %s in %.0f s\n", label, out, time() - t0)
end

"--report" in ARGS && include(joinpath(@__DIR__, "report_mixing_composition.jl"))
