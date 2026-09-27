#!/usr/bin/env julia
# Like-for-like P-wave multiplet parameters (M, S, T, L, V) from GI's printed
# masses and caption angles and from GIModel's final masses and angles, using
# the same inversion of GI Eqs. (23)-(26):
#   3P0 = M + S/4 - 2T - 2L,  3P2 = M + S/4 - T/5 + L,
#   J=1 block [[M - 3S/4, V], [V, M + S/4 + T - L]]  (1P1, 3P1 basis),
# with the block rebuilt from its two eigenvalues and Q_low = cos θ 1P1 + sin θ 3P1.
# GI rounding: masses uniform ±5 MeV, angles ±0.5°, Monte Carlo 68% intervals.
#   julia --project=GIPaper/scripts GIPaper/scripts/investigate_multiplet_parameters.jl
using GIModel, CSV, DataFrames, Printf, Random, Statistics
using GIPaper: quark_for

params, mq = load_parameters_and_quark_masses(default_parameters_path())
root = dirname(@__DIR__)
out = joinpath(root, "docs", "investigations", "mixing_composition")
mkpath(out)

# (label, f1, f2, n, GI masses MeV: 3P0, J=1 low, J=1 high, 3P2; GI θ or nothing for equal mass)
# Equal mass: "low/high" are the printed 1P1 and 3P1 masses in that order.
const SETS = [
    ("uu 1P", :u, :u, 1, (1090, 1220, 1240, 1310), nothing),
    ("uu 2P", :u, :u, 2, (1780, 1780, 1820, 1820), nothing),
    ("ss 1P", :s, :s, 1, (1360, 1470, 1480, 1530), nothing),
    ("cc 1P", :c, :c, 1, (3440, 3520, 3510, 3550), nothing),
    ("cc 2P", :c, :c, 2, (3920, 3960, 3950, 3980), nothing),
    ("bb 1P", :b, :b, 1, (9850, 9880, 9880, 9900), nothing),
    ("bb 2P", :b, :b, 2, (10230, 10250, 10250, 10260), nothing),
    ("us 1P", :u, :s, 1, (1240, 1340, 1380, 1430), 34.0),
    ("us 2P", :u, :s, 2, (1890, 1900, 1930, 1940), 15.0),
    ("cu 1P", :c, :u, 1, (2400, 2440, 2490, 2500), -41.0),
    ("cs 1P", :c, :s, 1, (2480, 2530, 2570, 2590), -44.0),
]

function invert(E0, lo, hi, E2, θ)
    if isnothing(θ)                       # equal mass: lo = 1P1, hi = 3P1, V = 0
        Es, Et, V = lo, hi, 0.0
    else
        c, s = cosd(θ), sind(θ)
        Es, Et, V = c^2 * lo + s^2 * hi, s^2 * lo + c^2 * hi, c * s * (lo - hi)
    end
    T = (3(Et - E0) - (E2 - E0)) / 7.2
    L = (Et - E0) - 3T
    S = (E2 + T / 5 - L) - Es
    (M = Es + 3S / 4, S = S, T = T, L = L, V = V, gap = Es - Et)
end

function model_masses(f1, f2, n, θGI)
    levels = [BasisState(n, "P", 3, 0), BasisState(n, "P", 1, 1), BasisState(n, "P", 3, 1), BasisState(n, "P", 3, 2)]
    spec = compute_spectrum(params, Meson(quark_for(mq, f1), quark_for(mq, f2)); levels)
    E(m, J) = 1000spectrum_state(spec, n, "P", m, J).mass_GeV
    isnothing(θGI) && return (E(3, 0), E(1, 1), E(3, 1), E(3, 2), nothing)
    s1, s3 = spectrum_state(spec, n, "P", 1, 1), spectrum_state(spec, n, "P", 3, 1)
    m1 = only(filter(m -> m.mechanism == "antisymmetric_spin_orbit", s1.mixings))
    m3 = only(filter(m -> m.mechanism == "antisymmetric_spin_orbit", s3.mixings))
    b = m1.result.block.basis
    is = findfirst(x -> x.label == s1.label, b); it = findfirst(x -> x.label == s3.label, b)
    phys = s1.mass_GeV <= s3.mass_GeV ? m1 : m3
    cs, ct = phys.components[is], phys.components[it]
    ph = cs < 0 ? -1 : 1
    (E(3, 0), 1000min(s1.mass_GeV, s3.mass_GeV), 1000max(s1.mass_GeV, s3.mass_GeV), E(3, 2), atand(ph * ct, ph * cs))
end

rng = MersenneTwister(1985)
rows = NamedTuple[]
for (label, f1, f2, n, (E0, lo, hi, E2), θ) in SETS
    gi = invert(E0, lo, hi, E2, θ)
    draws = [invert(E0 + 10rand(rng) - 5, lo + 10rand(rng) - 5, hi + 10rand(rng) - 5, E2 + 10rand(rng) - 5,
                    isnothing(θ) ? nothing : θ + rand(rng) - 0.5) for _ in 1:20000]
    q(f) = (x = sort([getfield(d, f) for d in draws]); (x[3200], x[16800]))
    me = invert(model_masses(f1, f2, n, θ)...)
    for f in (:S, :T, :L, :V, :gap)
        lo68, hi68 = q(f)
        push!(rows, (set = label, parameter = String(f), GI = getfield(gi, f), GI_lo68 = lo68, GI_hi68 = hi68,
            model = getfield(me, f), inside = lo68 <= getfield(me, f) <= hi68))
    end
end
df = DataFrame(rows)
CSV.write(joinpath(out, "multiplet_parameters.csv"), df)
@printf("%-6s %-4s %8s %17s %8s  %s\n", "set", "par", "GI", "GI 68% (rounding)", "GIModel", "")
for r in eachrow(df)
    @printf("%-6s %-4s %+8.1f  [%+6.1f, %+6.1f]  %+8.1f  %s\n", r.set, r.parameter, r.GI, r.GI_lo68, r.GI_hi68, r.model, r.inside ? "" : "<-- outside")
end
