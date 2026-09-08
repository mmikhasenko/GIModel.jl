#!/usr/bin/env julia
# =============================================================================
# Density clouds: the qqbar separation density |psi(r)|^2 as sampled points
# =============================================================================
# One dot per Monte-Carlo draw from |psi(r)|^2 -- the "picture of an atom" that
# shows the density itself rather than an isosurface hiding everything inside
# it. What is being drawn is the *relative-coordinate* density of the quark and
# the antiquark: r is their separation, not a position in the lab, and the plot
# is not a charge distribution.
#
# The pieces:
#   (1) psi(r) = sum_c c_c (u_c(r)/r) |(L_c S_c) J mJ>   from `physical_components`
#   (2) sum over the unobserved quark spin m_s, at fixed mJ. Summing over mJ
#       instead would give sum_m |Y_Lm|^2 = (2L+1)/4pi -- a sphere for every
#       state in the spectrum, which is why mJ is an argument here.
#   (3) exact radial marginal p(r) = sum_c c_c^2 u_c(r)^2 (orthonormality of the
#       Y_Lm kills every cross term), so r is drawn by CDF inversion and only
#       cos(theta) needs rejection. phi is uniform: all L in one spin group
#       share m_L = mJ - m_s, hence one common exp(i m_L phi).
#
# Interference between different L at the same S -- the 3S1/3D1 tensor mixing --
# survives step (2) and is visible as a mJ-dependent distortion of the cloud.
#
# The angular machinery is checked against the radial marginal at run time
# (`check_angular_normalization`): 2 pi int f(x) dx must reproduce
# sum_c c_c^2 u_c(r)^2 at every radius, which is a joint test of the
# Clebsch-Gordan and spherical-harmonic conventions used here.
#
# Run:  julia examples/density_3d.jl
# Writes: figures/density_cloud_panels.png, figures/density_cloud_2S.png
# =============================================================================

using Pkg
Pkg.activate(@__DIR__; io = devnull)     # examples/Project.toml — GIModel + CairoMakie
Pkg.instantiate(; io = devnull)

using Printf
using Random
using GIModel
using CairoMakie

CairoMakie.activate!(type = "png", px_per_unit = 2)

const PARAMS_PATH = joinpath(dirname(@__DIR__), "data", "parameters.provisional.toml")
const HBARC_FM = 0.19733                 # r[fm] = HBARC_FM * r[GeV^-1]
const RMAX_INV_GEV = 30.0
const NRADIAL = 1201

# -----------------------------------------------------------------------------
# Angular algebra (kept explicit; both conventions are verified at run time)
# -----------------------------------------------------------------------------

"""
    legendre_Y(l, m, x)

The theta-dependent factor of `Y_lm`: `sqrt((2l+1)/4pi (l-|m|)!/(l+|m|)!) P_l^|m|(x)`
at `x = cos(theta)`, Condon-Shortley phase included. The `exp(i m phi)` is
dropped because every term that is summed coherently here carries the same `m`.
"""
function legendre_Y(l::Integer, m::Integer, x::Real)
    m = abs(m)
    m > l && return 0.0
    pmm = 1.0
    if m > 0
        somx2 = sqrt(max(0.0, 1 - x^2))
        fact = 1.0
        for _ = 1:m
            pmm *= -fact * somx2
            fact += 2.0
        end
    end
    p = pmm
    if l > m
        pmmp1 = x * (2m + 1) * pmm
        p = pmmp1
        for ll = (m+2):l
            p = ((2ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
            pmm = pmmp1
            pmmp1 = p
        end
    end
    norm = sqrt((2l + 1) / (4pi) * factorial(l - m) / factorial(l + m))
    return norm * p
end

"""
    clebsch_gordan(j1, m1, j2, m2, J, M)

`<j1 m1 j2 m2 | J M>` by the Racah formula, integer arguments only (orbital
angular momentum coupled to the total quark spin `S = 0` or `1`).
"""
function clebsch_gordan(j1::Integer, m1::Integer, j2::Integer, m2::Integer,
                        J::Integer, M::Integer)
    m1 + m2 == M || return 0.0
    (abs(m1) <= j1 && abs(m2) <= j2 && abs(M) <= J) || return 0.0
    (abs(j1 - j2) <= J <= j1 + j2) || return 0.0
    pref = sqrt((2J + 1) * factorial(J + j1 - j2) * factorial(J - j1 + j2) *
                factorial(j1 + j2 - J) / factorial(j1 + j2 + J + 1))
    pref *= sqrt(factorial(J + M) * factorial(J - M) *
                 factorial(j1 - m1) * factorial(j1 + m1) *
                 factorial(j2 - m2) * factorial(j2 + m2))
    total = 0.0
    for k = 0:(j1+j2-J)
        a, b = j1 - m1 - k, j2 + m2 - k
        c, d = J - j2 + m1 + k, J - j1 - m2 + k
        (a >= 0 && b >= 0 && c >= 0 && d >= 0) || continue
        total += (-1)^k / (factorial(k) * factorial(j1 + j2 - J - k) *
                           factorial(a) * factorial(b) * factorial(c) * factorial(d))
    end
    return pref * total
end

# -----------------------------------------------------------------------------
# A physical state, resampled onto one shared radial grid
# -----------------------------------------------------------------------------

"""
One `(L, S)` channel of a physical state, with its **coefficient-weighted sum**
`U(r) = sum_n c_n u_n(r)` over the radial levels the mixing brought in.

Components sharing `(L, S)` and differing only in `n` are one channel: they
carry the same angular function, so they add as amplitudes at every `r` and only
their `r`-integral is orthogonal. Squaring them separately is a real error, and
the one `check_angular_normalization` catches.
"""
struct CloudChannel
    L::Int
    S::Int
    U::Vector{Float64}
end

"""Probability carried by this channel, `int U(r)^2 dr`."""
channel_weight(c::CloudChannel, h::Real) = sum(abs2, c.U) * h

"""Sample any `RadialWave` onto the shared plotting grid, normalized on it."""
u_samples(w::OscillatorWave, rgrid::AbstractVector{<:Real}) = sample_wave(w, rgrid).u

function u_samples(w::MeshWave, rgrid::AbstractVector{<:Real})
    out = similar(collect(float(rgrid)))
    for (i, r) in pairs(rgrid)
        if r <= w.r[1]
            out[i] = w.u[1] * r / w.r[1]          # u(0) = 0, linear into the first node
        elseif r >= w.r[end]
            out[i] = 0.0
        else
            k = searchsortedlast(w.r, r)
            t = (r - w.r[k]) / (w.r[k+1] - w.r[k])
            out[i] = (1 - t) * w.u[k] + t * w.u[k+1]
        end
    end
    h = float(rgrid[2]) - float(rgrid[1])
    out ./= sqrt(sum(abs2, out) * h)
    return out
end

"""
    cloud_state(spec, label, rgrid) -> (J, channels)

Collapse a physical state into `CloudChannel`s on `rgrid`, one per distinct
`(L, S)`. Every component of an intra-meson mixed state shares `J`; `L` and `S`
need not agree, and several radial levels may land in the same channel.
"""
function cloud_state(spec, label::AbstractString, rgrid::AbstractVector{<:Real})
    pieces = physical_components(spec, label)
    J = first(pieces).basis.J
    all(p.basis.J == J for p in pieces) ||
        error("$label mixes different J; the mJ argument would be ambiguous")

    accumulated = Dict{Tuple{Int,Int},Vector{Float64}}()
    for piece in pieces
        key = (GIModel.L_SYMBOLS[piece.basis.L_label], (piece.basis.multiplicity - 1) ÷ 2)
        U = get!(() -> zeros(length(rgrid)), accumulated, key)
        U .+= piece.coefficient .* u_samples(piece.wave, rgrid)
    end
    channels = [CloudChannel(L, S, U) for ((L, S), U) in sort(collect(accumulated))]
    return J, channels
end

"""
    angular_profile(channels, k, J, mJ, x)

`f(x)` at `x = cos(theta)` and radial index `k`: the density with the common
`1/r^2` dropped, summed incoherently over spin `S` and its projection `m_s` and
coherently over `L` within each spin group. Different `L` at the same `S` share
`m_L = mJ - m_s`, so their interference is real and survives the `m_s` sum.
"""
function angular_profile(channels::Vector{CloudChannel}, k::Int,
                         J::Int, mJ::Int, x::Real)
    total = 0.0
    for S = 0:1
        any(c.S == S for c in channels) || continue
        for ms = (-S):S
            mL = mJ - ms
            amp = 0.0
            for c in channels
                c.S == S && abs(mL) <= c.L || continue
                amp += c.U[k] * clebsch_gordan(c.L, mL, S, ms, J, mJ) *
                       legendre_Y(c.L, mL, x)
            end
            total += amp^2
        end
    end
    return total
end

"""Exact radial marginal `p(r) = sum_{L,S} U_{LS}(r)^2`, up to normalization."""
radial_weights(channels::Vector{CloudChannel}) =
    [sum(c.U[k]^2 for c in channels) for k in eachindex(first(channels).U)]

"""
Check the angular machinery against the exact radial marginal: orthonormality of
the `Y_lm` forces `2 pi int_{-1}^{1} f(x) dx == sum_c c_c^2 u_c(r)^2` at every
radius. Returns the worst relative deviation over the radii carrying weight.
"""
function check_angular_normalization(channels, J, mJ, weights; nx = 2001)
    xs = range(-1, 1; length = nx)
    dx = step(xs)
    peak = maximum(weights)
    worst = 0.0
    for k in eachindex(weights)
        weights[k] > 1e-6 * peak || continue
        values = [angular_profile(channels, k, J, mJ, x) for x in xs]
        integral = 2pi * dx * (sum(values) - (first(values) + last(values)) / 2)
        worst = max(worst, abs(integral - weights[k]) / weights[k])
    end
    return worst
end

# -----------------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------------

"""
    sample_cloud(channels, J, mJ, rgrid, npoints; cutaway, rng)

Draw `npoints` positions from `|psi|^2` in fm. `r` comes from the exact radial
marginal by CDF inversion; `cos(theta)` by rejection against a per-radius
envelope; `phi` is uniform. `cutaway` removes the quadrant facing the camera
(`x < 0, y < 0` at [`AZIMUTH`](@ref)) so the radial nodes inside are visible.
"""
function sample_cloud(channels::Vector{CloudChannel}, J::Int, mJ::Int,
                      rgrid::AbstractVector{<:Real}, npoints::Int;
                      cutaway::Bool = true, rng = Random.default_rng())
    cdf = cumsum(radial_weights(channels))
    cdf ./= cdf[end]

    xs_probe = range(-1, 1; length = 257)
    envelope = [1.05 * maximum(angular_profile(channels, k, J, mJ, x) for x in xs_probe)
                for k in eachindex(rgrid)]

    h = float(rgrid[2]) - float(rgrid[1])
    px, py, pz = Float64[], Float64[], Float64[]
    sizehint!.((px, py, pz), npoints)
    while length(px) < npoints
        k = searchsortedfirst(cdf, rand(rng))
        envelope[k] > 0 || continue
        r = (float(rgrid[k]) + (rand(rng) - 0.5) * h) * HBARC_FM
        cosθ = 0.0
        while true
            cosθ = 2rand(rng) - 1
            rand(rng) * envelope[k] <= angular_profile(channels, k, J, mJ, cosθ) && break
        end
        ϕ = 2pi * rand(rng)
        sinθ = sqrt(max(0.0, 1 - cosθ^2))
        x, y, z = r * sinθ * cos(ϕ), r * sinθ * sin(ϕ), r * cosθ
        cutaway && x < 0 && y < 0 && continue   # the quadrant facing the camera
        push!(px, x); push!(py, y); push!(pz, z)
    end
    return px, py, pz
end

# -----------------------------------------------------------------------------
# Rendering
# -----------------------------------------------------------------------------

const AZIMUTH = 1.18pi
const ELEVATION = 0.22pi
const CLOUD_COLOR = RGBAf(0.80, 0.85, 0.06, 0.50)   # the "picture of an atom" yellow-green

"""Painter's algorithm: CairoMakie draws in array order, so sort back to front."""
function depth_sorted(px, py, pz)
    view = (cos(ELEVATION) * cos(AZIMUTH), cos(ELEVATION) * sin(AZIMUTH), sin(ELEVATION))
    order = sortperm([view[1] * x + view[2] * y + view[3] * z
                      for (x, y, z) in zip(px, py, pz)])
    return px[order], py[order], pz[order]
end

"""
    cloud_axis(slot, title, half_width)

A bare `Axis3`: no box, no ticks, no grid. Every panel gets the same limits, so
the clouds are directly comparable in size, and the only scale cue is the bar
drawn by [`scale_bar!`](@ref).
"""
function cloud_axis(slot, title, half_width)
    ax = Axis3(slot; aspect = :data, title = title, titlesize = 15, titlegap = 2,
               azimuth = AZIMUTH, elevation = ELEVATION, perspectiveness = 0.30,
               viewmode = :fitzoom, protrusions = 0)
    limits!(ax, -half_width, half_width, -half_width, half_width,
            -half_width, half_width)
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

"""Draw a `length_fm` bar along x at the bottom-front corner of the box."""
function scale_bar!(ax, half_width; length_fm = 0.5)
    x0, y0, z0 = -0.9half_width, 0.9half_width, -0.95half_width
    lines!(ax, [Point3f(x0, y0, z0), Point3f(x0 + length_fm, y0, z0)];
           color = :gray30, linewidth = 2.5)
    text!(ax, [Point3f(x0 + length_fm / 2, y0, z0)]; text = ["$(length_fm) fm"],
          color = :gray30, fontsize = 12, align = (:center, :top), offset = (0, -4))
    return ax
end

function draw_cloud!(ax, px, py, pz; markersize)
    sx, sy, sz = depth_sorted(px, py, pz)
    scatter!(ax, sx, sy, sz; markersize = markersize, color = CLOUD_COLOR,
             strokewidth = 0)
    return ax
end

"""Half-box that contains `quantile`-fraction of the cloud, rounded up for ticks."""
function cloud_extent(px, py, pz; quantile = 0.99)
    radii = sort(sqrt.(px .^ 2 .+ py .^ 2 .+ pz .^ 2))
    r = radii[clamp(ceil(Int, quantile * length(radii)), 1, length(radii))]
    return ceil(r * 5) / 5
end

# -----------------------------------------------------------------------------
# Charmonium clouds
# -----------------------------------------------------------------------------

params, quark_masses = load_parameters_and_quark_masses(PARAMS_PATH)
ccbar = Meson(quark_masses, :c, :c)
levels = spectrum_levels(2; L_labels = ("S", "P", "D"))
spec = compute_spectrum(params, ccbar; levels = levels, solver = OscillatorSolver())

rgrid = range(0.0, RMAX_INV_GEV; length = NRADIAL)
rng = Random.MersenneTwister(20260908)

# label, mJ, what the panel is there to show
const PANELS = [
    ("1^3S_1", 0, "L = 0: no angular structure"),
    ("2^3S_1", 0, "one radial node, opened by the cut"),
    ("1^1P_1", 0, "|Y_10|^2 dumbbell along z"),
    ("1^1P_1", 1, "same level, mJ = 1: a torus"),
    ("1^3P_2", 0, "J = 2 out of L = S = 1: prolate"),
    ("1^3D_1", 0, "3D1, incl. its 3S1 tensor admixture"),
]

const NPOINTS = 60_000

println("Charmonium density clouds  (m_c = $(round(quark_masses["c"]; digits = 3)) GeV)\n")
@printf("%-9s %-4s %8s %10s %10s   %s\n",
        "state", "mJ", "M [GeV]", "r_rms [fm]", "ang.check", "channel weights")

clouds = Vector{Any}(undef, length(PANELS))
for (i, (label, mJ, _)) in pairs(PANELS)
    J, channels = cloud_state(spec, label, rgrid)
    weights = radial_weights(channels)
    residual = check_angular_normalization(channels, J, mJ, weights)
    clouds[i] = sample_cloud(channels, J, mJ, rgrid, NPOINTS; rng = rng)

    h = step(rgrid)
    rrms = sqrt(sum(weights .* collect(rgrid) .^ 2) / sum(weights)) * HBARC_FM
    breakdown = join([@sprintf("%.4f %s%d", channel_weight(c, h),
                               GIModel.L_LABELS[c.L], 2c.S + 1)
                      for c in channels], "  ")
    @printf("%-9s %-4d %8.3f %10.3f %10.1e   %s\n",
            label, mJ, spectrum_state(spec, label).mass_GeV, rrms, residual, breakdown)
end

# One box for all six panels, so the 1S/2S size difference is the real one.
const HALF = maximum(cloud_extent(cloud...) for cloud in clouds)

fig = Figure(size = (1500, 1020), backgroundcolor = :white)
for (i, (label, mJ, caption)) in pairs(PANELS)
    row, col = fldmod1(i, 3)
    ax = cloud_axis(fig[row, col], "$label,  mJ = $mJ\n$caption", HALF)
    draw_cloud!(ax, clouds[i]...; markersize = 2.6)
    scale_bar!(ax, HALF)
end
Label(fig[0, 1:3],
      "Charmonium: quark-antiquark separation density |psi(r)|^2, one dot per sample" *
      "  —  common scale, x > 0 & y > 0 quadrant cut away";
      fontsize = 19, padding = (0, 0, 4, 0))
rowgap!(fig.layout, 4)

figdir = joinpath(@__DIR__, "figures")
mkpath(figdir)
panels_path = joinpath(figdir, "density_cloud_panels.png")
save(panels_path, fig)

# A single large 2^3S_1: the radial node is the point of the cutaway.
J_hero, channels_hero = cloud_state(spec, "2^3S_1", rgrid)
px, py, pz = sample_cloud(channels_hero, J_hero, 0, rgrid, 160_000; rng = rng)
half_hero = cloud_extent(px, py, pz)
hero = Figure(size = (950, 950), backgroundcolor = :white)
ax = cloud_axis(hero[1, 1], "charmonium 2^3S_1  (psi(2S)),  mJ = 0", half_hero)
draw_cloud!(ax, px, py, pz; markersize = 2.8)
scale_bar!(ax, half_hero)
hero_path = joinpath(figdir, "density_cloud_2S.png")
save(hero_path, hero)

println("\nwrote $(relpath(panels_path, dirname(@__DIR__)))")
println("wrote $(relpath(hero_path, dirname(@__DIR__)))")
