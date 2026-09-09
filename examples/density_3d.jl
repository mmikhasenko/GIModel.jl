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
#   (3) exact radial marginal p(r) = sum_{L,S} U_LS(r)^2, where U_LS is the
#       coefficient-weighted sum of the u(r) that mixing put into that channel.
#       So r is drawn by CDF inversion and only cos(theta) needs rejection; phi
#       is uniform, since all L in one spin group share m_L = mJ - m_s and
#       therefore one common exp(i m_L phi).
#
# Interference between different L at the same S -- the 3S1/3D1 tensor mixing --
# survives step (2) and is visible as a mJ-dependent distortion of the cloud.
#
# Two run-time checks, both exact identities rather than reference numbers:
#   `check_angular_normalization`  2 pi int f(x) dx == p(r) at every radius,
#                                  from orthonormality of the Y_Lm
#   `check_mj_sum_isotropy`        summing the density over mJ must give a
#                                  sphere, which pins each Clebsch-Gordan
#                                  coefficient individually
#
# Run:  julia examples/density_3d.jl
# Writes: figures/density_cloud_panels.png, figures/density_cloud_2S.png
# =============================================================================

using Pkg
Pkg.activate(@__DIR__; io = devnull)     # examples/Project.toml — GIModel + GLMakie
Pkg.instantiate(; io = devnull)

using Printf
using Random
using GIModel
using GLMakie
using LinearAlgebra: normalize

# Real spheres with lighting and ambient occlusion, so `ssao` on the screen too.
GLMakie.activate!(visible = false, ssao = true, fxaa = true, px_per_unit = 2)

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

"""
    check_mj_sum_isotropy(channels, J; nx = 401)

Sum the density over all `mJ` and demand it come out spherical. `sum_mJ |J mJ><J mJ|`
commutes with rotations, so this holds for any state of definite `J`, including
an `L`-superposition -- and unlike [`check_angular_normalization`](@ref), which
only constrains the Clebsch-Gordan coefficients through `sum_ms CG^2 = 1`, it
pins down each coefficient individually. Returns the peak-to-peak spread of the
summed profile, relative to its mean.
"""
function check_mj_sum_isotropy(channels, J; nx = 401)
    k = argmax(radial_weights(channels))
    profile = [sum(angular_profile(channels, k, J, mJ, x) for mJ = (-J):J)
               for x in range(-1, 1; length = nx)]
    return (maximum(profile) - minimum(profile)) / (sum(profile) / nx)
end

# -----------------------------------------------------------------------------
# Sampling
# -----------------------------------------------------------------------------

"""
    sample_cloud(channels, J, mJ, rgrid, npoints; keep, rng)

Draw `npoints` accepted positions from `|psi|^2` in fm. `r` comes from the exact
radial marginal by CDF inversion, `cos(theta)` by rejection against a per-radius
envelope, `phi` uniform. `keep(x, y, z)` selects the region to render — see
[`whole_cloud`](@ref) and [`open_doors`](@ref).
"""
function sample_cloud(channels::Vector{CloudChannel}, J::Int, mJ::Int,
                      rgrid::AbstractVector{<:Real}, npoints::Int;
                      keep = whole_cloud, rng = Random.default_rng())
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
        keep(x, y, z) || continue
        push!(px, x); push!(py, y); push!(pz, z)
    end
    return px, py, pz
end

# -----------------------------------------------------------------------------
# Rendering
# -----------------------------------------------------------------------------

const AZIMUTH = 1.25pi
const ELEVATION = 0.13pi            # near the side: polar structure reads best
const EYE_DISTANCE = 2.05           # in units of the cloud's half-extent

const BACKGROUND = RGBf(1.0, 1.0, 1.0)
const CLOUD_COLOR = RGBf(0.85, 0.90, 0.10)


# Remove a camera-facing wedge to expose two radial faces. Interior particles
# remain present; a geometric depth tint separates them from the exposed faces.
# This tint is an illustrative depth cue, not a value of the wave function.
const DOOR_HALF_ANGLE = deg2rad(60)     # a 120-degree opening facing the camera

"""Keep every sample: the cloud as the object it is."""
whole_cloud(x, y, z) = true

"""
    open_doors(x, y, z)

Drop the wedge within [`DOOR_HALF_ANGLE`](@ref) of the camera's azimuth. Both cut
planes contain the z axis, so this is one wedge hinged on the vertical and opened
toward the viewer.
"""
function open_doors(x, y, z)
    Δ = mod(atan(y, x) - AZIMUTH + pi, 2pi) - pi
    return abs(Δ) > DOOR_HALF_ANGLE
end

"""
    headlamp(eye, half_width)

A point light at the camera, a dim fill, and ambient illumination. The short
range preserves depth falloff without saturating the yellow-green particles.
Sphere lighting and SSAO supply local texture; `cutaway_colors` supplies the
larger-scale cavity depth cue. These lights do not cast physical volume shadows.
"""
function headlamp(eye, half_width)
    fill = normalize(-eye + Vec3f(0.45, -0.25, -0.35))
    return [
        AmbientLight(RGBf(0.10, 0.10, 0.12)),
        PointLight(RGBf(8.0, 7.68, 6.72), eye, 6.0 * half_width),
        DirectionalLight(RGBf(0.22, 0.22, 0.20), fill),
    ]
end

"""
    cloud_scene(slot, half_width)

A bare `LScene` looking at the origin. The eye distance is tied to `half_width`,
so passing the same `half_width` to every panel puts them all on one scale.
"""
function cloud_scene(slot, half_width)
    # Camera settings are explicit so adding plots cannot reframe each panel.
    distance = EYE_DISTANCE * half_width
    eye = Vec3f(distance * cos(ELEVATION) * cos(AZIMUTH),
                distance * cos(ELEVATION) * sin(AZIMUTH),
                distance * sin(ELEVATION))
    ls = LScene(slot; show_axis = false,
                scenekw = (lights = headlamp(eye, half_width),
                           backgroundcolor = BACKGROUND, clear = true,
                           ssao = Makie.SSAO(radius = 0.18f0, bias = 0.003f0,
                                             blur = 3)))
    cam3d!(ls.scene; center = false, projectiontype = Makie.Perspective, fov = 55)
    update_cam!(ls.scene, eye, Vec3f(0, 0, 0), Vec3f(0, 0, 1))
    return ls
end

"""
    cutaway_colors(px, py; depth_scale = 0.16, shadow_floor = 0.24)

Illustrative depth tint, based only on distance behind the nearest exposed
half-plane (in fm). The half-plane ends at the z axis: points whose perpendicular
projection lies beyond that edge use distance to the axis. Brightness decays
smoothly from 1 at the cut to `shadow_floor` in the interior. No radial node locations or
wave-function values enter this shading, and no particles are moved or removed.
"""
function cutaway_colors(px, py; depth_scale = 0.16, shadow_floor = 0.24,
                        base_color = CLOUD_COLOR)
    0 <= shadow_floor <= 1 || throw(ArgumentError("shadow_floor must be in [0, 1]"))
    depth_scale > 0 || throw(ArgumentError("depth_scale must be positive"))
    return map(px, py) do x, y
        Δ = abs(mod(atan(y, x) - AZIMUTH + pi, 2pi) - pi)
        depth = hypot(x, y) * sin(clamp(Δ - DOOR_HALF_ANGLE, 0, pi / 2))
        brightness = shadow_floor + (1 - shadow_floor) * exp(-depth / depth_scale)
        RGBf(brightness * base_color.r, brightness * base_color.g,
             brightness * base_color.b)
    end
end

"""Draw the samples as lit spheres of radius `radius` fm. Opaque: the depth
buffer resolves occlusion, so no painter's-algorithm sorting is needed."""
function draw_cloud!(ls, px, py, pz; radius, color = CLOUD_COLOR)
    meshscatter!(ls, px, py, pz; markersize = radius, color = color,
                 shading = true, ssao = true,
                 diffuse = Vec3f(0.90), specular = Vec3f(0.55), shininess = 96.0f0)
    return ls
end

"""Draw a `length_fm` bar below the cloud; the only scale cue in a bare scene."""
function scale_bar!(ls, half_width; length_fm = 0.5, fontsize = 15)
    # In the image plane through r = 0, so the bar is horizontal and its length
    # is not foreshortened. Perspective scale is referenced to that plane.
    right = Vec3f(-sin(AZIMUTH), cos(AZIMUTH), 0)
    up = Vec3f(-sin(ELEVATION) * cos(AZIMUTH),
               -sin(ELEVATION) * sin(AZIMUTH), cos(ELEVATION))
    start = -0.65half_width * right - 0.94half_width * up
    stop = start + length_fm * right
    lines!(ls, [Point3f(start), Point3f(stop)]; color = :gray25, linewidth = 3)
    text!(ls, [Point3f((start + stop) / 2)]; text = ["$(length_fm) fm"],
          color = :gray25, fontsize = fontsize, align = (:center, :top), offset = (0, -6))
    return ls
end

"""Half-box that contains `quantile`-fraction of the cloud."""
function cloud_extent(px, py, pz; quantile = 0.99)
    radii = sort(sqrt.(px .^ 2 .+ py .^ 2 .+ pz .^ 2))
    r = radii[clamp(ceil(Int, quantile * length(radii)), 1, length(radii))]
    return ceil(r * 5) / 5
end

# -----------------------------------------------------------------------------
# Charmonium clouds
# -----------------------------------------------------------------------------

display_label(label) = replace(label, "^1" => "¹", "^3" => "³", "_0" => "₀",
                               "_1" => "₁", "_2" => "₂", "_3" => "₃")

function density_demo()
    params, quark_masses = load_parameters_and_quark_masses(PARAMS_PATH)
    ccbar = Meson(quark_masses, :c, :c)
    levels = spectrum_levels(2; L_labels = ("S", "P", "D"))
    spec = compute_spectrum(params, ccbar; levels = levels, solver = OscillatorSolver())

    rgrid = range(0.0, RMAX_INV_GEV; length = NRADIAL)
    rng = Random.MersenneTwister(20260908)

    # label, mJ, what the panel is there to show
    PANELS = [
        ("1^3S_1", 0, "Ground state · compact core"),
        ("2^3S_1", 0, "Radial excitation · inner core and outer shell"),
        ("1^1P_1", 0, "Polar lobes · |Y₁₀|²"),
        ("1^1P_1", 1, "Equatorial torus · |Y₁₁|²"),
        ("1^3P_2", 0, "Spin-coupled P wave · prolate"),
        ("1^3D_1", 0, "D wave · includes S-wave tensor mixing"),
    ]

    NPOINTS = 90_000
    # Ball size against point count is the whole readability trade-off: big enough
    # that the light gives each sample a highlight, sparse enough that the gaps let
    # you see the shells behind. Solid coverage hides everything but the silhouette.
    SPHERE_RADIUS = 0.010         # fm

    println("Charmonium density clouds  (m_c = $(round(quark_masses["c"]; digits = 3)) GeV)\n")
    @printf("%-9s %-4s %8s %10s %10s %10s   %s\n",
            "state", "mJ", "M [GeV]", "r_rms [fm]", "norm.chk", "isotropy", "channel weights")

    clouds = Vector{Any}(undef, length(PANELS))
    for (i, (label, mJ, _)) in pairs(PANELS)
        J, channels = cloud_state(spec, label, rgrid)
        weights = radial_weights(channels)
        residual = check_angular_normalization(channels, J, mJ, weights)
        isotropy = check_mj_sum_isotropy(channels, J)
        @assert residual < 1e-6 "Angular normalization failed for $label"
        @assert isotropy < 1e-12 "mJ-sum isotropy failed for $label"
        clouds[i] = sample_cloud(channels, J, mJ, rgrid, NPOINTS;
                                 keep = open_doors, rng = rng)

        h = step(rgrid)
        rrms = sqrt(sum(weights .* collect(rgrid) .^ 2) / sum(weights)) * HBARC_FM
        breakdown = join([@sprintf("%.4f %s%d", channel_weight(c, h),
                                   GIModel.L_LABELS[c.L], 2c.S + 1)
                          for c in channels], "  ")
        @printf("%-9s %-4d %8.3f %10.3f %10.1e %10.1e   %s\n",
                label, mJ, spectrum_state(spec, label).mass_GeV, rrms, residual, isotropy,
                breakdown)
    end

    # One eye distance for all six panels, so the 1S/2S size difference is the real one.
    HALF = maximum(cloud_extent(cloud...) for cloud in clouds)

    fig = Figure(size = (1500, 1120), backgroundcolor = BACKGROUND)
    for (i, (label, mJ, caption)) in pairs(PANELS)
        row, col = fldmod1(i, 3)
        ls = cloud_scene(fig[row, col], HALF)
        draw_cloud!(ls, clouds[i]...; radius = SPHERE_RADIUS,
                    color = cutaway_colors(clouds[i][1], clouds[i][2]))
        scale_bar!(ls, HALF)
        Label(fig[row, col, Top()], "$(display_label(label))   ·   mⱼ = $mJ\n$caption";
              fontsize = 16, font = :bold, padding = (0, 0, 2, 0), justification = :center)
    end
    Label(fig[0, 1:3], "CHARMONIUM  /  Separation probability density";
          fontsize = 25, font = :bold, padding = (0, 0, 0, 14))
    Label(fig[3, 1:3],
          "120° cutaway  ·  Common camera and scale  ·  Particle concentration represents |ψ(r)|²\n" *
          "Brightness indicates depth behind the cut faces; it is not a density colour scale.";
          fontsize = 15, color = :gray35, padding = (0, 0, 8, 8))
    rowgap!(fig.layout, 1, 26)

    figdir = joinpath(@__DIR__, "figures")
    mkpath(figdir)
    panels_path = joinpath(figdir, "density_cloud_panels.png")
    save(panels_path, fig)

    # The cutaway is an exact subset of the full draw, preserving particle density.
    J_hero, channels_hero = cloud_state(spec, "2^3S_1", rgrid)
    shut = sample_cloud(channels_hero, J_hero, 0, rgrid, 180_000; rng = rng)
    keep_hero = open_doors.(shut...)
    opened = map(v -> v[keep_hero], shut)
    half_hero = cloud_extent(shut...)

    hero = Figure(size = (1700, 1080), backgroundcolor = BACKGROUND)
    for (col, (points, caption, cut)) in pairs([
            (shut, "01   Full probability cloud", false),
            (opened, "02   Cutaway · radial depletion", true)])
        ls = cloud_scene(hero[1, col], half_hero)
        colors = cut ? cutaway_colors(points[1], points[2]) : CLOUD_COLOR
        draw_cloud!(ls, points...; radius = SPHERE_RADIUS, color = colors)
        scale_bar!(ls, half_hero)
        Label(hero[1, col, Top()], caption; fontsize = 19, font = :bold,
              padding = (0, 0, 6, 0))
    end
    Label(hero[0, 1:2], "Inside charmonium 2³S₁  ·  ψ(2S)  ·  mⱼ = 0";
          fontsize = 28, font = :bold, padding = (0, 0, 0, 14))

    # A quantitative companion to the shaded view. This is the full mixed state's
    # radial marginal, normalized per fm, not a histogram of the displayed cutaway.
    weights_hero = radial_weights(channels_hero)
    r_fm = collect(rgrid) .* HBARC_FM
    probability = weights_hero ./ (sum(weights_hero) * step(rgrid) * HBARC_FM)
    S_channel = only(filter(c -> c.L == 0, channels_hero))
    node_index = findfirst(k -> S_channel.U[k] * S_channel.U[k+1] < 0,
                           2:(length(rgrid)-1))
    # `findfirst` above returns the index within the range, whose first value is 2.
    node_index === nothing && error("Expected a radial node in the 2S channel")
    k_node = node_index + 1
    node_fm = (r_fm[k_node] * abs(S_channel.U[k_node+1]) +
               r_fm[k_node+1] * abs(S_channel.U[k_node])) /
              (abs(S_channel.U[k_node]) + abs(S_channel.U[k_node+1]))
    profile = Axis(hero[2, 1]; xlabel = "Separation r [fm]",
                   ylabel = "p(r) [fm⁻¹]", backgroundcolor = BACKGROUND,
                   xgridvisible = false, ygridcolor = (:gray60, 0.15),
                   xlabelsize = 17, ylabelsize = 17, xticklabelsize = 14, yticklabelsize = 14)
    band!(profile, r_fm, zeros(length(r_fm)), probability; color = (CLOUD_COLOR, 0.22))
    lines!(profile, r_fm, probability; color = RGBf(0.37, 0.43, 0.04), linewidth = 3)
    vlines!(profile, [node_fm]; color = :gray35, linestyle = :dash, linewidth = 1.5)
    xlims!(profile, 0, half_hero)
    ylims!(profile, 0, 1.08maximum(probability))
    hidespines!(profile, :t, :r)
    Label(hero[2, 2],
          "THE CALCULATED RADIAL STRUCTURE\n\n" *
          @sprintf("S-wave node at r ≈ %.3f fm\n", node_fm) *
          "The curve includes the tensor-mixed channels.\n" *
          "p(r) = r² ∫ |ψ(r)|² dΩ;  ∫ p(r) dr = 1.\n\n" *
          "The cutaway keeps the same sampled positions.\n" *
          "Particle concentration shows probability density;\n" *
          "brightness is an illustrative cue for depth.";
          fontsize = 17, color = :gray25, justification = :left, halign = :center,
          tellwidth = false)
    colsize!(hero.layout, 1, Relative(0.5))
    colsize!(hero.layout, 2, Relative(0.5))
    rowsize!(hero.layout, 2, Fixed(220))
    rowgap!(hero.layout, 1, 24)
    Label(hero[3, 1:2], "Quark–antiquark relative coordinate  ·  GI model with provisional parameters  ·  120° wedge removed";
          fontsize = 15, color = :gray40, padding = (0, 0, 6, 10))
    hero_path = joinpath(figdir, "density_cloud_2S.png")
    save(hero_path, hero)
    @printf("\n2S S-channel node: %.6f fm; radial probability integral: %.12f\n",
            node_fm, sum(probability) * step(rgrid) * HBARC_FM)

    println("\nwrote $(relpath(panels_path, dirname(@__DIR__)))")
    println("wrote $(relpath(hero_path, dirname(@__DIR__)))")

end # density_demo

abspath(PROGRAM_FILE) == abspath(@__FILE__) && density_demo()
