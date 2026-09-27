# Figures for the landing page, written to docs/src/public/home/.
# Run by docs/render.jl in the docs/quarto environment:
#
#   julia --project=docs/quarto docs/quarto/home_figures.jl

using GIModel
using CairoMakie
using Random

const OUT = joinpath(dirname(@__DIR__), "src", "public", "home")
mkpath(OUT)

params, mq = load_parameters_and_quark_masses(default_parameters_path())
solver = OscillatorSolver()

# --- 1. A wavefunction as a cloud of points ------------------------------------
#
# Positions are drawn from |ψ(r)|² = |u(r)/r|² |Y_20(θ)|² for charmonium 1¹D₂
# (m = 0), and each point is coloured by the sign of ψ.

charmonium = compute_spectrum(params, Meson(mq, :c, :c);
    levels = [BasisState(1, "D", 1, 2)], solver)
wave = sample_wave(radial_wave(charmonium, "1^1D_2"), range(0, 12; length = 2401))

rng = MersenneTwister(1985)
cdf = cumsum(wave.u .^ 2); cdf ./= cdf[end]
draw_r() = wave.r[searchsortedfirst(cdf, rand(rng))]
ylm(c) = 3c^2 - 1                        # Y_20 up to normalization
function draw_cos()
    while true
        c = 2rand(rng) - 1
        rand(rng) * 4 <= ylm(c)^2 && return c
    end
end

npoints = 30_000
points = Point3f[]
signs = Float32[]
for _ in 1:npoints
    r, c, φ = draw_r(), draw_cos(), 2π * rand(rng)
    s = sqrt(1 - c^2)
    push!(points, Point3f(r * s * cos(φ), r * s * sin(φ), r * c))
    push!(signs, sign(ylm(c)))
end

fig = Figure(size = (800, 800), backgroundcolor = :transparent, figure_padding = 0)
ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 0.35π, elevation = 0.12π,
    backgroundcolor = :transparent, protrusions = 0, viewmode = :fit)
hidedecorations!(ax); hidespines!(ax)
limits!(ax, -4.3, 4.3, -4.3, 4.3, -4.3, 4.3)
colors = [s > 0 ? RGBAf(0.36, 0.48, 1.00, 0.55) : RGBAf(1.00, 0.58, 0.18, 0.55) for s in signs]
scatter!(ax, points; color = colors, markersize = 5, strokewidth = 0)
save(joinpath(OUT, "wavefunction_cloud.png"), fig; px_per_unit = 1)

# --- 2. The b c̄ spectrum, in the layout of Fig. 9 of the paper ----------------

levels = vcat(spectrum_levels(3; L_labels = ("S",)), spectrum_levels(2; L_labels = ("P",)),
    spectrum_levels(1; L_labels = ("D",)))
bc = compute_spectrum(params, Meson(mq, :b, :c); levels, solver)

# One column per J^P; the two J = L states of each multiplet share a column,
# because singlet–triplet mixing makes them physical superpositions.
columns = [
    ("0⁻", [("S", 1, 0)]), ("1⁻", [("S", 3, 1), ("D", 3, 1)]),
    ("0⁺", [("P", 3, 0)]), ("1⁺", [("P", 1, 1), ("P", 3, 1)]), ("2⁺", [("P", 3, 2)]),
    ("2⁻", [("D", 1, 2), ("D", 3, 2)]), ("3⁻", [("D", 3, 3)]),
]
fig = Figure(size = (900, 560), backgroundcolor = :white, fontsize = 18)
ax = Axis(fig[1, 1]; ylabel = "mass  [GeV]",
    title = rich("B", subscript("c"), " mesons (b c̄), Godfrey–Isgur model"),
    xticks = (1:length(columns), first.(columns)), xgridvisible = false,
    ylabelsize = 20, titlesize = 22)
for (x, (_, keys)) in enumerate(columns)
    for s in bc.states
        (s.L, s.multiplicity, s.J) in keys || continue
        colour = s.L == "S" ? :royalblue : s.L == "P" ? :darkorange : :seagreen
        lines!(ax, [x - 0.32, x + 0.32], fill(s.mass_GeV, 2); color = colour, linewidth = 3)
    end
end
threshold = 5.279 + 1.870          # B D
hlines!(ax, [threshold]; color = :gray50, linestyle = :dash)
text!(ax, 0.55, threshold + 0.02; text = "B D threshold", color = :gray40, fontsize = 15)
xlims!(ax, 0.4, length(columns) + 0.6)
axislegend(ax, [LineElement(color = c, linewidth = 3) for c in (:royalblue, :darkorange, :seagreen)],
    ["S wave", "P wave", "D wave"]; position = :rb, framevisible = false)
save(joinpath(OUT, "bc_spectrum.png"), fig; px_per_unit = 1.5)

println("wrote ", readdir(OUT))
