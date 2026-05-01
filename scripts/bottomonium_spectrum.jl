#!/usr/bin/env julia
# Prints GI-model bottomonium masses (GeV) for Fig. 8 states, saves a PNG ladder plot (reference vs model).
# SVG ladders matching the paper figures live in scripts/plot_spectrum_digitizations.py (Fig. 8 -> fig08_bottomonia).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie
using DataFrames

root = dirname(@__DIR__)
using GIModel

"""
Spectroscopic label J^PC from orbital L, spin multiplicity (1 singlet / 3 triplet), and total J.
Uses P = (-1)^(L+1) and C = (-1)^(L+S) for a fermion–antifermion pair.
"""
function spectroscopic_jpc(L::AbstractString, multiplicity::Int, J::Int)
    ℓ_dict = Dict('S' => 0, 'P' => 1, 'D' => 2, 'F' => 3, 'G' => 4)
    ℓ = ℓ_dict[only(L)]
    Sspin = multiplicity == 3 ? 1 : 0
    pch = iseven(ℓ + 1) ? '+' : '-'
    cch = iseven(ℓ + Sspin) ? '+' : '-'
    return string(J, pch, cch)
end

# Matches scripts/plot_spectrum_digitizations.py CONFIG["08"]["order"] (GI Fig. 8 bottomonium ladder).
const FIG8_BOTTOMONIUM_JPC = [
    "0-+",
    "1--",
    "1+-",
    "0++",
    "1++",
    "2++",
    "2-+",
    "2--",
    "3--",
    "3+-",
    "3++",
    "4++",
]

"""Closed polygon for an axis-aligned rectangle with circular corners (CCW, data coordinates)."""
function rounded_rect_polygon(
    xmin::Float64,
    ymin::Float64,
    xmax::Float64,
    ymax::Float64,
    r::Float64,
)
    r = min(r, (xmax - xmin) / 2, (ymax - ymin) / 2)
    n = 14
    pts = Point2f[]
    function quad_arc(cx::Float64, cy::Float64, θ0::Float64, θ1::Float64)
        for i in 1:(n-1)
            θ = θ0 + (θ1 - θ0) * i / n
            push!(pts, Point2f(cx + r * cos(θ), cy + r * sin(θ)))
        end
    end
    push!(pts, Point2f(xmin + r, ymin))
    push!(pts, Point2f(xmax - r, ymin))
    quad_arc(xmax - r, ymin + r, Float64(3π / 2), Float64(2π))
    push!(pts, Point2f(xmax, ymax - r))
    quad_arc(xmax - r, ymax - r, 0.0, Float64(π / 2))
    push!(pts, Point2f(xmax - r, ymax))
    push!(pts, Point2f(xmin + r, ymax))
    quad_arc(xmin + r, ymax - r, Float64(π / 2), Float64(π))
    push!(pts, Point2f(xmin, ymin + r))
    quad_arc(xmin + r, ymin + r, Float64(π), Float64(3π / 2))
    push!(pts, pts[1])
    return pts
end

"""Ladder plot: horizontal J^PC columns like Fig. 8; reference vs computed masses."""
function plot_bottomonium_ref_vs_model(rows; path::AbstractString)
    n = length(rows)
    xs_center = zeros(Float64, n)
    y_ref = zeros(Float64, n)
    y_pred = zeros(Float64, n)
    @inbounds for i in 1:n
        row = rows[i]
        jpc = spectroscopic_jpc(row.L, row.multiplicity, row.J)
        col = findfirst(==(jpc), FIG8_BOTTOMONIUM_JPC)
        col === nothing &&
            error("state row $(i): reconstructed JPC=$(repr(jpc)) not in Fig. 8 column list")
        xs_center[i] = Float64(col)
        y_ref[i] = row.reference_GeV
        y_pred[i] = row.predicted_GeV
    end

    δ = 0.085
    x_pad = 0.14
    y_pad = 0.035
    min_box_w = 0.32
    min_box_h = 0.048
    corner_r = 0.11
    fill_alpha = 0.26

    sectors = Dict{Tuple{Int,String},Vector{Int}}()
    for i in 1:n
        key = (rows[i].n, rows[i].L)
        push!(get!(sectors, key, Int[]), i)
    end
    sector_list = sort!(collect(keys(sectors)); by=kl -> (kl[2], kl[1]))
    palette = Makie.wong_colors()

    fig = Figure(size=(1000, 720))
    ax = Axis(
        fig[1, 1];
        xlabel="JPC",
        ylabel="Mass (GeV)",
        title="Bottomonium: reference vs GI model (Fig. 8 states)",
        xticks=(1:length(FIG8_BOTTOMONIUM_JPC), FIG8_BOTTOMONIUM_JPC),
        xticklabelrotation=deg2rad(40),
        yminorticksvisible=true,
        yminorgridvisible=true,
    )
    ylims!(ax, 9.15, 11.32)

    for (si, key) in enumerate(sector_list)
        idxs = sectors[key]
        xmin = minimum(xs_center[i] - δ for i in idxs) - x_pad
        xmax = maximum(xs_center[i] + δ for i in idxs) + x_pad
        ymin = minimum(min(y_ref[i], y_pred[i]) for i in idxs) - y_pad
        ymax = maximum(max(y_ref[i], y_pred[i]) for i in idxs) + y_pad
        if xmax - xmin < min_box_w
            mid = (xmin + xmax) / 2
            xmin, xmax = mid - min_box_w / 2, mid + min_box_w / 2
        end
        if ymax - ymin < min_box_h
            mid = (ymin + ymax) / 2
            ymin, ymax = mid - min_box_h / 2, mid + min_box_h / 2
        end
        poly = rounded_rect_polygon(xmin, ymin, xmax, ymax, corner_r)
        base = palette[mod1(si, length(palette))]
        poly!(
            ax,
            poly;
            color=(base, fill_alpha),
            strokecolor=(base, 0.62),
            strokewidth=0.5,
            shading=NoShading,
        )
        n_radial, Lorb = key
        text!(
            ax,
            xmax,
            ymin;
            text=string(n_radial, Lorb),
            align=(:center, :center),
            fontsize=15,
            color=:black,
        )
    end

    for i in 1:n
        lines!(
            ax,
            [xs_center[i] - δ, xs_center[i] + δ],
            [y_ref[i], y_pred[i]];
            color=(:gray42, 0.45),
            linewidth=0.85,
        )
    end
    scatter!(
        ax,
        xs_center .- δ,
        y_ref;
        label="reference",
        markersize=10,
        color=(:steelblue4, 0.9),
        strokewidth=0.5,
        strokecolor=:white,
    )
    scatter!(
        ax,
        xs_center .+ δ,
        y_pred;
        label="computed",
        markersize=10,
        color=(:darkorange, 0.95),
        strokewidth=0.5,
        strokecolor=:white,
    )
    axislegend(ax; position=:rt)
    save(path, fig)
    return fig
end

params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
ref = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_bottomonium.csv"))
computed = compute_sector(params, ref, "b"; kinetic=:relativistic)
rows = compare(
    computed,
    ref;
    contact_hyperfine=true,
    use_fine_structure=params.fine_structure,
)

plot_path = joinpath(@__DIR__, "bottomonium_spectrum_ref_vs_model.png")
plot_bottomonium_ref_vs_model(rows; path=plot_path)
println("\nWrote spectrum plot: ", plot_path)

df = DataFrame(rows)
df = transform(
    df,
    [:n, :multiplicity, :L, :J] =>
        ByRow((n, mult, L, J) -> string(n, '^', mult, L, '_', J)) => :state,
    :predicted_GeV => ByRow(x -> round(x; digits=3)) => :model_GeV,
    :reference_GeV => ByRow(x -> round(x; digits=3)) => :ref_GeV,
    :residual_MeV => ByRow(x -> round(x; digits=1)) => :delta_MeV,
)
df = select(df, :state, :model_GeV, :ref_GeV, :delta_MeV)

println("Bottomonium (predicted masses, GeV)\n")
show(stdout, df; allrows=true, show_row_number=false)
println()
