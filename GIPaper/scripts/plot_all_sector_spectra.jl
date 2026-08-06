#!/usr/bin/env julia
# Save reference-vs-model ladder plots for all digitized GI sectors, before and
# after the currently assigned same-J antisymmetric spin-orbit mixing.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie
using DataFrames
using LaTeXStrings
using Printf

root = dirname(@__DIR__)
using GIModel
using GIPaper

const OUTDIR = joinpath(@__DIR__, "spectrum_plots")

# Wong2 colorblind-safe palette (Okabe-Ito, black last) and a boxed,
# Computer-Modern theme to match the report style.
const WONG2 = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
]
const PAPER_COLOR = "#0072B2"     # wong2 blue
const COMPUTED_COLOR = "#D55E00"  # wong2 vermillion

set_theme!(
    merge(
        theme_latexfonts(),
        Theme(
            fontsize = 16,
            Axis = (
                xticksmirrored = true,
                yticksmirrored = true,
                xtickalign = 1,
                ytickalign = 1,
                xminortickalign = 1,
                yminortickalign = 1,
                topspinevisible = true,
                rightspinevisible = true,
                xgridvisible = true,
                ygridvisible = true,
            ),
        ),
    ),
)

const SPECTRUM_FILES = [
    ("isovector", joinpath(root, "data", "reference_spectrum_isovector.csv")),
    ("strange", joinpath(root, "data", "reference_spectrum_strange.csv")),
    ("isoscalar", joinpath(root, "data", "reference_spectrum_isoscalar.csv")),
    ("charmonium", joinpath(root, "data", "reference_spectrum_charmonium.csv")),
    ("charmed", joinpath(root, "data", "reference_spectrum_charmed.csv")),
    ("bottomonium", joinpath(root, "data", "reference_spectrum_bottomonium.csv")),
    ("b_flavored", joinpath(root, "data", "reference_spectrum_b_flavored.csv")),
]

const L_ORDER = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)

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
        for i = 1:(n-1)
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

function row_key(row)
    return (String(row.sector), String(row.state), row.n, row.multiplicity, String(row.L), row.J)
end

function state_label(row)
    return @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
end

function sector_label(row)
    return String(row.sector)
end

function spectroscopic_column(row)
    Lval = get(L_ORDER, String(row.L), 0)
    Sspin = row.multiplicity == 3 ? 1 : 0
    pch = iseven(Lval + 1) ? '+' : '-'
    if isapprox(row.m1_GeV, row.m2_GeV; rtol = 0.0, atol = 1e-12)
        cch = iseven(Lval + Sspin) ? '+' : '-'
        return string(row.J, pch, cch)
    end
    return string(row.J, pch)
end

function panel_labels(rows)
    labels = unique(sector_label(row) for row in rows)
    return sort!(collect(labels))
end

function j_labels(rows)
    return unique(spectroscopic_column(row) for row in rows)
end

"""Render a `J^{PC}` column string such as `"1--"` as a LaTeX tick label."""
function latexify_column(s::AbstractString)
    m = match(r"^(\d+)([+-]+)$", s)
    m === nothing && return latexstring(s)
    return latexstring(m.captures[1], "^{", m.captures[2], "}")
end

function axis_title(label)
    return replace(label, "_" => " ")
end

function mixed_rows(rows, params)
    bykey = Dict(row_key(row) => i for (i, row) in enumerate(rows))
    row_by_key = Dict(row_key(row) => row for row in rows)
    pred = [row.predicted_GeV for row in rows]
    touched = falses(length(rows))
    solve_cache = Dict{Tuple{ConstituentMasses,String},Any}()

    groups = Dict{Tuple{String,Int,String,Int,Float64,Float64},Vector{Any}}()
    for row in rows
        Lval = get(L_ORDER, String(row.L), -1)
        Lval > 0 || continue
        row.J == Lval || continue
        key = (String(row.sector), row.n, String(row.L), row.J, row.m1_GeV, row.m2_GeV)
        push!(get!(groups, key, Any[]), row)
    end

    for (_, group) in groups
        singlets = [row for row in group if row.multiplicity == 1]
        triplets = [row for row in group if row.multiplicity == 3]
        (isempty(singlets) || isempty(triplets)) && continue
        for singlet in singlets, triplet in triplets
            skey = row_key(singlet)
            tkey = row_key(triplet)
            masses = ConstituentMasses(singlet.m1_GeV, singlet.m2_GeV)
            sol = get!(solve_cache, (masses, String(singlet.L))) do
                ev, vecs, r = channel_solution(
                    params, masses, L_ORDER[String(singlet.L)];
                    nlevels = 6,
                    solver = FiniteDifferenceSolver(kinetic = :relativistic),
                )
                (eigenvalues_GeV = ev, vectors = vecs, r = r)
            end
            singlet.n <= size(sol.vectors, 2) || continue
            wave = RadialWaveOnUniformMesh(sol.vectors[:, singlet.n], sol.r)
            offdiag = spin_orbit_mixing_components(
                params,
                masses,
                String(singlet.L),
                wave;
                k_spin_orbit = params.fine_structure.k_spin_orbit,
            ).total
            isapprox(offdiag, 0.0; atol = 1e-12, rtol = 0.0) && continue
            mix = same_j_mixing(singlet.predicted_GeV, triplet.predicted_GeV, offdiag)
            si = bykey[skey]
            ti = bykey[tkey]
            refs = [singlet.reference_GeV, triplet.reference_GeV]
            if refs[1] <= refs[2]
                pred[si], pred[ti] = mix.masses[1], mix.masses[2]
            else
                pred[si], pred[ti] = mix.masses[2], mix.masses[1]
            end
            touched[si] = true
            touched[ti] = true
        end
    end

    return [
        merge(
            row_by_key[row_key(row)],
            (
                predicted_GeV = pred[i],
                residual_MeV = 1000 * (pred[i] - row.reference_GeV),
                mixed = touched[i],
            ),
        ) for (i, row) in enumerate(rows)
    ]
end

function plot_sector_grid(all_rows, path; title, mixed_mode = false)
    panels = panel_labels(all_rows)
    ncols = length(panels) <= 2 ? length(panels) : 2
    nrows = cld(length(panels), ncols)
    max_labels = maximum(length(j_labels([row for row in all_rows if sector_label(row) == panel])) for panel in panels)
    panel_width = length(panels) == 1 ? max(760, 35 * max_labels + 180) : 520
    panel_height = length(panels) == 1 ? 520 : 390
    fig = Figure(size = (panel_width * ncols, panel_height * nrows))
    palette = WONG2
    single_panel = length(panels) == 1

    for (pi, panel) in enumerate(panels)
        prow = (pi - 1) ÷ ncols + 1
        pcol = (pi - 1) % ncols + 1
        rows = [row for row in all_rows if sector_label(row) == panel]
        labels = j_labels(rows)
        xmap = Dict(label => i for (i, label) in enumerate(labels))
        xs_center = [Float64(xmap[spectroscopic_column(row)]) for row in rows]
        y_ref = [row.reference_GeV for row in rows]
        y_pred = [row.predicted_GeV for row in rows]

        ax = Axis(
            fig[prow, pcol];
            title = single_panel ? "" : axis_title(panel),
            xlabel = L"J^{P(C)}",
            ylabel = L"\mathrm{Mass~(GeV)}",
            xticks = (1:length(labels), latexify_column.(labels)),
            xticklabelrotation = deg2rad(40),
            yminorticksvisible = true,
            yminorgridvisible = true,
        )
        ylo = minimum(vcat(y_ref, y_pred)) - 0.12
        yhi = maximum(vcat(y_ref, y_pred)) + 0.12
        ylims!(ax, ylo, yhi)
        xlims!(ax, 0.45, length(labels) + 0.55)

        groups = Dict{Tuple{Int,String},Vector{Int}}()
        for i in eachindex(rows)
            push!(get!(groups, (rows[i].n, String(rows[i].L)), Int[]), i)
        end
        group_keys = sort!(collect(keys(groups)); by = key -> (get(L_ORDER, key[2], 99), key[1]))
        δ = 0.085
        for (gi, key) in enumerate(group_keys)
            idxs = groups[key]
            xmin = minimum(xs_center[i] - δ for i in idxs) - 0.14
            xmax = maximum(xs_center[i] + δ for i in idxs) + 0.14
            ymin = minimum(min(y_ref[i], y_pred[i]) for i in idxs) - 0.035
            ymax = maximum(max(y_ref[i], y_pred[i]) for i in idxs) + 0.035
            if xmax - xmin < 0.32
                mid = (xmin + xmax) / 2
                xmin, xmax = mid - 0.16, mid + 0.16
            end
            if ymax - ymin < 0.048
                mid = (ymin + ymax) / 2
                ymin, ymax = mid - 0.024, mid + 0.024
            end
            base = palette[mod1(gi, length(palette))]
            poly!(
                ax,
                rounded_rect_polygon(xmin, ymin, xmax, ymax, 0.11);
                color = (base, 0.22),
                strokecolor = (base, 0.58),
                strokewidth = 0.5,
                shading = NoShading,
            )
            text!(
                ax,
                xmax,
                ymin;
                text = latexstring(key[1], key[2]),
                align = (:center, :center),
                fontsize = 14,
                color = :black,
            )
        end

        for i in eachindex(rows)
            lines!(
                ax,
                [xs_center[i] - δ, xs_center[i] + δ],
                [y_ref[i], y_pred[i]];
                color = (:gray42, 0.45),
                linewidth = 0.85,
            )
        end
        scatter!(
            ax,
            xs_center .- δ,
            y_ref;
            label = "paper",
            markersize = 9,
            color = (PAPER_COLOR, 0.9),
            strokewidth = 0.45,
            strokecolor = :white,
        )
        scatter!(
            ax,
            xs_center .+ δ,
            y_pred;
            label = mixed_mode ? "computed + mix" : "computed",
            markersize = 9,
            color = (COMPUTED_COLOR, 0.95),
            strokewidth = 0.45,
            strokecolor = :white,
        )
        if mixed_mode
            touched_x = [
                xs_center[i] + δ for i in eachindex(rows) if
                hasproperty(rows[i], :mixed) && rows[i].mixed
            ]
            touched_y = [
                y_pred[i] for i in eachindex(rows) if
                hasproperty(rows[i], :mixed) && rows[i].mixed
            ]
            if !isempty(touched_x)
                scatter!(
                    ax,
                    touched_x,
                    touched_y;
                    label = "mixed rows",
                    marker = :xcross,
                    markersize = 12,
                    color = :black,
                )
            end
        end
        pi == 1 && axislegend(ax; position = :rt)
    end
    # Multi-panel grids keep per-axis sector titles for identification; the
    # single-panel report figures carry no title (described in the caption).
    if !single_panel
        Label(fig[0, :], title; fontsize = 22, font = :bold)
    end
    save(path, fig)
end

function safe_slug(label)
    return replace(lowercase(label), r"[^a-z0-9]+" => "_")
end

function compute_all_rows(params, mq)
    all_plain = NamedTuple[]
    all_mixed = NamedTuple[]
    for (label, path) in SPECTRUM_FILES
        reference = GIPaper.load_reference_spectrum(path)
        rows = compare_reference(
            params,
            mq,
            reference;
            contact_hyperfine = true,
            use_fine_structure = params.fine_structure.enabled,
            kinetic = :relativistic,
        )
        append!(all_plain, [merge(row, (mixed = false,)) for row in rows])
        append!(all_mixed, mixed_rows(rows, params))
        println(@sprintf("%-12s %3d rows", label, length(rows)))
    end
    return all_plain, all_mixed
end

mkpath(OUTDIR)
params, mq = load_parameters_and_quark_masses(joinpath(dirname(root), "data", "parameters.provisional.toml"))
plain_rows, mixed_rows_all = compute_all_rows(params, mq)

before_path = joinpath(OUTDIR, "all_sectors_before_mixing.png")
after_path = joinpath(OUTDIR, "all_sectors_after_available_mixing.png")

plot_sector_grid(
    plain_rows,
    before_path;
    title = "All GI sectors: paper labels vs computed masses (before mixing)",
)
plot_sector_grid(
    mixed_rows_all,
    after_path;
    title = "All GI sectors: paper labels vs computed masses (after available same-J mixing)",
    mixed_mode = true,
)

for panel in panel_labels(plain_rows)
    plain_panel = [row for row in plain_rows if sector_label(row) == panel]
    mixed_panel = [row for row in mixed_rows_all if sector_label(row) == panel]
    before_panel_path = joinpath(OUTDIR, string(safe_slug(panel), "_before_mixing.png"))
    after_panel_path = joinpath(OUTDIR, string(safe_slug(panel), "_after_available_mixing.png"))
    plot_sector_grid(
        plain_panel,
        before_panel_path;
        title = string(axis_title(panel), ": before mixing"),
    )
    plot_sector_grid(
        mixed_panel,
        after_panel_path;
        title = string(axis_title(panel), ": after available same-J mixing"),
        mixed_mode = true,
    )
    println("wrote ", before_panel_path)
    println("wrote ", after_panel_path)
end

println("wrote ", before_path)
println("wrote ", after_path)
