#!/usr/bin/env julia
# # Mixing Studies
#
# This literate notebook consolidates the current mixing studies in the
# Godfrey-Isgur reproduction. It can be read top-to-bottom as documentation or
# executed directly:
#
# ```bash
# julia --project=GIPaper GIPaper/docs/mixing_studies.jl
# ```
#
# The notebook regenerates `docs/residual_reports/mixing_mechanism_demo.md` and
# its PNG figures. It keeps comparison separate from computation: fixed-sector
# radial channels are solved first, mixing is assigned from cached
# `SectorComputation` data, and reference residuals are evaluated last.

using Pkg
Pkg.activate(joinpath(dirname(@__DIR__), "scripts"))

using CairoMakie
using LinearAlgebra
using Printf

root = dirname(@__DIR__)
using GIModel
using GIPaper

const OUTDIR = joinpath(root, "docs", "residual_reports", "figures")

# ## Helpers
#
# The helpers below load sectors, run `compare_reference` with selected mechanism
# toggles, and recover the concrete two-state or four-state blocks used in the figures.

function row_label(row)
    return @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
end

function sector_rows(params, mq, file; ngrid = 180, rmax = 18.0, aso = true, tensor = true, annihilation = :none)
    refs = GIPaper.load_reference_spectrum(joinpath(root, "data", file))
    return compare_reference(
        params,
        mq,
        refs;
        contact_hyperfine = true,
        use_fine_structure = true,
        antisymmetric_spin_orbit_mixing = aso,
        tensor_mixing = tensor,
        isoscalar_pseudoscalar_annihilation = annihilation,
        strange_mass_GeV = annihilation == :none ? nothing : mq["s"],
        ngrid = ngrid,
        rmax = rmax,
        kinetic = :relativistic,
    )
end

function group_first_pair(rows, field::Symbol, marker::String)
    selected = [row for row in rows if getproperty(row, field) == marker]
    groups = Dict{Any,Vector{Any}}()
    for row in selected
        if field == :same_j_mixing_scheme
            key = (row.sector, row.n, row.L, row.J, row.m1_GeV, row.m2_GeV)
        else
            key = (row.sector, row.J, row.m1_GeV, row.m2_GeV, round(row.tensor_offdiag_GeV; digits = 12))
        end
        push!(get!(groups, key, Any[]), row)
    end
    pairs = [sort(group; by = row -> row.reference_GeV) for group in values(groups) if length(group) == 2]
    isempty(pairs) && error("no pair found for $marker")
    return first(sort(pairs; by = pair -> sum(row.reference_GeV for row in pair)))
end

function matrix_for_pair(pair, kind::Symbol)
    if kind == :aso
        diag = [row.same_j_unmixed_GeV for row in pair]
        off = first(pair).same_j_offdiag_GeV
    elseif kind == :tensor
        diag = [row.tensor_unmixed_GeV for row in pair]
        off = first(pair).tensor_offdiag_GeV
    else
        error("unsupported pair kind $kind")
    end
    return [diag[1] off; off diag[2]]
end

function pair_effect(pair, kind::Symbol)
    before = kind == :aso ? [row.same_j_unmixed_GeV for row in pair] :
             kind == :tensor ? [row.tensor_unmixed_GeV for row in pair] :
             error("unsupported pair kind $kind")
    after = [row.predicted_GeV for row in pair]
    labels = [row_label(row) for row in pair]
    return labels, before, after
end

function isoscalar_solution(params, mq, rows; ngrid = 180, rmax = 18.0)
    qm = ConstituentMasses(mq["q"], mq["q"])
    sm = ConstituentMasses(mq["s"], mq["s"])
    q_solution = channel_solution(
        params, qm, 0;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax, kinetic = :relativistic),
    )
    s_solution = channel_solution(
        params, sm, 0;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax, kinetic = :relativistic),
    )
    ordered = sort(rows; by = row -> (row.n, row.reference_GeV))
    s_wave = radial_wave(s_solution, 1)
    s_contact = contact_hyperfine_nonperturbative_states(params, sm, "S", 1, s_wave.r, 2)
    s_levels = isnothing(s_contact) ? Float64[] : s_contact.eigenvalues_GeV
    diag = [
        ordered[1].isoscalar_annihilation_unmixed_GeV,
        s_levels[1],
        ordered[3].isoscalar_annihilation_unmixed_GeV,
        s_levels[2],
    ]
    basis = [
        # n nbar rows are the coherent (u ubar + d dbar)/sqrt(2) channel; s sbar are not.
        pseudoscalar_annihilation_basis_input(BasisState(1, "S", 1, 0; label = "1 n nbar", flavors = (:q, :q)), mq["q"], diag[1], radial_wave(q_solution, 1); isoscalar_coherent = true),
        pseudoscalar_annihilation_basis_input(BasisState(1, "S", 1, 0; label = "1 s sbar", flavors = (:s, :s)), mq["s"], diag[2], radial_wave(s_solution, 1); isoscalar_coherent = false),
        pseudoscalar_annihilation_basis_input(BasisState(2, "S", 1, 0; label = "2 n nbar", flavors = (:q, :q)), mq["q"], diag[3], radial_wave(q_solution, 2); isoscalar_coherent = true),
        pseudoscalar_annihilation_basis_input(BasisState(2, "S", 1, 0; label = "2 s sbar", flavors = (:s, :s)), mq["s"], diag[4], radial_wave(s_solution, 2); isoscalar_coherent = false),
    ]
    p1 = isoscalar_pseudoscalar_annihilation_solution(PaperP1Annihilation(), params, basis)
    return p1.block.matrix
end

# ## Plots
#
# The three figures answer different questions:
#
# 1. Which order does mixing apply in?
# 2. What matrix block is being diagonalized?
# 3. How far do representative masses move?

function draw_pipeline(path)
    fig = Figure(size = (1700, 330), backgroundcolor = RGBf(0.98, 0.98, 0.965))
    ax = Axis(fig[1, 1]; limits = (0, 15, 0, 3.4), aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    boxes = [
        (1.4, "1 fixed sectors", "solve L,S,J channels"),
        (4.4, "2 antisymmetric LS", "^1L_J <-> ^3L_J"),
        (7.4, "3 tensor", "L=J-1 <-> L=J+1"),
        (10.4, "4 annihilation", "flavor/radial block"),
        (13.4, "5 compare", "physical rows + residuals"),
    ]
    colors = [RGBf(0.22, 0.43, 0.67), RGBf(0.85, 0.38, 0.24), RGBf(0.28, 0.55, 0.39), RGBf(0.58, 0.38, 0.68), RGBf(0.26, 0.26, 0.26)]
    for (i, (x, title, subtitle)) in enumerate(boxes)
        poly!(ax, Rect2f(x - 1.12, 1.25, 2.24, 1.1); color = (colors[i], 0.18), strokecolor = colors[i], strokewidth = 3)
        text!(ax, x, 1.95; text = title, align = (:center, :center), fontsize = 24, color = colors[i], font = :bold)
        text!(ax, x, 1.55; text = subtitle, align = (:center, :center), fontsize = 16, color = RGBf(0.18, 0.18, 0.18))
        if i < length(boxes)
            lines!(ax, [x + 1.18, x + 2.04], [1.8, 1.8]; linewidth = 3, color = RGBf(0.25, 0.25, 0.25))
            scatter!(ax, [x + 2.04], [1.8]; marker = :rtriangle, markersize = 22, color = RGBf(0.25, 0.25, 0.25))
        end
    end
    text!(
        ax,
        7.5,
        0.45;
        text = "All mixing acts after fixed-sector computation, using cached radial solutions. Reference comparison is last.",
        align = (:center, :center),
        fontsize = 20,
        color = RGBf(0.15, 0.15, 0.15),
    )
    save(path, fig)
end

function draw_matrix_blocks(path, aso_matrix, tensor_matrix, ann_matrix)
    fig = Figure(size = (1500, 520), backgroundcolor = RGBf(0.98, 0.98, 0.965))
    specs = [
        ("Antisymmetric spin-orbit", aso_matrix, ["^1P_1", "^3P_1"]),
        ("Tensor", tensor_matrix, ["^3S_1", "^3D_1"]),
        ("Isoscalar annihilation P1", ann_matrix, ["1 nn", "1 ss", "2 nn", "2 ss"]),
    ]
    for (i, (title, mat, labels)) in enumerate(specs)
        ax = Axis(
            fig[1, i];
            title = title,
            xticks = (1:length(labels), labels),
            yticks = (1:length(labels), labels),
            aspect = DataAspect(),
        )
        heatmap!(ax, mat; colormap = :viridis)
        for a in axes(mat, 1), b in axes(mat, 2)
            text!(ax, b, a; text = @sprintf("%.3f", mat[a, b]), align = (:center, :center), fontsize = 16, color = :white)
        end
        xlims!(ax, 0.5, length(labels) + 0.5)
        ylims!(ax, 0.5, length(labels) + 0.5)
    end
    Label(
        fig[2, 1:3],
        "Diagonal entries are unmixed masses in GeV. Off-diagonal entries are the mechanism-specific coupling in the selected block.",
        fontsize = 20,
        tellwidth = false,
    )
    save(path, fig)
end

function draw_effects(path, aso_pair, tensor_pair, ann_rows)
    panels = [
        ("Antisymmetric spin-orbit", pair_effect(aso_pair, :aso)...),
        ("Tensor", pair_effect(tensor_pair, :tensor)...),
        (
            "Calibrated isoscalar annihilation",
            ["eta(1)", "eta'(1)", "eta(2)", "eta'(2)"],
            [row.isoscalar_annihilation_unmixed_GeV for row in sort(ann_rows; by = row -> row.reference_GeV)],
            [row.predicted_GeV for row in sort(ann_rows; by = row -> row.reference_GeV)],
        ),
    ]
    fig = Figure(size = (1500, 620), backgroundcolor = RGBf(0.98, 0.98, 0.965))
    for (i, (title, labels, before, after)) in enumerate(panels)
        ax = Axis(
            fig[1, i];
            title = title,
            ylabel = i == 1 ? "mass (GeV)" : "",
            xticks = (1:length(labels), labels),
            xticklabelrotation = 0.18,
            ygridvisible = true,
        )
        ymin = minimum(vcat(before, after)) - 0.08
        ymax = maximum(vcat(before, after)) + 0.08
        ylims!(ax, ymin, ymax)
        xlims!(ax, 0.45, length(labels) + 0.55)
        for j in eachindex(labels)
            lines!(ax, [j - 0.13, j + 0.13], [before[j], after[j]]; linewidth = 3, color = RGBf(0.2, 0.2, 0.2))
            scatter!(ax, [j - 0.13], [before[j]]; markersize = 18, color = RGBf(0.45, 0.45, 0.45), label = j == 1 && i == 1 ? "before" : nothing)
            scatter!(ax, [j + 0.13], [after[j]]; markersize = 20, color = RGBf(0.1, 0.45, 0.75), label = j == 1 && i == 1 ? "after" : nothing)
            text!(
                ax,
                j,
                max(before[j], after[j]) + 0.018;
                text = @sprintf("%+.1f MeV", 1000 * (after[j] - before[j])),
                align = (:center, :bottom),
                fontsize = 14,
                color = RGBf(0.15, 0.15, 0.15),
            )
        end
        if i == 1
            axislegend(ax; position = :lt)
        end
    end
    Label(
        fig[2, 1:3],
        "Small post-diagonalization blocks barely move heavy/open-flavor levels; the isoscalar pseudoscalar annihilation block is intentionally large.",
        fontsize = 20,
        tellwidth = false,
    )
    save(path, fig)
end

function write_markdown(path)
    open(path, "w") do io
        println(io, "# Mixing Mechanism Demo")
        println(io)
        println(io, "Generated by `julia docs/mixing_studies.jl`.")
        println(io, numerics_provenance(FiniteDifferenceSolver(ngrid = 180, rmax = 18.0)))
        println(io)
        println(io, "## Application Order")
        println(io)
        println(io, "![Mixing application order](figures/mixing_order.png)")
        println(io)
        println(io, "The fixed-sector solve comes first. `compare` then applies the assigned post-diagonalization blocks in this order: antisymmetric spin-orbit, tensor, and isoscalar annihilation. Reference residuals are computed only after those assignments.")
        println(io)
        println(io, "## Matrix Blocks")
        println(io)
        println(io, "![Mixing block matrices](figures/mixing_blocks.png)")
        println(io)
        println(io, "The first two mechanisms are ordinary two-state blocks. Isoscalar pseudoscalar annihilation is a flavor/radial block over `1 nn`, `1 ss`, `2 nn`, and `2 ss`. The comparison layer selects the mechanism; cached radial solutions provide the needed matrix elements.")
        println(io)
        println(io, "## Mass Movement")
        println(io)
        println(io, "![Mixing mass movement](figures/mixing_effects.png)")
        println(io)
        println(io, "The arrows show unmixed-to-mixed movement for representative blocks. Antisymmetric spin-orbit and tensor mixing are small corrections in these examples. The calibrated isoscalar pseudoscalar annihilation path is large by design because it absorbs the GI Table-III/Fig.-5 pseudoscalar assignment.")
        println(io)
        println(io, "## Concrete Blocks Used")
        println(io)
        println(io, "- Antisymmetric spin-orbit example: `charmed` sector, `1^1P_1` / `1^3P_1`.")
        println(io, "- Tensor example: `isovector` sector, `2^3S_1` / `1^3D_1`.")
        println(io, "- Isoscalar annihilation example: `isoscalar` `^1S_0` block over `1 nn`, `1 ss`, `2 nn`, `2 ss`.")
    end
end

# ## Execute Study
#
# The selected examples are intentionally representative rather than exhaustive:
# a heavy-light/open-flavor pair for antisymmetric spin-orbit, a light isovector
# triplet pair for tensor mixing, and the light isoscalar pseudoscalar block for
# annihilation.

function main()
    mkpath(OUTDIR)
    params, mq = load_parameters_and_quark_masses(default_parameters_path())

    charmed_rows = sector_rows(params, mq, "reference_spectrum_charmed.csv"; aso = true, tensor = false)
    aso_pair = group_first_pair(charmed_rows, :same_j_mixing_scheme, "antisymmetric_spin_orbit")
    aso_matrix = matrix_for_pair(aso_pair, :aso)

    isovector_rows = sector_rows(params, mq, "reference_spectrum_isovector.csv"; aso = false, tensor = true)
    tensor_pair = group_first_pair(isovector_rows, :tensor_mixing_scheme, "tensor_mixing")
    tensor_matrix = matrix_for_pair(tensor_pair, :tensor)

    isoscalar_rows = sector_rows(
        params,
        mq,
        "reference_spectrum_isoscalar.csv";
        aso = false,
        tensor = false,
        annihilation = :calibrated_p1,
    )
    ann_rows = [row for row in isoscalar_rows if row.isoscalar_annihilation_scheme == "calibrated_p1"]
    ann_matrix = isoscalar_solution(params, mq, ann_rows)

    draw_pipeline(joinpath(OUTDIR, "mixing_order.png"))
    draw_matrix_blocks(joinpath(OUTDIR, "mixing_blocks.png"), aso_matrix, tensor_matrix, ann_matrix)
    draw_effects(joinpath(OUTDIR, "mixing_effects.png"), aso_pair, tensor_pair, ann_rows)
    write_markdown(joinpath(root, "docs", "residual_reports", "mixing_mechanism_demo.md"))

    println("wrote ", joinpath(root, "docs", "residual_reports", "mixing_mechanism_demo.md"))
    println("wrote figures to ", OUTDIR)
end

main()
