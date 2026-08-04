### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 95d0273c-66e1-4f98-ba9b-0690edca44cb
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(@__DIR__)          # examples/Project.toml — GIModel + CairoMakie + PlutoUI
    Pkg.instantiate()

    using GIModel
    using CairoMakie
    using PlutoUI
    using LaTeXStrings
    using Printf
end

# ╔═╡ cca85821-185b-4126-9239-897d0c388fae
md"""
# Charmonium ``\to`` bottomonium

The Godfrey–Isgur medium is flavour-blind: the string tension ``b``, the Coulomb
strength, the smearing width ``\sigma_0`` and the running ``\alpha_s(Q^2)`` see
nothing but the **constituent masses**. A heavy quark is therefore a free dial —
turn ``m_Q`` from the charm value up to the bottom value and the whole charmonium
spectrum walks over into the bottomonium one, with no reparametrization anywhere.

Drag the slider and watch three things:

* the levels **climb by ``\approx 2\,\Delta m_Q``** (or ``\Delta m_Q`` for heavy-light) — that is just the rest mass;
* the **splittings shrink**: hyperfine ``\propto 1/m_Q``, so ``J/\psi-\eta_c`` collapses towards ``\Upsilon-\eta_b``;
* the **pattern of the excitations barely moves** once the quark masses are subtracted — the confining potential is the same string in both systems.
"""

# ╔═╡ 23cc7a63-7ad5-4c1a-b3ad-9ea61a5b31a4
md"""
## Controls

| | |
|---:|:---|
| **heavy-quark mass** ``m_Q`` | $(@bind mQ PlutoUI.Slider(1.20:0.02:5.40; default = 1.628, show_value = true)) GeV |
| **flavour content** | $(@bind pairing PlutoUI.Select([:QQ => "Q Q̄  —  quarkonium (both quarks heavy)", :Qq => "Q q̄  —  heavy-light (partner is u/d)", :Qs => "Q s̄  —  heavy-strange"])) |
| **mass axis** | $(@bind lock_axis PlutoUI.CheckBox(default = true)) locked to the full sweep (uncheck to zoom on the current spectrum) |
| **subtract ``m_1 + m_2``** | $(@bind subtract PlutoUI.CheckBox(default = false)) plot the binding part ``M - m_1 - m_2`` instead of ``M`` |

`Slider` is qualified because Makie exports one too.

The physical values are ``m_c = 1.628`` GeV and ``m_b = 4.977`` GeV (dashed markers in the right-hand panels).
"""

# ╔═╡ e63b63b6-c494-4825-88f9-6010465d8c8e
md"""
!!! tip "What the switch changes"
	``Q\bar Q`` moves **both** constituent masses, so the levels climb twice as
	fast and the same-``J`` antisymmetric spin-orbit mixing is switched off by
	charge conjugation (`is_equal_flavor`). ``Q\bar q`` and ``Q\bar s`` keep the
	partner fixed: only the heavy quark gets heavier, the light quark still
	orbits at the same scale, and the unequal masses turn the ``{}^1P_1``/``{}^3P_1``
	mixing back on. With the quark masses subtracted, the heavy-light curves
	flatten out onto ``\bar\Lambda = M - m_Q`` — heavy-quark symmetry appearing
	on screen.
"""

# ╔═╡ b4646a9d-96c9-486e-99e7-593cb99f02dd
md"""
## The same spectrum in ``J^P``

``{}^{2S+1}L_J`` is model bookkeeping; ``J^P`` is what an experiment sees. Below,
each ``L``-block is collapsed onto its ``J^P`` columns — ``{}^1P_1`` and
``{}^3P_1`` end up in the same ``1^+`` column, drawn as two half-bars — and the
mass axis simply follows the current spectrum, so the ground state and the top
level stay put and only the pattern between them moves.

The two panels answer the same question in two ways: what does the *shape* of
the spectrum do as the quark gets heavy?
"""

# ╔═╡ d4b35dbb-0297-421a-9c79-e0d6e5fe0c6d
md"""
!!! note "Why the 1⁺ pair is drawn side by side"
	Those two states are exactly what the antisymmetric spin-orbit term mixes —
	but only when the quarks have different masses. For ``Q\bar Q`` charge
	conjugation keeps them apart (``{}^1P_1`` is ``1^{+-}``, ``{}^3P_1`` is
	``1^{++}``) and the model switches the mixing off; flip the flavour switch to
	``Q\bar q`` and the pair starts to repel. Note also that ``1^-`` occurs twice,
	once in the ``S``-block and once in the ``D``-block: those are the two states
	the tensor term mixes (``2{}^3S_1`` with ``1{}^3D_1``).
"""

# ╔═╡ 71dbbdf3-984b-4db3-b060-a967d4f6699b
md"""
## Record the animation

Ticking the box sweeps ``m_Q`` from ``m_c`` to ``m_b`` and back and writes a GIF
next to this notebook (`examples/figures/`). It re-solves the model once per
frame, so it takes about a minute.
"""

# ╔═╡ 4b603125-af67-4329-bf5b-b592efe21db8
md"""**render the GIF** $(@bind make_movie PlutoUI.CheckBox(default = false)) of $(@bind movie_kind PlutoUI.Select([:levels => "the level scheme", :shape => "the J^P shape"]))"""

# ╔═╡ a6baa3dc-868f-4bc0-9d6f-fff6432749e2
md"""
---
## Machinery

Everything below is the plumbing: model setup, the mass scan that feeds the
right-hand panels, and the figure builder.
"""

# ╔═╡ 9ed0f5e7-7753-426a-aee5-1620a6bf43ab
md"""### Model setup"""

# ╔═╡ 04b38bd7-3852-4b6b-b37a-3b0bb072c7b3
params, quark_masses = load_parameters_and_quark_masses(
    joinpath(pkgdir(GIModel), "data", "parameters.provisional.toml"),
);

# ╔═╡ ac4a3e21-7e1c-47b3-a1b9-8adeb730fc93
# Coarser than the defaults (450 / 24.0) and indistinguishable from them here:
# masses agree to well under an MeV, and one solve drops to ~0.4 s, which is what
# makes the slider usable. Heavy quarkonium is compact, so it needs points near
# the origin rather than reach.
solver = RadialSolver(; ngrid = 360, rmax = 20.0, nlevels_per_channel = 3)

# ╔═╡ a7f2ee1f-1687-4a0f-8bec-21ebbde5c613
levels = spectrum_levels(2; L_labels = ("S", "P", "D"))

# ╔═╡ 5b4da1bc-7fc3-418c-955b-38338a1555d9
m_charm, m_bottom = quark_masses["c"], quark_masses["b"]

# ╔═╡ ac12e568-2a86-41b8-8dcd-b9ab051e1aa9
begin
    # The heavy quark is the dial; the partner is either a copy of it (quarkonium)
    # or a fixed light/strange quark (heavy-light).
    partner_quark(p) =
        p === :Qq ? LightQuark(quark_masses["q"]) :
        p === :Qs ? StrangeQuark(quark_masses["s"]) : nothing

    function meson_at(p, mass_GeV)
        Q = HeavyQuark{:up}(mass_GeV, :Q)
        partner = partner_quark(p)
        return isnothing(partner) ? Meson(Q, Q) : Meson(Q, partner)
    end

    spectrum_at(p, mass_GeV) =
        compute_spectrum(params, meson_at(p, mass_GeV); levels = levels, solver = solver)
end

# ╔═╡ 6948f467-4c25-4c84-bda1-dda84dc0070e
spec = spectrum_at(pairing, mQ)

# ╔═╡ 6701758f-9916-4b9d-b46b-300da03b2b31
spec

# ╔═╡ c83f7247-2e3e-47ab-9ae9-5d8895269dfe
md"""### Scan over the heavy-quark mass"""

# ╔═╡ 5bf97951-0a3b-4127-8542-be6cca806f0f
mQ_grid = range(1.2, 5.4; length = 22)

# ╔═╡ 2bc03fbb-f73c-4abb-bf46-eca8db565837
"""
    trajectories(pairing, grid)

Every level of `levels` solved on the whole `grid` of heavy-quark masses. Rows
follow `levels`, columns follow `grid`; `total` is `m_1 + m_2` per column and
`hf` the ``1{}^3S_1 - 1{}^1S_0`` splitting in MeV.
"""
function trajectories(p, grid)
    M = Matrix{Float64}(undef, length(levels), length(grid))
    total = Vector{Float64}(undef, length(grid))
    for (j, mass_GeV) in enumerate(grid)
        s = spectrum_at(p, mass_GeV)
        for (i, state) in enumerate(s.states)
            M[i, j] = state.mass_GeV
        end
        total[j] = s.meson.constituent_masses.m1_GeV + s.meson.constituent_masses.m2_GeV
    end
    row(label) = findfirst(l -> l.label == label, levels)
    hf = 1e3 .* (M[row("1^3S_1"), :] .- M[row("1^1S_0"), :])
    return (grid = collect(grid), mass = M, total = total, hf = hf)
end

# ╔═╡ 57b97451-30ec-4a96-8d0d-5368b65fdaf5
# A scan is ~20 solves, so keep the finished ones: flipping the flavour switch
# back and forth then costs nothing. Re-run this cell to clear the cache after
# editing `trajectories` or the solver.
traj_cache = Dict{Any,Any}()

# ╔═╡ aaeace5c-3aac-49a6-bd7e-0abff975d92e
traj = get!(() -> trajectories(pairing, mQ_grid), traj_cache, (pairing, mQ_grid))

# ╔═╡ b0942b1c-bbd6-4c1b-b507-c5f63f21c9e7
md"""### Figure"""

# ╔═╡ 00a2e3af-31a1-4d9f-8880-a944670d2adb
begin
    L_ORDER = Dict("S" => 0, "P" => 1, "D" => 2)
    L_COLOR = Dict("S" => "#0072B2", "P" => "#D55E00", "D" => "#009E73")  # Okabe-Ito
    NOW_COLOR = "#CC79A7"

    gi_theme = merge(
        theme_latexfonts(),
        Theme(
            fontsize = 16,
            Axis = (
                xticksmirrored = true, yticksmirrored = true,
                xtickalign = 1, ytickalign = 1,
                topspinevisible = true, rightspinevisible = true,
            ),
        ),
    )

    L_LETTER = ("S", "P", "D")

    pairing_tex(p) =
        p === :QQ ? "Q\\bar{Q}\\;\\mathrm{quarkonium}" :
        p === :Qq ? "Q\\bar{q}\\;\\mathrm{heavy-light}" :
        "Q\\bar{s}\\;\\mathrm{heavy-strange}"
    column_tex(k) = latexstring("{}^{$(k[2])}\\mathrm{$(L_LETTER[k[1]+1])}_{$(k[3])}")
    level_title(mass_GeV) =
        latexstring("\\mathrm{levels~at~} m_Q = ", @sprintf("%.2f", mass_GeV), "\\,\\mathrm{GeV}")

    # Parity of a q qbar level is (-1)^(L+1); within one L-block it is common to
    # all levels, so a J^P column is just "the J's of this block".
    parity_char(L) = iseven(L + 1) ? "+" : "-"
    jp_tex(k) = latexstring(k[2], "^{", parity_char(k[1]), "}")
    state_tex(l) = latexstring(l.n, "\\,^{", l.multiplicity, "}\\mathrm{", l.L_label, "}_{", l.J, "}")

    # Frame that hugs the current spectrum, with headroom on top for the
    # block labels. Extra margin above than below, hence not `extrema` alone.
    function spectrum_span(ys)
        lo, hi = extrema(ys)
        span = hi - lo
        return (lo - 0.05 * span, hi + 0.15 * span)
    end
end

# ╔═╡ 2f63bd7b-ec01-4103-93ec-9c98a1632eb8
"""
    snapshot(pairing, mQ, subtract)

One solve, reduced to what the moving parts of the figure need: the level
positions (mass, or mass minus the constituent masses) and the 1S hyperfine
splitting in MeV.
"""
function snapshot(p, mass_GeV, subtract)
    s = spectrum_at(p, mass_GeV)
    cm = s.meson.constituent_masses
    shift = subtract ? cm.m1_GeV + cm.m2_GeV : 0.0
    ys = [state.mass_GeV - shift for state in s.states]
    hf = 1e3 * (spectrum_state(s, "1^3S_1").mass_GeV - spectrum_state(s, "1^1S_0").mass_GeV)
    return (ys = ys, hf = hf)
end

# ╔═╡ 75afdc5e-d4a4-4c45-b2a7-6478774c078f
"""
    transition_figure(pairing, mQ::Observable, traj; subtract, lock_axis, size)

Left: the level scheme at the current ``m_Q``, one column per ``{}^{2S+1}L_J``,
bars stacked by radial quantum number. Right: the same levels as functions of
``m_Q`` with the current value marked, and the 1S hyperfine splitting below.

`mQ` is an `Observable` so the very same figure serves the slider (one static
value) and `record_sweep` (a value that moves) — each update costs one solve.
"""
function transition_figure(p, mQ_obs::Observable, traj;
    subtract = false, lock_axis = true, size = (1220, 760))

    Y = subtract ? traj.mass .- traj.total' : traj.mass
    ylo, yhi = extrema(Y)
    pad = 0.06 * (yhi - ylo)
    ylo, yhi = ylo - pad, yhi + pad
    ylab = subtract ? L"M - m_1 - m_2\;\;[\mathrm{GeV}]" : L"M\;\;[\mathrm{GeV}]"

    snap = @lift(snapshot(p, $mQ_obs, subtract))
    columns = sort(unique((L_ORDER[l.L_label], l.multiplicity, l.J) for l in levels))
    xof = Dict(k => i for (i, k) in enumerate(columns))

    with_theme(gi_theme) do
        fig = Figure(; size = size)
        Label(fig[0, 1:2],
            latexstring("\\mathrm{Godfrey-Isgur~spectrum:}\\;", pairing_tex(p));
            fontsize = 21, padding = (0, 0, 0, 4))

        # --- level scheme ----------------------------------------------------
        ax1 = Axis(fig[1:2, 1];
            title = @lift(level_title($mQ_obs)),
            xticks = (1:length(columns), column_tex.(columns)),
            ylabel = ylab, xgridvisible = false,
            yminorticksvisible = true, yminorgridvisible = true,
            yautolimitmargin = (0.06, 0.16))
        lock_axis && ylims!(ax1, ylo, yhi)
        xlims!(ax1, 0.4, length(columns) + 0.6)

        xspan = length(columns) + 0.2   # width of the xlims! window above
        for L in ("S", "P", "D")
            xs = [xof[k] for k in columns if k[1] == L_ORDER[L]]
            isempty(xs) && continue
            vspan!(ax1, minimum(xs) - 0.5, maximum(xs) + 0.5; color = (L_COLOR[L], 0.07))
            text!(ax1, Point2f(((minimum(xs) + maximum(xs)) / 2 - 0.4) / xspan, 0.985);
                text = latexstring("\\mathrm{$L}\\text{-wave}"), space = :relative,
                align = (:center, :top), fontsize = 15, color = L_COLOR[L])
        end
        for (i, l) in enumerate(levels)
            x = xof[(L_ORDER[l.L_label], l.multiplicity, l.J)]
            color = L_COLOR[l.L_label]
            style = l.n == 1 ? :solid : (:dot, :dense)
            lines!(ax1, @lift([Point2f(x - 0.38, $snap.ys[i]), Point2f(x + 0.38, $snap.ys[i])]);
                color = color, linewidth = 4, linestyle = style)
            text!(ax1, @lift(Point2f(x - 0.38, $snap.ys[i])); text = string(l.n),
                align = (:left, :bottom), offset = (2, 3), fontsize = 12, color = color)
        end

        # --- level trajectories ----------------------------------------------
        ax2 = Axis(fig[1, 2];
            title = L"\mathrm{charm}\;\longrightarrow\;\mathrm{bottom}",
            xlabel = L"m_Q\;\;[\mathrm{GeV}]", ylabel = ylab, yminorticksvisible = true)
        ylims!(ax2, ylo, yhi)
        xlims!(ax2, first(traj.grid), last(traj.grid))
        for (mark, name) in ((m_charm, "c"), (m_bottom, "b"))
            vlines!(ax2, mark; color = (:gray30, 0.55), linestyle = :dash, linewidth = 1)
            text!(ax2, mark, yhi; text = latexstring("m_$name"), align = (:left, :top),
                offset = (3, -3), fontsize = 14, color = :gray30)
        end
        for (i, l) in enumerate(levels)
            lines!(ax2, traj.grid, Y[i, :]; color = (L_COLOR[l.L_label], 0.85),
                linewidth = 1.6, linestyle = l.n == 1 ? :solid : (:dot, :dense))
        end
        vlines!(ax2, @lift([$mQ_obs]); color = NOW_COLOR, linewidth = 2)
        scatter!(ax2, @lift(Point2f.($mQ_obs, $snap.ys)); color = NOW_COLOR,
            markersize = 7, strokewidth = 0.5, strokecolor = :white)

        # --- hyperfine splitting ---------------------------------------------
        ax3 = Axis(fig[2, 2];
            title = L"1S\;\mathrm{hyperfine~splitting}",
            xlabel = L"m_Q\;\;[\mathrm{GeV}]",
            ylabel = L"M(1^3S_1) - M(1^1S_0)\;\;[\mathrm{MeV}]")
        xlims!(ax3, first(traj.grid), last(traj.grid))
        for mark in (m_charm, m_bottom)
            vlines!(ax3, mark; color = (:gray30, 0.55), linestyle = :dash, linewidth = 1)
        end
        lines!(ax3, traj.grid, traj.hf; color = L_COLOR["S"], linewidth = 2.5)
        scatter!(ax3, @lift([Point2f($mQ_obs, $snap.hf)]); color = NOW_COLOR,
            markersize = 11, strokewidth = 0.5, strokecolor = :white)

        Legend(fig[3, 1:2],
            [LineElement(color = :gray25, linewidth = 3),
                LineElement(color = :gray25, linewidth = 3, linestyle = (:dot, :dense)),
                LineElement(color = NOW_COLOR, linewidth = 3)],
            [L"n = 1\;\mathrm{(ground~radial)}", L"n = 2\;\mathrm{(first~radial~excitation)}",
                L"\mathrm{current}\;m_Q"];
            orientation = :horizontal, framevisible = false,
            patchsize = (28, 10), colgap = 24, tellheight = true)

        colsize!(fig.layout, 1, Relative(0.55))
        rowsize!(fig.layout, 2, Relative(0.34))
        fig
    end
end

# ╔═╡ d54530a1-76d0-482e-ad3b-168a4be74e45
transition_figure(pairing, Observable(mQ), traj; subtract, lock_axis)

# ╔═╡ 736d4588-ce0c-4573-bafc-d1398ecf7a53
"""
    shape_figure(pairing, mQ::Observable, traj; subtract, size)

The same levels seen through the quantum numbers that survive: within each
``L``-block the columns are ``J^P``, so ``{}^1P_1`` and ``{}^3P_1`` share the
``1^+`` column (side by side — for quarkonium they nearly coincide), and the
mass axis is rescaled to the current spectrum. Ground state and top level are
then pinned to the frame and only the *pattern* between them can move.

Right panel: the same rescaling as a function of ``m_Q``, so the bottom and top
curves are flat at 0 and 1 by construction.
"""
function shape_figure(p, mQ_obs::Observable, traj; subtract = false, size = (1220, 620))
    snap = @lift(snapshot(p, $mQ_obs, subtract))
    ylab = subtract ? L"M - m_1 - m_2\;\;[\mathrm{GeV}]" : L"M\;\;[\mathrm{GeV}]"

    columns = sort(unique((L_ORDER[l.L_label], l.J) for l in levels))
    xof = Dict(k => i for (i, k) in enumerate(columns))
    # A column holding both a singlet and a triplet holds the pair the same-J
    # spin-orbit term mixes; give them half-columns instead of stacking them at
    # the same x, where quarkonium would draw them on top of each other.
    split = Dict(k => length(unique(l.multiplicity for l in levels
        if (L_ORDER[l.L_label], l.J) == k)) > 1 for k in columns)
    function bar_span(l)
        x = xof[(L_ORDER[l.L_label], l.J)]
        split[(L_ORDER[l.L_label], l.J)] || return (x - 0.42, x + 0.42)
        return l.multiplicity == 1 ? (x - 0.44, x - 0.04) : (x + 0.04, x + 0.44)
    end

    # Level positions rescaled onto the span of the spectrum they belong to.
    Ynorm = similar(traj.mass)
    for j in axes(traj.mass, 2)
        lo, hi = extrema(view(traj.mass, :, j))
        Ynorm[:, j] .= (view(traj.mass, :, j) .- lo) ./ (hi - lo)
    end

    with_theme(gi_theme) do
        fig = Figure(; size = size)
        Label(fig[0, 1:2],
            latexstring("\\mathrm{Spectrum~shape~in~} J^P\\mathrm{:}\\;", pairing_tex(p));
            fontsize = 21, padding = (0, 0, 0, 4))

        # --- J^P level scheme, frame following the spectrum ------------------
        ax1 = Axis(fig[1, 1];
            title = @lift(level_title($mQ_obs)),
            xticks = (1:length(columns), jp_tex.(columns)),
            xlabel = L"J^P", ylabel = ylab, xgridvisible = false,
            yminorticksvisible = true, yminorgridvisible = true)
        xlims!(ax1, 0.4, length(columns) + 0.6)
        span = @lift(spectrum_span($snap.ys))
        ylims!(ax1, span[]...)
        on(lim -> ylims!(ax1, lim[1], lim[2]), span)

        xspan = length(columns) + 0.2
        for (Lval, L) in enumerate(L_LETTER)
            xs = [xof[k] for k in columns if k[1] == Lval - 1]
            isempty(xs) && continue
            vspan!(ax1, minimum(xs) - 0.5, maximum(xs) + 0.5; color = (L_COLOR[L], 0.07))
            text!(ax1, Point2f(((minimum(xs) + maximum(xs)) / 2 - 0.4) / xspan, 0.985);
                text = latexstring("\\mathrm{$L}\\text{-wave}"), space = :relative,
                align = (:center, :top), fontsize = 15, color = L_COLOR[L])
        end
        for (i, l) in enumerate(levels)
            xlo, xhi = bar_span(l)
            color = L_COLOR[l.L_label]
            style = l.n == 1 ? :solid : (:dot, :dense)
            lines!(ax1, @lift([Point2f(xlo, $snap.ys[i]), Point2f(xhi, $snap.ys[i])]);
                color = color, linewidth = 4, linestyle = style)
            text!(ax1, @lift(Point2f((xlo + xhi) / 2, $snap.ys[i])); text = state_tex(l),
                align = (:center, :bottom), offset = (0, 3), fontsize = 11, color = color)
        end

        # --- the same rescaling across the sweep -----------------------------
        ax2 = Axis(fig[1, 2];
            title = L"\mathrm{charm}\;\longrightarrow\;\mathrm{bottom}",
            xlabel = L"m_Q\;\;[\mathrm{GeV}]",
            ylabel = L"(M - M_\mathrm{min})\,/\,(M_\mathrm{max} - M_\mathrm{min})",
            yminorticksvisible = true)
        xlims!(ax2, first(traj.grid), last(traj.grid))
        ylims!(ax2, -0.05, 1.14)
        for (mark, name) in ((m_charm, "c"), (m_bottom, "b"))
            vlines!(ax2, mark; color = (:gray30, 0.55), linestyle = :dash, linewidth = 1)
            text!(ax2, mark, 1.14; text = latexstring("m_$name"), align = (:left, :top),
                offset = (3, -3), fontsize = 14, color = :gray30)
        end
        for (i, l) in enumerate(levels)
            lines!(ax2, traj.grid, Ynorm[i, :]; color = (L_COLOR[l.L_label], 0.85),
                linewidth = 1.6, linestyle = l.n == 1 ? :solid : (:dot, :dense))
        end
        vlines!(ax2, @lift([$mQ_obs]); color = NOW_COLOR, linewidth = 2)
        scatter!(ax2, @lift(let ys = $snap.ys, (lo, hi) = extrema(ys)
            Point2f.($mQ_obs, (ys .- lo) ./ (hi - lo))
        end); color = NOW_COLOR, markersize = 7, strokewidth = 0.5, strokecolor = :white)

        Legend(fig[2, 1:2],
            [LineElement(color = :gray25, linewidth = 3),
                LineElement(color = :gray25, linewidth = 3, linestyle = (:dot, :dense)),
                LineElement(color = NOW_COLOR, linewidth = 3)],
            [L"n = 1", L"n = 2", L"\mathrm{current}\;m_Q"];
            orientation = :horizontal, framevisible = false,
            patchsize = (28, 10), colgap = 24, tellheight = true)

        colsize!(fig.layout, 1, Relative(0.58))
        fig
    end
end

# ╔═╡ 5feec0bb-db26-434f-a357-708173793b2e
shape_figure(pairing, Observable(mQ), traj; subtract)

# ╔═╡ f1a0682f-2a0a-4563-96cc-8b9fe80ebd93
"""
    record_sweep(builder, pairing, traj; suffix, nframes, framerate, kwargs...) -> path

Drive the `Observable` of a figure `builder` (`transition_figure` or
[`shape_figure`](@ref)) from ``m_c`` to ``m_b`` and back, one model solve per
frame, into `examples/figures/*.gif`. Extra keywords go to the builder.
"""
function record_sweep(builder, p, traj; suffix = "", nframes = 48, framerate = 12, kwargs...)
    dir = mkpath(joinpath(@__DIR__, "figures"))
    path = joinpath(dir, "charm_to_bottom_$(p)$(suffix).gif")
    mQ_obs = Observable(m_charm)
    fig = builder(p, mQ_obs, traj; size = (980, 620), kwargs...)
    sweep = vcat(
        range(m_charm, m_bottom; length = nframes), fill(m_bottom, 8),
        range(m_bottom, m_charm; length = nframes), fill(m_charm, 8),
    )
    record(fig, path, sweep; framerate = framerate) do mass_GeV
        mQ_obs[] = mass_GeV
    end
    return path
end

# ╔═╡ 4289bec0-6ba8-4673-b3a4-0b83c2605748
if make_movie
    let path = movie_kind === :shape ?
               record_sweep(shape_figure, pairing, traj;
            suffix = "_shape", subtract = subtract) :
               record_sweep(transition_figure, pairing, traj;
            suffix = subtract ? "_binding" : "", subtract = subtract, lock_axis = true)
        LocalResource(path)
    end
else
    md"_(unticked — nothing rendered)_"
end

# ╔═╡ a4a6a37b-82ae-44fb-9453-64f25d93b3a8
TableOfContents(; title = "Contents")

# ╔═╡ Cell order:
# ╟─cca85821-185b-4126-9239-897d0c388fae
# ╠═95d0273c-66e1-4f98-ba9b-0690edca44cb
# ╟─d54530a1-76d0-482e-ad3b-168a4be74e45
# ╟─e63b63b6-c494-4825-88f9-6010465d8c8e
# ╟─b4646a9d-96c9-486e-99e7-593cb99f02dd
# ╟─5feec0bb-db26-434f-a357-708173793b2e
# ╟─23cc7a63-7ad5-4c1a-b3ad-9ea61a5b31a4
# ╟─d4b35dbb-0297-421a-9c79-e0d6e5fe0c6d
# ╠═6701758f-9916-4b9d-b46b-300da03b2b31
# ╟─71dbbdf3-984b-4db3-b060-a967d4f6699b
# ╟─4b603125-af67-4329-bf5b-b592efe21db8
# ╟─4289bec0-6ba8-4673-b3a4-0b83c2605748
# ╟─a6baa3dc-868f-4bc0-9d6f-fff6432749e2
# ╟─9ed0f5e7-7753-426a-aee5-1620a6bf43ab
# ╠═04b38bd7-3852-4b6b-b37a-3b0bb072c7b3
# ╠═ac4a3e21-7e1c-47b3-a1b9-8adeb730fc93
# ╠═a7f2ee1f-1687-4a0f-8bec-21ebbde5c613
# ╠═5b4da1bc-7fc3-418c-955b-38338a1555d9
# ╠═ac12e568-2a86-41b8-8dcd-b9ab051e1aa9
# ╠═6948f467-4c25-4c84-bda1-dda84dc0070e
# ╟─c83f7247-2e3e-47ab-9ae9-5d8895269dfe
# ╠═5bf97951-0a3b-4127-8542-be6cca806f0f
# ╠═2bc03fbb-f73c-4abb-bf46-eca8db565837
# ╠═57b97451-30ec-4a96-8d0d-5368b65fdaf5
# ╠═aaeace5c-3aac-49a6-bd7e-0abff975d92e
# ╟─b0942b1c-bbd6-4c1b-b507-c5f63f21c9e7
# ╠═00a2e3af-31a1-4d9f-8880-a944670d2adb
# ╠═2f63bd7b-ec01-4103-93ec-9c98a1632eb8
# ╠═75afdc5e-d4a4-4c45-b2a7-6478774c078f
# ╠═736d4588-ce0c-4573-bafc-d1398ecf7a53
# ╠═f1a0682f-2a0a-4563-96cc-8b9fe80ebd93
# ╟─a4a6a37b-82ae-44fb-9453-64f25d93b3a8
