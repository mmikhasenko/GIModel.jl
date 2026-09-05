#!/usr/bin/env julia
# Demonstration plots contrasting the finite-difference (FD) and
# harmonic-oscillator (HO) radial wavefunctions used by the two solver paths.
# FD returns mesh samples; HO returns a native analytic oscillator wave. This
# plotting script explicitly samples the latter on the FD grid only to render
# directly comparable curves. The calculation itself does not use that mesh.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie
using LaTeXStrings
using Printf

root = dirname(@__DIR__)
using GIModel

const OUTDIR = joinpath(@__DIR__, "spectrum_plots")
const HBARC = 0.1973269804  # GeV*fm: r[GeV^-1] -> r[fm]

const WONG2 = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
]
const FD_COLOR = "#0072B2"   # wong2 blue
const HO_COLOR = "#D55E00"   # wong2 vermillion

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

"""Unit-normalize the reduced wave on the mesh and fix the inner-lobe sign."""
function normalize_wave(u::AbstractVector, h::Real)
    norm = sqrt(sum(abs2, u) * h)
    v = norm > 0 ? u ./ norm : copy(u)
    cutoff = 0.1 * maximum(abs, v)
    i = findfirst(x -> abs(x) > cutoff, v)
    i === nothing && return v
    return v .* sign(v[i])
end

"""FD and HO reduced waves for one channel/level on the shared mesh, in fm.

`nlevels` matches the audit path (`compute_spectrum` uses 6): the HO variational
beta is chosen to minimize the highest requested level, so larger `nlevels`
degrades the low-level HO representation while FD is unaffected."""
function basis_waves(params, masses, L_label, level; nlevels = 6)
    Lval = GIModel.L_SYMBOLS[L_label]
    nl = max(nlevels, level)
    fd_solution = channel_solution(
        params, masses, Lval; solver = FiniteDifferenceSolver(), nlevels = nl)
    ho_solution = channel_solution(
        params, masses, Lval; solver = OscillatorSolver(), nlevels = nl)
    fd_wave = radial_wave(fd_solution, level)
    r = fd_wave.r
    ho_wave = sample_wave(radial_wave(ho_solution, level), r)
    h = r[2] - r[1]
    r_fm = r .* HBARC
    u_fd = normalize_wave(fd_wave.u, h)
    u_ho = normalize_wave(ho_wave.u, h)
    return (r_fm = r_fm, u_fd = u_fd, u_ho = u_ho,
        E_fd = fd_solution.eigenvalues_GeV[level], E_ho = ho_solution.eigenvalues_GeV[level])
end

# Charm and strange constituent masses (GeV), Table II.
const M_CC = ConstituentMasses(1.628, 1.628)
const M_STRANGE = ConstituentMasses(0.220, 0.419)

mkpath(OUTDIR)

# --- Figure 1: representation. Ground-state cc 1^3S_1: FD as mesh samples,
# HO as a smooth finite-basis curve. Both approximate the same function. ----
let
    w = basis_waves(GIModel.load_parameters(joinpath(dirname(root), "data", "parameters.provisional.toml")),
        M_CC, "S", 1)
    fig = Figure(size = (760, 520))
    ax = Axis(fig[1, 1];
        xlabel = L"r~\mathrm{(fm)}",
        ylabel = L"u(r)~\mathrm{(arb.)}",
        title = "",
    )
    xlims!(ax, 0, 1.6)
    # FD: line plus subsampled markers to convey the discrete mesh.
    lines!(ax, w.r_fm, w.u_fd; color = (FD_COLOR, 0.9), linewidth = 2, label = "FD (mesh)")
    idx = 1:8:length(w.r_fm)
    scatter!(ax, w.r_fm[idx], w.u_fd[idx]; color = FD_COLOR, markersize = 7)
    # HO: smooth analytic finite-basis reconstruction.
    lines!(ax, w.r_fm, w.u_ho; color = (HO_COLOR, 0.9), linewidth = 2,
        linestyle = :dash, label = "HO (basis)")
    axislegend(ax; position = :rt)
    save(joinpath(OUTDIR, "wavefn_representation.png"), fig)
    println("wrote wavefn_representation.png")
end

# --- Figure 2: (a) cc S-wave radial excitations; (b) an HO breakdown case. ---
let
    params = GIModel.load_parameters(joinpath(dirname(root), "data", "parameters.provisional.toml"))
    fig = Figure(size = (1180, 480))

    axa = Axis(fig[1, 1];
        xlabel = L"r~\mathrm{(fm)}", ylabel = L"u(r)~\mathrm{(arb.)}",
        title = "charmonium S-wave radial excitations")
    xlims!(axa, 0, 2.6)
    for (k, n) in enumerate((1, 2, 3))
        w = basis_waves(params, M_CC, "S", n)
        c = WONG2[k]
        lines!(axa, w.r_fm, w.u_fd; color = c, linewidth = 2)
        lines!(axa, w.r_fm, w.u_ho; color = c, linewidth = 2, linestyle = :dash)
    end
    # Legend proxies (solid = FD, dashed = HO; colors = n).
    lfd = [LineElement(color = :black, linewidth = 2)]
    lho = [LineElement(color = :black, linewidth = 2, linestyle = :dash)]
    ncol = [LineElement(color = WONG2[k], linewidth = 3) for k = 1:3]
    Legend(fig[1, 1], [lfd, lho, ncol[1], ncol[2], ncol[3]],
        ["FD", "HO", "1S", "2S", "3S"];
        tellheight = false, tellwidth = false,
        halign = :right, valign = :top, margin = (6, 6, 6, 6), nbanks = 1, framevisible = true)

    axb = Axis(fig[1, 2];
        xlabel = L"r~\mathrm{(fm)}", ylabel = L"\log_{10}|u(r)|",
        title = "large-r tail (log scale)")
    xlims!(axb, 0, 3.2)
    ylims!(axb, -8, 0.3)
    logabs(u) = log10.(abs.(u) .+ 1e-12)
    for (k, n) in enumerate((1, 2))
        w = basis_waves(params, M_CC, "S", n)
        c = WONG2[k]
        lines!(axb, w.r_fm, logabs(w.u_fd); color = c, linewidth = 2)
        lines!(axb, w.r_fm, logabs(w.u_ho); color = c, linewidth = 2, linestyle = :dash)
    end
    Legend(fig[1, 2], [lfd, lho, ncol[1], ncol[2]], ["FD", "HO", "1S", "2S"];
        tellheight = false, tellwidth = false,
        halign = :right, valign = :top, margin = (6, 6, 6, 6), framevisible = true)

    save(joinpath(OUTDIR, "wavefn_basis_difference.png"), fig)
    println("wrote wavefn_basis_difference.png")
end
