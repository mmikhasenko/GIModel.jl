#!/usr/bin/env julia
# Visualize the normalization bug (now fixed) in the off-diagonal same-J tensor
# mixing element, using the strange 3S1-3D1 pair. The diagonal energies/shifts
# are normalization-invariant (Rayleigh quotient) and agree between bases; the
# off-diagonal cross element dot(u_S, B K B u_D) is NOT normalization-invariant,
# so the physically-normalized HO reconstruction (||u||2 = 1/sqrt(h)) inflated it
# by ~1/h relative to the Euclidean-normalized FD eigenvectors. Dividing by the
# input norms restores agreement.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie
using LaTeXStrings
using LinearAlgebra
using Printf

root = dirname(@__DIR__)
using GIModel
const G = GIModel

const OUTDIR = joinpath(@__DIR__, "spectrum_plots")
const HBARC = 0.1973269804
const FD_COLOR = "#0072B2"   # wong2 blue
const HO_COLOR = "#D55E00"   # wong2 vermillion
const FIX_COLOR = "#009E73"  # wong2 green

set_theme!(merge(theme_latexfonts(), Theme(fontsize = 16, Axis = (
    xticksmirrored = true, yticksmirrored = true, xtickalign = 1, ytickalign = 1,
    xminortickalign = 1, yminortickalign = 1, topspinevisible = true,
    rightspinevisible = true, xgridvisible = true, ygridvisible = true))))

params = G.load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
ho = G.with_basis(params, HarmonicOscillatorBasis)
m = ConstituentMasses(0.220, 0.419)   # strange q-sbar

fixsign(v) = v .* sign(v[argmax(abs.(v))])

"""Return the off-diagonal 3S1-3D1 tensor integrand pieces on the mesh."""
function pieces(P)
    _, vS, r = channel_solution(P, m, 0; nlevels = 6)
    _, vD, _ = channel_solution(P, m, 2; nlevels = 6)
    uS = fixsign(vS[:, 1]); uD = fixsign(vD[:, 1])
    h = r[2] - r[1]
    side = G.gi_spin_dependent_side_exponent(params.epsilon_t)
    BS = G.momentum_relativization_matrix(m.m1_GeV, m.m2_GeV, side, eigen(G.p2_operator(params, m.m1_GeV, 0, r, h)))
    BD = G.momentum_relativization_matrix(m.m1_GeV, m.m2_GeV, side, eigen(G.p2_operator(params, m.m1_GeV, 2, r, h)))
    K = [G.tensor_kernel_smeared_coulomb(params, m, ri) for ri in r]
    integ_raw = (BS * uS) .* K .* (BD * uD)              # solver's bare cross product
    normfac = sqrt(dot(uS, uS)) * sqrt(dot(uD, uD))      # ||u_S||2 ||u_D||2
    (r_fm = r .* HBARC, uD = uD ./ sqrt(dot(uD, uD)), integ_raw = integ_raw, integ_fix = integ_raw ./ normfac)
end

fd = pieces(params)
ho_ = pieces(ho)

fig = Figure(size = (1180, 470))

axa = Axis(fig[1, 1]; xlabel = L"r~\mathrm{(fm)}", ylabel = L"u_{3D_1}(r)~\mathrm{(unit~norm)}",
    title = "the 3D₁ wavefunctions are identical")
xlims!(axa, 0, 3.0)
lines!(axa, fd.r_fm, fd.uD; color = (FD_COLOR, 0.95), linewidth = 2, label = "FD")
lines!(axa, ho_.r_fm, ho_.uD; color = (HO_COLOR, 0.95), linewidth = 2, linestyle = :dash, label = "HO")
axislegend(axa; position = :rt)

axb = Axis(fig[1, 2]; xlabel = L"r~\mathrm{(fm)}",
    ylabel = L"\mathrm{off\text{-}diagonal~tensor~integrand}",
    title = "off-diagonal element: bug vs fix")
xlims!(axb, 0, 2.0)
lines!(axb, fd.r_fm, fd.integ_raw; color = (FD_COLOR, 0.95), linewidth = 2, label = "FD")
lines!(axb, ho_.r_fm, ho_.integ_raw; color = (HO_COLOR, 0.95), linewidth = 2,
    linestyle = :dash, label = "HO (unnormalized, ×1/h)")
lines!(axb, ho_.r_fm, ho_.integ_fix; color = (FIX_COLOR, 0.95), linewidth = 2.5,
    linestyle = :dot, label = "HO (normalized fix)")
axislegend(axb; position = :rt)

mkpath(OUTDIR)
save(joinpath(OUTDIR, "wavefn_offdiagonal.png"), fig)
println("wrote wavefn_offdiagonal.png")
