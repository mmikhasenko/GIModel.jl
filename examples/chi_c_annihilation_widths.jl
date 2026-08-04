#!/usr/bin/env julia
# =============================================================================
# chi_c0 vs chi_c2 annihilation: two-gluon and two-photon widths
# =============================================================================
# Worked example on top of the Table VII machinery (`src/annihilation_widths.jl`).
# It answers: why is Gamma(chi_c0 -> gg) about six times Gamma(chi_c2 -> gg),
# when the lowest-order spin algebra alone says 15/4 = 3.75?
#
# Wave treatment matches `GIPaper/scripts/audit_table_vii.jl`: the paper-order
# finite-HO-basis full diagonalization of H_central + the calibrated spin-orbit
# + tensor operator, so each ^3P_J gets its own J-distorted radial wave.
#
# Three blocks:
#   (1) gg widths      — `gluonic_annihilation_amplitude`, channels :P0_2g/:P2_2g
#   (2) gamma-gamma    — Table VII(b) ^3P_2 kernel, plus its ^3P_0 partner
#   (3) decomposition  — the same P-wave solved with V_spin = 0, which isolates
#                        how much of the 0/2 ratio is spin algebra and how much
#                        is J-dependent distortion of the wave at the origin.
#
# Findings (see the summary the script prints):
#   Gamma(gg):     chi_c0 5.60 MeV, chi_c2 0.910 MeV, ratio 6.15
#   Gamma(gamgam): chi_c0 3.25 keV, chi_c2 0.623 keV, ratio 5.21
#   The gg ratio factorizes as 3.75 (spin) x 1.61 (|S_1|^2 distortion), the last
#   2% from alpha_s(M) at two different meson masses. Switching V_spin off
#   collapses the wave factor to 1 and the ratio to the bare 3.75.
#
# Run:  julia examples/chi_c_annihilation_widths.jl
# =============================================================================

using Pkg
Pkg.activate(@__DIR__; io = devnull)     # examples/Project.toml — GIModel + CairoMakie + PlutoUI
Pkg.instantiate(; io = devnull)

using Printf
using GIModel

const PARAMS_PATH = joinpath(dirname(@__DIR__), "data", "parameters.provisional.toml")
const NGRID, RMAX, NPTS, NB = 1200, 24.0, 900, 24

# The one internal reach in this file: `fine_structure_grid_operator` is exported
# and takes a mesh, but the mesh constructor itself is not. Reimplementing the
# convention here would risk drifting from the mesh `ho_full_distorted_states`
# builds from the same (ngrid, rmax) — the operator and the solve must agree.

# Phase convention of the Table VII audit: outermost antinode positive, so the
# wavefunction-at-origin flips sign once per radial node.
function fix_outer_antinode_positive!(u)
    peak = maximum(abs, u)
    i = findlast(x -> abs(x) > 0.2 * peak, u)
    (i !== nothing && u[i] < 0) && (u .*= -1)
    return u
end

params, mq = load_parameters_and_quark_masses(PARAMS_PATH)
params_ho = with_basis(params, HarmonicOscillatorBasis)
const MC = mq["c"]
const R, H = GIModel.radial_grid(NGRID, RMAX)

"""Ground-state cc̄ P-wave under the spin operator `V` (`V = 0` → central wave)."""
function pwave_level(V)
    vals, vecs, r = ho_full_distorted_states(params_ho, ConstituentMasses(MC, MC), 1, V;
        nlevels = 4, ngrid = NGRID, rmax = RMAX, nbasis = NB)
    u = fix_outer_antinode_positive!(copy(vecs[:, 1]))
    wave = RadialWaveOnUniformMesh(u, r)
    nrm = sum(abs2, u) * wave.h
    rrms = sqrt(sum(@. u^2 * wave.r^2) * wave.h / nrm)
    S1 = wavefunction_origin_smearing(wave, MC; L = 1, npoints = NPTS)
    return (M = vals[1], wave = wave, S1 = S1, rrms = rrms)
end

spin_operator(J) = fine_structure_grid_operator(params, ConstituentMasses(MC, MC), J, R, H; L = 1)

chi = Dict(J => pwave_level(spin_operator(J)) for J in (0, 1, 2))
central = pwave_level(zeros(NGRID, NGRID))

# --- (1) two gluons ----------------------------------------------------------
# Gamma(^3P_0 -> 2g) = 8pi as^2/(3 mQ^2) |S_1|^2
# Gamma(^3P_2 -> 2g) = 32pi as^2/(45 mQ^2)|S_1|^2      (ratio 15/4 at fixed S_1)
# ^3P_1 has no entry: Landau-Yang forbids J=1 -> two massless vectors.
function gluonic(st, channel)
    as = GIModel.alpha_s_q(st.M)
    amp = gluonic_annihilation_amplitude(channel, st.S1, as, MC)      # GeV^(1/2)
    return (as = as, amp = amp, width = amp^2)
end

gg0 = gluonic(chi[0], :P0_2g)
gg2 = gluonic(chi[2], :P2_2g)

# --- (2) two photons ---------------------------------------------------------
# Table VII(b) tabulates only the light ^3P_2 rows (A2/f/f'), so
# `two_photon_amplitude` ships :P and :P2 but no :P0. The ^3P_0 partner is the
# SAME kernel — same wave, same mock mass, same momentum weight [m p/E^2] — with
# sqrt(3) in place of sqrt(4/5) out front. So the ^3P_0 amplitude is obtained by
# running the :P2 path on the ^3P_0 wave and rescaling by the coefficient ratio
#   sqrt(3)/sqrt(4/5) = sqrt(15)/2,
# which is what makes Gamma(^3P_0)/Gamma(^3P_2) = 15/4 at fixed wave — the same
# 15/4 already carried by the gluonic prefactors above. (Overall sign is a
# rephasing of the state and does not enter the width.)
const QCC = 4 / 9                              # effective squared charge of cc̄
const P0_OVER_P2_COEFF = sqrt(3) / sqrt(4 / 5)

function two_photon(st; p0::Bool)
    A = two_photon_amplitude(:P2, st.wave, MC, st.M, QCC; npoints = NPTS)   # GeV^(1/2)
    p0 && (A *= P0_OVER_P2_COEFF)
    Mt = mock_meson_mass(st.wave, MC, MC; L = 1, npoints = NPTS)
    return (Mt = Mt, A = A, width = A^2)
end

yy0 = two_photon(chi[0]; p0 = true)
yy2 = two_photon(chi[2]; p0 = false)

# --- report ------------------------------------------------------------------
to_mev_amp(a) = a * sqrt(1000)              # GeV^(1/2) -> MeV^(1/2)
to_kev_amp(a) = a / sqrt(1e-6)              # GeV^(1/2) -> keV^(1/2)

@printf("m_c = %.4f GeV   (constituent, from %s)\n\n", MC, basename(PARAMS_PATH))

println("P-wave cc̄ solutions (ground state of each channel)")
@printf("  %-10s %8s %8s %14s\n", "wave", "M(GeV)", "S_1", "r_rms(GeV^-1)")
@printf("  %-10s %8.3f %8.4f %14.3f\n", "central", central.M, central.S1, central.rrms)
for J in (0, 1, 2)
    @printf("  %-10s %8.3f %8.4f %14.3f\n", "3P_$J", chi[J].M, chi[J].S1, chi[J].rrms)
end
println()

println("Two-gluon widths   [Gamma = prefactor(channel, alpha_s(M), m_c) * |S_1|^2]")
@printf("  %-8s alpha_s = %.3f   amp = %+7.3f MeV^1/2   Gamma = %8.3f MeV\n",
    "chi_c0", gg0.as, to_mev_amp(gg0.amp), 1e3 * gg0.width)
@printf("  %-8s alpha_s = %.3f   amp = %+7.3f MeV^1/2   Gamma = %8.3f MeV\n",
    "chi_c2", gg2.as, to_mev_amp(gg2.amp), 1e3 * gg2.width)
@printf("  paper Table VII(c): 2.5 and 0.88 MeV^1/2  (model/paper = %.2f, %.2f)\n\n",
    abs(to_mev_amp(gg0.amp)) / 2.5, abs(to_mev_amp(gg2.amp)) / 0.88)

println("Two-photon widths  [Table VII(b) kernel; ^3P_0 coefficient sqrt(3)]")
@printf("  %-8s M~ = %.3f GeV   amp = %+7.3f keV^1/2   Gamma = %8.3f keV\n",
    "chi_c0", yy0.Mt, to_kev_amp(yy0.A), 1e6 * yy0.width)
@printf("  %-8s M~ = %.3f GeV   amp = %+7.3f keV^1/2   Gamma = %8.3f keV\n\n",
    "chi_c2", yy2.Mt, to_kev_amp(yy2.A), 1e6 * yy2.width)

wave_factor = (chi[0].S1 / chi[2].S1)^2
println("Why the gg ratio is not 15/4")
@printf("  spin algebra (two transverse gluons)        %.3f\n", 15 / 4)
@printf("  |S_1(3P_0)|^2 / |S_1(3P_2)|^2 distortion    %.3f\n", wave_factor)
@printf("  product                                     %.3f\n", (15 / 4) * wave_factor)
@printf("  actual Gamma(gg) ratio                      %.3f   (residual: alpha_s(M))\n",
    gg0.width / gg2.width)
@printf("  same ratio with V_spin = 0 for both waves   %.3f\n",
    gluonic(central, :P0_2g).width / gluonic(central, :P2_2g).width)
@printf("  gamma-gamma ratio                           %.3f   (extra (M/M~)^(3/2) tilt)\n\n",
    yy0.width / yy2.width)

@printf("Gamma(gg)/Gamma(gamgam): chi_c0 %.0f, chi_c2 %.0f  [LO (2/9)(alpha_s/alpha)^2/e_c^4 = %.0f]\n",
    gg0.width / yy0.width, gg2.width / yy2.width,
    (2 / 9) * (gg0.as / ALPHA_EM)^2 / (2 / 3)^4)
println("The two slices use different relativistic weights ([m/E][p/E] vs [m p/E^2])")
println("and the gamma-gamma kernel carries (M/M~)^(3/2), so they do not cancel the")
println("way the nonrelativistic |R'(0)|^2 formulas do.")
