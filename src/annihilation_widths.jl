# =============================================================================
# Table VII annihilation widths — gluonic (QQ̄ → gluons) slice.
# =============================================================================
# Lowest-order QCD annihilation of heavy quarkonia into gluons (Table VII part
# (c), page 28). Every width is proportional to |S_L(Ψ)|², the Eq. (17)
# smeared wavefunction-at-origin already implemented for the pseudoscalar
# annihilation block (`_sL_smearing_factor`). This slice therefore has **zero
# free parameters**: the model wavefunction, the constituent quark mass, and the
# running coupling α_s(M) fix the amplitudes outright.

"""
    wavefunction_origin_smearing(radial, mass_GeV; L=0, npoints=900)

The Eq. (17) smeared wavefunction-at-origin

    S_L(Ψ) = (2π)^(-3/2) ∫ d³p (4π)^(-1/2) Φ(p) [p/E]^L (m/E)

for a reduced radial wave `radial` (u(r) = r·R(r)) of a constituent quark of
mass `mass_GeV`, with `Φ(p)` the momentum-space wave from the jₗ transform and
`E = √(m² + p²)`. The wave is normalized to `∫u²dr = 1` internally, so the
result carries the absolute GeV^(3/2) scale of a wavefunction at the origin; its
sign follows the radial phase of `radial` (which alternates with the number of
radial nodes, matching the Table VII sign pattern within a quarkonium family).
"""
function wavefunction_origin_smearing(
    radial::RadialWaveOnUniformMesh,
    mass_GeV::Real;
    L::Integer = 0,
    npoints::Integer = 900,
)
    u = radial.u
    nrm = sqrt(sum(abs2, u) * radial.h)
    nrm > 0 || throw(ArgumentError("wavefunction_origin_smearing: zero-norm radial wave"))
    wave = RadialWaveOnUniformMesh(u ./ nrm, radial.r)
    input = PseudoscalarAnnihilationBasisInput("", float(mass_GeV), 0.0, wave)
    return _sL_smearing_factor(FDMomentumIntegralSmearing(npoints), input, L)
end

# The four lowest-order gluonic channels of Table VII(c). Each key maps a decay
# channel to the QCD width prefactor P so that Γ = P · |S_L|²:
#   :S0_2g  Γ(¹S₀→2g) = 8π α_s²/(3 m_Q²) |S₀|²
#   :S1_3g  Γ(³S₁→3g) = 40(π²−9)/(81 m_Q²) α_s³ |S₀|²
#   :P2_2g  Γ(³P₂→2g) = 32π α_s²/(45 m_Q²) |S₁|²
#   :P0_2g  Γ(³P₀→2g) = 8π α_s²/(3 m_Q²) |S₁|²
# (α_s = α_s(μ) with μ the mass of the decaying meson; S_L is S₀ for the S-wave
# channels and S₁ for the P-wave channels.)
const GLUONIC_CHANNELS = (:S0_2g, :S1_3g, :P2_2g, :P0_2g)

function gluonic_width_prefactor(channel::Symbol, alpha_s::Real, mQ::Real)
    a = float(alpha_s)
    m2 = float(mQ)^2
    channel === :S0_2g && return 8π * a^2 / (3 * m2)
    channel === :S1_3g && return 40 * (π^2 - 9) / (81 * m2) * a^3
    channel === :P2_2g && return 32π * a^2 / (45 * m2)
    channel === :P0_2g && return 8π * a^2 / (3 * m2)
    throw(ArgumentError("unknown gluonic channel $channel (expected one of $(GLUONIC_CHANNELS))"))
end

"""
    gluonic_annihilation_amplitude(channel, S_L, alpha_s, mQ)

Signed lowest-order gluonic annihilation amplitude in GeV^(1/2); the amplitude
squared is the width Γ (this matches the Table VII "predicted amplitude" column,
which tabulates √Γ with the sign of `S_L`). `channel` ∈ `GLUONIC_CHANNELS`;
`S_L` is `wavefunction_origin_smearing` at `L=0` for the S-wave channels
(`:S0_2g`, `:S1_3g`) and `L=1` for the P-wave channels (`:P2_2g`, `:P0_2g`);
`alpha_s = α_s(M)` at the decaying-meson mass; `mQ` the constituent quark mass.
"""
gluonic_annihilation_amplitude(channel::Symbol, S_L::Real, alpha_s::Real, mQ::Real) =
    sqrt(gluonic_width_prefactor(channel, alpha_s, mQ)) * float(S_L)

"""
    gluonic_annihilation_width(channel, S_L, alpha_s, mQ)

Lowest-order gluonic annihilation width Γ (GeV) = amplitude². See
[`gluonic_annihilation_amplitude`](@ref).
"""
gluonic_annihilation_width(channel::Symbol, S_L::Real, alpha_s::Real, mQ::Real) =
    gluonic_width_prefactor(channel, alpha_s, mQ) * float(S_L)^2
