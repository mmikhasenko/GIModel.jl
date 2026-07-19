# =============================================================================
# Table VII annihilation amplitudes — gluonic (QQ̄ → gluons) and leptonic
# (qq̄ → W → lν, qq̄ → γ* → l⁺l⁻) slices.
# =============================================================================
# Table VII (page 28) tabulates annihilation observables of the model
# wavefunctions. Both slices here run on the same Eq. (17)-style smeared
# momentum integral over Φ(p) (the jₗ transform of the radial wave) and have
# **zero free parameters** beyond the solved wavefunctions:
#   - gluonic widths ∝ |S_L(Ψ)|² (`_sL_smearing_factor`),
#   - leptonic decay constants via the mock-meson factors P_P, P'_A1, V_V, V'_V
#     of Table VII(a) / Appendix D (Eqs. D4-D6), which are the same kernel with
#     different [m/E]-type weights and mock/physical mass prefactors.

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

# -----------------------------------------------------------------------------
# Leptonic decay constants — the Table VII(a) mock-meson factors (Eqs. D4-D6).
# -----------------------------------------------------------------------------

# ∫u²dr = 1 copy of a radial wave (the factors carry absolute scale)
function _unit_norm_wave(radial::RadialWaveOnUniformMesh)
    nrm = sqrt(sum(abs2, radial.u) * radial.h)
    nrm > 0 || throw(ArgumentError("zero-norm radial wave"))
    return RadialWaveOnUniformMesh(radial.u ./ nrm, radial.r)
end

# K[w] = (2π)^(-3/2) ∫d³p (4π)^(-1/2) Φ_L(p) w(p): the Eq. (17) kernel with a
# configurable momentum weight (w ≡ m/E reproduces `_sL_smearing_factor`'s
# weight at L=0). `wave` must already be unit-normalized.
function _mock_momentum_kernel(wave::RadialWaveOnUniformMesh, L::Integer, w; npoints::Integer = 900)
    npts = max(npoints, 32)
    pmax = π / wave.h
    dp = pmax / (npts - 1)
    accum = 0.0
    for k in 1:npts
        p = (k - 1) * dp
        weight = (k == 1 || k == npts) ? 0.5 : 1.0
        accum += weight * p^2 * _momentum_radial_wave(wave, p, L) * w(p)
    end
    return sqrt(2 / π) * accum * dp / sqrt(4π)
end

"""
    mock_meson_mass(radial, m1_GeV, m2_GeV; L=0, npoints=900)

The mock-meson mass `M̃ = <E₁> + <E₂>` (GeV): the free quark-pair energy
averaged over the momentum-space wavefunction `|Φ_L(p)|²` of the radial wave.
This is the `M̃` appearing in the Table VII(a) leptonic-factor prefactors.
"""
function mock_meson_mass(
    radial::RadialWaveOnUniformMesh,
    m1_GeV::Real,
    m2_GeV::Real;
    L::Integer = 0,
    npoints::Integer = 900,
)
    wave = _unit_norm_wave(radial)
    npts = max(npoints, 32)
    pmax = π / wave.h
    dp = pmax / (npts - 1)
    accE, accN = 0.0, 0.0
    for k in 1:npts
        p = (k - 1) * dp
        weight = (k == 1 || k == npts) ? 0.5 : 1.0
        Φ2 = _momentum_radial_wave(wave, p, L)^2
        accE += weight * p^2 * Φ2 * (sqrt(m1_GeV^2 + p^2) + sqrt(m2_GeV^2 + p^2))
        accN += weight * p^2 * Φ2
    end
    accN > 0 || throw(ArgumentError("mock_meson_mass: vanishing momentum norm"))
    return accE / accN
end

# kind => (orbital L of the wavefunction, needs equal masses)
const LEPTONIC_FACTOR_KINDS = Dict(
    :P_P => (0, false),    # ¹S₀ pseudoscalar, weight √(m₁m₂/E₁E₂)
    :V_V => (0, false),    # ³S₁ vector,       weight √(m₁m₂/E₁E₂)
    :Vp_V => (2, true),    # ³D₁ vector,       weight (m/E)(1 − m/E)
    :Pp_A1 => (1, true),   # ³P₁ axial,        weight m·p/E²
)

"""
    leptonic_decay_factor(kind, radial, m1_GeV, m2_GeV, M_GeV; npoints=900)

The dimensionless Table VII(a) mock-meson leptonic factor (Eqs. D4-D6):

    P_P   = M⁻¹ M̃^(-1/2) K[√(m₁m₂/E₁E₂)]     (¹S₀ wave, L=0)
    V_V   = M⁻² M̃^(+1/2) K[√(m₁m₂/E₁E₂)]     (³S₁ wave, L=0)
    V'_V  = M⁻² M̃^(+1/2) K[(m/E)(1 − m/E)]   (³D₁ wave, L=2, equal masses)
    P'_A1 = M⁻² M̃^(+1/2) K[m·p/E²]           (³P₁ wave, L=1, equal masses)

with `K[w] = (2π)^(-3/2) ∫d³p (4π)^(-1/2) Φ_L(p) w(p)`, `M` the meson mass and
`M̃` the mock mass [`mock_meson_mass`](@ref) of the same wave. The tabulated
amplitude is a quark-charge/flavor coefficient times this factor (e.g.
`f_π/M_π = 2√3 P_π`, `f_ψ = (16/3)^(1/2) V_ψ`); the paper's implicit `m/E`
exponent is unity, as stated below the formula block. `radial` must be the
wave of the orbital the kind expects; its overall sign propagates to the
factor (fix a phase convention upstream for sign comparisons).
"""
function leptonic_decay_factor(
    kind::Symbol,
    radial::RadialWaveOnUniformMesh,
    m1_GeV::Real,
    m2_GeV::Real,
    M_GeV::Real;
    npoints::Integer = 900,
)
    haskey(LEPTONIC_FACTOR_KINDS, kind) ||
        throw(ArgumentError("unknown leptonic factor kind $kind (expected one of $(keys(LEPTONIC_FACTOR_KINDS)))"))
    L, equal_only = LEPTONIC_FACTOR_KINDS[kind]
    m1, m2, M = float(m1_GeV), float(m2_GeV), float(M_GeV)
    (equal_only && !isapprox(m1, m2)) &&
        throw(ArgumentError("$kind is defined for equal constituent masses (got $m1, $m2)"))
    M > 0 || throw(ArgumentError("meson mass must be positive"))

    wave = _unit_norm_wave(radial)
    Mtilde = mock_meson_mass(wave, m1, m2; L = L, npoints = npoints)
    w = if kind === :P_P || kind === :V_V
        p -> sqrt(m1 * m2 / (sqrt(m1^2 + p^2) * sqrt(m2^2 + p^2)))
    elseif kind === :Vp_V
        p -> (E = sqrt(m1^2 + p^2); (m1 / E) * (1 - m1 / E))
    else # :Pp_A1
        p -> m1 * p / (m1^2 + p^2)
    end
    K = _mock_momentum_kernel(wave, L, w; npoints = npoints)
    prefactor = kind === :P_P ? 1 / (M * sqrt(Mtilde)) : sqrt(Mtilde) / M^2
    return prefactor * K
end

# -----------------------------------------------------------------------------
# Two-photon amplitudes — the Table VII(b) γγ formulas (page 28).
# -----------------------------------------------------------------------------

"""Electromagnetic fine-structure constant used in the Table VII(b) γγ formulas."""
const ALPHA_EM = 1 / 137.036

# raw momentum moment ∫ p² Φ_L(p) w(p) dp over the (unit-normalized) wave
function _momentum_moment(wave::RadialWaveOnUniformMesh, L::Integer, w; npoints::Integer = 900)
    npts = max(npoints, 32)
    pmax = π / wave.h
    dp = pmax / (npts - 1)
    accum = 0.0
    for k in 1:npts
        p = (k - 1) * dp
        weight = (k == 1 || k == npts) ? 0.5 : 1.0
        accum += weight * p^2 * _momentum_radial_wave(wave, p, L) * w(p)
    end
    return accum * dp
end

"""
    two_photon_amplitude(kind, radial, m_GeV, M_GeV, q_eff; npoints=900)

Two-photon annihilation amplitude (GeV^(1/2), amplitude² = Γ) of Table VII(b):

    A(P→γγ)   = √6   q_eff (α/m) (M/M̃)^(3/2) (1/2π)     ∫d³p φ_P(p) [m/E]
    A(³P₂→γγ) = −√(4/5) q_eff (α/m) (M/M̃)^(3/2) (2/π)^(1/2) ∫dp p² Φ(p) [m·p/E²]

for `kind` `:P` (S-wave pseudoscalar, `radial` the ¹S₀ wave) or `:P2` (the
³P₂ tensor, `radial` the P-wave). `m_GeV` is the constituent quark mass (equal
masses), `M_GeV` the meson mass, `M̃` the mock mass [`mock_meson_mass`](@ref),
and `q_eff = Σ aᵢ eᵢ²` the state's effective squared charge (the flavor
amplitude weighted sum of quark charges; e.g. `(e_u²−e_d²)/√2` for a π⁰-like
isovector, `4/9` for cc̄). For the `:P` S-wave, `∫d³p φ_P(m/E) = √(4π) ∫p²Φ(m/E)dp`.
The sign follows `q_eff` and the wave's phase convention.
"""
function two_photon_amplitude(
    kind::Symbol,
    radial::RadialWaveOnUniformMesh,
    m_GeV::Real,
    M_GeV::Real,
    q_eff::Real;
    npoints::Integer = 900,
)
    m, M = float(m_GeV), float(M_GeV)
    M > 0 || throw(ArgumentError("meson mass must be positive"))
    wave = _unit_norm_wave(radial)
    if kind === :P
        Mtilde = mock_meson_mass(wave, m, m; L = 0, npoints = npoints)
        moment = _momentum_moment(wave, 0, p -> m / sqrt(m^2 + p^2); npoints = npoints)
        integral = sqrt(4π) * moment                    # ∫d³p φ_P(m/E), S-wave
        return sqrt(6) * q_eff * (ALPHA_EM / m) * (M / Mtilde)^1.5 * (1 / (2π)) * integral
    elseif kind === :P2
        # ³P₂ is a P-wave (orbital L=1); the extra p in the [m·p/E²] weight is the
        # P-wave momentum factor.
        Mtilde = mock_meson_mass(wave, m, m; L = 1, npoints = npoints)
        moment = _momentum_moment(wave, 1, p -> m * p / (m^2 + p^2); npoints = npoints)
        return -sqrt(4 / 5) * q_eff * (ALPHA_EM / m) * (M / Mtilde)^1.5 * sqrt(2 / π) * moment
    end
    throw(ArgumentError("unknown two-photon kind $kind (expected :P or :P2)"))
end

# -----------------------------------------------------------------------------
# Charge radii — the Table VII(d) formula (page 28).
# -----------------------------------------------------------------------------

# <(m/E)^power>_φ = ∫ p² |Φ_0(p)|² (m/E)^power dp over the unit-normalized wave
function _rel_momentum_average(wave::RadialWaveOnUniformMesh, m::Real, power::Real; npoints::Integer = 900)
    npts = max(npoints, 32)
    pmax = π / wave.h
    dp = pmax / (npts - 1)
    accum = 0.0
    for k in 1:npts
        p = (k - 1) * dp
        weight = (k == 1 || k == npts) ? 0.5 : 1.0
        Φ = _momentum_radial_wave(wave, p, 0)
        accum += weight * p^2 * Φ^2 * (m / sqrt(m^2 + p^2))^power
    end
    return accum * dp
end

"""Conversion (ħc)² : an r² in GeV⁻² is `HBARC_FM2 · r²` in fm²."""
const HBARC_FM2 = 0.19733^2

"""
    charge_radius_squared(radial, m1_GeV, e1, m2_GeV, e2; f=0.2, npoints=900)

Mean-square charge radius `r_E²` (GeV⁻²) of a qq̄ meson, Table VII(d):

    r_E² = Σᵢ eᵢ [ <rᵢ²> + (3/4mᵢ²) ∫d³p |φ(p)|² (mᵢ/Eᵢ)^{2f} ]

`rᵢ = (mⱼ/M)·r` is quark i's position from the meson center of mass, so
`<rᵢ²> = (mⱼ/M)² <r²>` with `<r²>` the position-space expectation over the
`radial` (¹S₀) wave; the second term is the relativistic smearing of the
quark-position operator, `f` the fit exponent (0.2, fitted to the π⁺). `e1`,
`e2` are the quark charges (e.g. `+2/3`, `+1/3` for the u and d̄ of π⁺;
`−1/3`, `+1/3` for the d and s̄ of K⁰). Multiply by [`HBARC_FM2`](@ref) for fm².
"""
function charge_radius_squared(
    radial::RadialWaveOnUniformMesh,
    m1_GeV::Real,
    e1::Real,
    m2_GeV::Real,
    e2::Real;
    f::Real = 0.2,
    npoints::Integer = 900,
)
    m1, m2 = float(m1_GeV), float(m2_GeV)
    Mtot = m1 + m2
    Mtot > 0 || throw(ArgumentError("total constituent mass must be positive"))
    u, h = radial.u, radial.h
    nrm = sum(abs2, u) * h
    nrm > 0 || throw(ArgumentError("charge_radius_squared: zero-norm radial wave"))
    r2 = sum(@. u^2 * radial.r^2) * h / nrm
    wave = RadialWaveOnUniformMesh(u ./ sqrt(nrm), radial.r)
    term(mi, ei, mj) =
        ei * ((mj / Mtot)^2 * r2 + (3 / (4 * mi^2)) * _rel_momentum_average(wave, mi, 2f; npoints = npoints))
    return term(m1, e1, m2) + term(m2, e2, m1)
end
