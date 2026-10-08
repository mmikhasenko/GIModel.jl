# Radiative observables built on the shared mock-meson overlap kernels.

"""
    m1_recoil_moment(singlet, triplet, m, coefficient, q)

Equal-flavor hindered M1 moment, in nuclear magnetons, including the
`-q² E₂/(24m)` term retained in Table VI footnote c. Waves are radial waves,
`m` and photon momentum `q` are in GeV. The charge coefficient multiplies
both the direct and recoil terms.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
excited = OscillatorWave(0, 0.5, [0.0, 1.0])
@assert isfinite(QMT.m1_recoil_moment(wave, excited, 1.628, 4/3, 0.2))
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
- `mock_meson_radial_moment` — Appendix-D electric radial moment.
"""
function m1_recoil_moment(
    singlet::RadialWave, triplet::RadialWave, m::Real,
    coefficient::Real, q::Real;
    magnetic_exponent::Real = ELECTROMAGNETIC_DEFAULTS.magnetic_exponent,
    electric_exponent::Real = ELECTROMAGNETIC_DEFAULTS.electric_exponent,
)
    m > 0 || throw(ArgumentError("constituent mass must be positive"))
    q >= 0 || throw(ArgumentError("photon momentum must be nonnegative"))
    sp, tp = momentum_wave(singlet, 0), momentum_wave(triplet, 0)
    direct = m1_transition_moment(
        sp, tp, m, m, [(coefficient, m)]; exponent = magnetic_exponent,
    )
    E2 = mock_meson_radial_moment(singlet, triplet,
        mock_mean_energy(sp, m), mock_mean_energy(tp, m), m;
        n = 2, exponent = electric_exponent)
    return direct - coefficient * q^2 / (24m) * E2 * NUCLEON_MASS_GEV
end

"""
    photon_recoil_form_factor(q; beta=0.40)

Amplitude form factor exp(-q²/(16β²)) for emission with recoil absorbed by a
light-quark system (Table VI footnote g). Both q and beta are in GeV.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.photon_recoil_form_factor(0.0) == 1.0
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
function photon_recoil_form_factor(q::Real; beta::Real = ELECTROMAGNETIC_DEFAULTS.recoil_beta_GeV)
    q >= 0 || throw(ArgumentError("photon momentum must be nonnegative"))
    beta > 0 || throw(ArgumentError("beta must be positive"))
    return exp(-q^2 / (16beta^2))
end

"""
    m1_radiative_width(moment, q; parent_spin=1)

M1 partial width in GeV from a moment in nuclear magnetons and photon momentum
in GeV. The spin average gives α μ² q³/(3 M_N²) for V -> P gamma and three
times that for P -> V gamma.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.m1_radiative_width(1.0, 0.1) > 0 # GeV
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
function m1_radiative_width(moment::Real, q::Real; parent_spin::Integer = 1)
    parent_spin in (0, 1) || throw(ArgumentError("parent_spin must be 0 or 1"))
    q >= 0 || throw(ArgumentError("photon momentum must be nonnegative"))
    return ALPHA_EM * abs2(moment) * q^3 /
           ((2parent_spin + 1) * NUCLEON_MASS_GEV^2)
end

"""
    neutral_m1_charge(flavor; isovector_left=false, isovector_right=false)

Pure-flavor M1 charge coefficient, before any physical-state mixing. For the
normalized nonstrange states, equal isospin gives e_u+e_d=1/3 and opposite
isospin gives e_u-e_d=1. A single strange/charm/bottom flavor gives twice its
quark charge. Perfect eta mixing factors must not be included here: the
physical-state composition supplies them.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.neutral_m1_charge(:c) == 4/3
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
function neutral_m1_charge(
    flavor::Symbol; isovector_left::Bool = false, isovector_right::Bool = false,
)
    if flavor == :q
        return xor(isovector_left, isovector_right) ? 1.0 : 1 / 3
    end
    (isovector_left || isovector_right) &&
        throw(ArgumentError("an isovector state must be nonstrange"))
    flavor == :c && return 4 / 3
    flavor in (:s, :b) && return -2 / 3
    throw(ArgumentError("unsupported neutral flavor $flavor"))
end

# Orbital matrix element <Lf mf| r_hat_mu |Li mi> of the unit vector.
function _unit_vector_element(Lf::Int, mf::Int, Li::Int, mi::Int, mu::Int)
    mf == mi + mu || return 0.0
    return sqrt((2Li + 1) / (2Lf + 1)) *
           Float64(CG(Li, 0, 1, 0, Lf, 0)) * Float64(CG(Li, mi, 1, mu, Lf, mf))
end

_photon_spin_element(Sf, mSf, Si, mSi, ::Nothing) = Float64(Sf == Si && mSf == mSi)
_photon_spin_element(Sf, mSf, Si, mSi, mu::Int) = _spin_factor(
    _QuarkEmission(), _coupled_spin(Sf, mSf), _coupled_spin(Si, mSi), mu)

"""
Photon angular rate between (Li,Si)Ji and (Lf,Sf)Jf: summed over final and
averaged over initial projections. Each operator is an (orbital, spin)
component pair; `nothing` is the spin identity.
"""
function _photon_angular_rate(Lf::Int, Sf::Int, Jf::Int, Li::Int, Si::Int, Ji::Int, operators)
    total = 0.0
    for Mi in -Ji:Ji, Mf in -Jf:Jf, (orbital, spin) in operators
        amplitude = 0.0
        for mLi in -Li:Li, mSi in -Si:Si, mLf in -Lf:Lf, mSf in -Sf:Sf
            mLi + mSi == Mi && mLf + mSf == Mf || continue
            amplitude += Float64(CG(Li, mLi, Si, mSi, Ji, Mi)) *
                         Float64(CG(Lf, mLf, Sf, mSf, Jf, Mf)) *
                         _unit_vector_element(Lf, mLf, Li, mLi, orbital) *
                         _photon_spin_element(Sf, mSf, Si, mSi, spin)
        end
        total += abs2(amplitude)
    end
    return total / (2Ji + 1)
end

# E1 is spin blind: the three dipole components with the spin identity.
const _E1_OPERATORS = Tuple((mu, nothing) for mu in -1:1)
# Spin-flip term sigma.(q x eps)(q.r) with q along z: r_hat_0 times sigma_{-lambda}.
const _SPIN_FLIP_OPERATORS = ((0, -1), (0, 1))

"""
    e1_angular_coefficient(J_P; singlet=false, parent_is_S=false)

Angular factor multiplying the neutral M1 charge coefficient times q E1. It is
sqrt(R/3), where R is the E1 angular rate from L⊗S coupling, so that
Γ = (4/3) α e² q³ R |E1|². E1 does not act on spin: P -> S gives 1/3 for both
triplet P_J -> S1 and singlet P1 -> S0, and S -> P carries an extra
sqrt((2J_P+1)/(2J_S+1)).

GI 1985 Table VI prints sqrt(2) q/3 for B -> pi gamma, 3sqrt(2) times the
q/9 of the equivalent A2 -> rho gamma row. That coefficient is not used here.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.e1_angular_coefficient(1; singlet=true) ≈ 1/3
```

## Related

- `e1_transition_amplitude` — assemble an E1 amplitude.
"""
function e1_angular_coefficient(J_P::Integer; singlet::Bool = false, parent_is_S::Bool = false)
    if singlet
        J_P == 1 || throw(ArgumentError("singlet P state must have J=1"))
    else
        J_P in 0:2 || throw(ArgumentError("triplet P state must have J in 0:2"))
    end
    S = singlet ? 0 : 1
    P, Sw = (1, S, Int(J_P)), (0, S, S)
    final, initial = parent_is_S ? (P, Sw) : (Sw, P)
    return sqrt(_photon_angular_rate(final..., initial..., _E1_OPERATORS) / 3)
end

# Denominator of the spin-flip P_J -> S0 amplitude, in units of m: sqrt(8/T),
# i.e. sqrt(120) for J=2 and sqrt(72) for J=1. The magnetization term
# (ê/2m) sigma.(q x eps*)(-i q.r_i) with r_i = r/2 and Gamma = 2 alpha q <|M|^2>
# gives 8; `MultipolePhotonEmission`, which reproduces the standard M1 and E1
# widths and the Karl-Meshkov-Rosner M2/E1 ratios, confirms it. GI 1985
# Table VI prints sqrt(60) m and 6m, sqrt(2) larger in amplitude.
function _spin_flip_denominator(J_P::Integer)
    rate = _photon_angular_rate(0, 0, 0, 1, 1, Int(J_P), _SPIN_FLIP_OPERATORS)
    rate > 0 || throw(ArgumentError("spin-flip P$(J_P) -> S0 photon emission vanishes"))
    return sqrt(8 / rate)
end

"""
    spin_flip_photon_amplitude(wave_S, wave_P, terms, J_P, q)

Spin-flip P -> S photon amplitude in MeV^(1/2): M2 for P2 -> S0,
E1 for P1 -> S0. The power of q alone does not determine the multipole. Each term is an emitting
quark (charge coefficient, constituent mass) pair; q is in GeV. The P2 and
P1 denominators are sqrt(120)m and sqrt(72)m, respectively; GI 1985 prints
sqrt(60)m and 6m, which overstate the rate by a factor of 2. Charge coefficients
include the relative sign of antiquark emission for this multipole and, for
unequal masses, the centre-of-mass weight 2 m_other/(m_1+m_2) of the emitter.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
pw = OscillatorWave(1, 0.5, [1.0])
@assert isfinite(QMT.spin_flip_photon_amplitude(wave, pw, [(4/3, 1.628)], 2, 0.2))
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
function spin_flip_photon_amplitude(
    wave_S::RadialWave, wave_P::RadialWave, terms, J_P::Integer, q::Real;
    exponent::Real = ELECTROMAGNETIC_DEFAULTS.electric_exponent,
)
    J_P in (1, 2) || throw(ArgumentError("spin-flip prescription requires J_P=1 or 2"))
    q >= 0 || throw(ArgumentError("photon momentum must be nonnegative"))
    factor = _spin_flip_denominator(J_P)
    mom_S, mom_P = momentum_wave(wave_S, 0), momentum_wave(wave_P, 1)
    return sum(terms) do (c, m)
        m > 0 || throw(ArgumentError("constituent mass must be positive"))
        e1_transition_amplitude(wave_S, mom_S, wave_P, mom_P, m,
            q -> c * q^2 / (factor * m), 1.0, 0.0;
            q = q, exponent = exponent)
    end
end
