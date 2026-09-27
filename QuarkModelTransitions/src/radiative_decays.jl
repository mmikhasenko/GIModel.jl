# Radiative observables built on the shared mock-meson overlap kernels.

"""
    m1_recoil_moment(singlet, triplet, m, coefficient, q)

Equal-flavor hindered M1 moment, in nuclear magnetons, including the
`-q² E₂/(24m)` term retained in Table VI footnote c. Waves are radial waves,
`m` and photon momentum `q` are in GeV. The charge coefficient multiplies
both the direct and recoil terms.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
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
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
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
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
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
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
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

"""
    e1_angular_coefficient(J_P; singlet=false, parent_is_S=false)

Angular factor multiplying the neutral M1 charge coefficient times q E1.
Triplet P -> S transitions carry 1/3, with sqrt((2J_P+1)/3) for the inverse
S -> P direction. The singlet P1 -> S0 factor is sqrt(2).

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
@assert QMT.e1_angular_coefficient(1; singlet=true) == sqrt(2.0)
```

## Related

- `e1_transition_amplitude` — assemble an E1 amplitude.
"""
function e1_angular_coefficient(J_P::Integer; singlet::Bool = false, parent_is_S::Bool = false)
    if singlet
        J_P == 1 || throw(ArgumentError("singlet P state must have J=1"))
        parent_is_S && throw(ArgumentError("inverse singlet E1 is not provided"))
        return sqrt(2.0)
    end
    J_P in 0:2 || throw(ArgumentError("triplet P state must have J in 0:2"))
    return (parent_is_S ? sqrt((2J_P + 1) / 3) : 1.0) / 3
end

"""
    spin_flip_photon_amplitude(wave_S, wave_P, terms, J_P, q)

Spin-flip P -> S photon amplitude in MeV^(1/2): M2 for P2 -> S0,
E1 for P1 -> S0. The power of q alone does not determine the multipole. Each term is an emitting
quark (charge coefficient, constituent mass) pair; q is in GeV. The P2 and
P1 denominators are sqrt(60)m and 6m, respectively. Charge coefficients
include the relative sign of antiquark emission for this multipole.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
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
    factor = J_P == 2 ? sqrt(60.0) : 6.0
    mom_S, mom_P = momentum_wave(wave_S, 0), momentum_wave(wave_P, 1)
    return sum(terms) do (c, m)
        m > 0 || throw(ArgumentError("constituent mass must be positive"))
        e1_transition_amplitude(wave_S, mom_S, wave_P, mom_P, m,
            q -> c * q^2 / (factor * m), 1.0, 0.0;
            q = q, exponent = exponent)
    end
end
