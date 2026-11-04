# Appendix-D mock-meson transition matrix elements (page 40), promoted from the
# Table VI photon-decay audit into src with tests.
#
# A mock meson is the weakly-bound |q q̄⟩ state whose momentum-space radial wave
# Φ_L(p) is the spherical-Bessel transform of the reduced radial wave u(r), and
# whose mock mass is M̃ = ⟨E₁⟩ + ⟨E₂⟩. The two kernels are:
#
#   I_i(x,y)   = √(4 M̃ₓ M̃_y)/(M̃ₓ + M̃_y) · ∫ dp p² Φₓ(p) Φ_y(p) (1/mᵢ)(mᵢ/Eᵢ)^f
#   Eₙⁱ(x,y)   = |mᵢ / √(⟨Eᵢ⟩ₓ ⟨Eᵢ⟩_y)|^g · ∫ dr uₓ(r) u_y(r) rⁿ
#
# with Eᵢ = √(mᵢ² + p²) the relativistic quark energy and f, g the paper's fitted
# exponents (0.7 for I_i, 0.5 for Eₙⁱ). I_i drives the M1 magnetic moments; Eₙⁱ
# the E1/M2 multipole amplitudes and the hindered-transition recoil term.

"""
Mean relativistic quark energy `⟨E⟩ = ∫ p² Φ² √(m²+p²) dp` over a mock wave.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
mw = momentum_wave(wave, 0)
@assert QMT.mock_mean_energy(mw, 1.628) > 1.628
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
mock_mean_energy(mw::MomentumWave, m::Real) =
    momentum_expect(mw, p -> sqrt(float(m)^2 + p^2))

"""
Mock mass `M̃ = ⟨E₁⟩ + ⟨E₂⟩` of a mock wave with constituent masses `m1, m2`.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
mw = momentum_wave(wave, 0)
@assert QMT.mock_wave_mass(mw, 1.628, 1.628) > 2 * 1.628
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
mock_wave_mass(mw::MomentumWave, m1::Real, m2::Real) =
    mock_mean_energy(mw, m1) + mock_mean_energy(mw, m2)

"""
    mock_meson_overlap(mwx, mwy, m_emit; Mx, My, exponent=0.7) -> I_i (GeV⁻¹)

Appendix-D M1 transition overlap `I_i(x,y)` for emitting quark mass `m_emit`,
with `Mx`, `My` the mock masses of the two states (see `mock_wave_mass`):

    I_i = √(4 Mx My)/(Mx + My) · ∫ dp p² Φₓ Φ_y (1/m_emit)(m_emit/E)^exponent .

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
mw = momentum_wave(wave, 0)
mass = QMT.mock_wave_mass(mw, 1.628, 1.628)
@assert QMT.mock_meson_overlap(mw, mw, 1.628; Mx=mass, My=mass) > 0
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
function mock_meson_overlap(mwx::MomentumWave, mwy::MomentumWave, m_emit::Real;
                            Mx::Real, My::Real, exponent::Real = ELECTROMAGNETIC_DEFAULTS.magnetic_exponent)
    m = float(m_emit)
    pref = sqrt(4 * Mx * My) / (Mx + My)
    radial = momentum_overlap(
        mwx,
        mwy,
        p -> (1 / m) * (m / sqrt(m^2 + p^2))^exponent,
    )
    return pref * radial
end

"""
    mock_meson_radial_moment(wx, wy, Ex, Ey, m_emit; n=1, exponent=0) -> Eₙⁱ (GeV⁻ⁿ)

Appendix-D E1/M2 radial moment `Eₙⁱ(x,y)` on the position-space reduced waves,
with `Ex`, `Ey` the mean quark energies (see `mock_mean_energy`):

    Eₙⁱ = |m_emit / √(Ex Ey)|^exponent · ∫ dr uₓ(r) u_y(r) rⁿ .

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
mw = momentum_wave(wave, 0)
energy = QMT.mock_mean_energy(mw, 1.628)
@assert QMT.mock_meson_radial_moment(wave, wave, energy, energy, 1.628; n=2) > 0
```

## Related

- `e1_transition_amplitude` — assemble an E1 amplitude.
"""
function mock_meson_radial_moment(wx::RadialWave, wy::RadialWave,
                                  Ex::Real, Ey::Real, m_emit::Real;
                                  n::Integer = 1, exponent::Real = ELECTROMAGNETIC_DEFAULTS.electric_exponent)
    # This kernel IS `radial_overlap` with f(r) = r^n, times the Appendix-D
    # energy prefactor. Going through the interface normalizes both waves, which
    # this used to leave to the caller without saying so -- the same unstated
    # requirement that `_rel_momentum_average` carried. For the normalized waves
    # every solve now returns, the divisor is 1 and no number moves.
    radial = radial_overlap(wx, wy, x -> x^n)
    return abs(float(m_emit) / sqrt(Ex * Ey))^exponent * radial
end

# --- Radiative transition assembly (Eq. 22, Table VI) ------------------------
# The kernels above are the Appendix-D matrix elements; these three assemble
# them into the multipole amplitudes the paper tabulates. Promoted out of the
# Table VI audit script, where they sat next to the paper-comparison rows.

"""
Proton mass in GeV — the unit of the M1 moments, which Table VI prints as μ/μ_N.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.NUCLEON_MASS_GEV ≈ 0.93827
```

## Related

- `m1_transition_moment` — assemble a moment in nuclear magnetons.
"""
const NUCLEON_MASS_GEV = 0.93827

"""
    photon_momentum(M_parent_GeV, M_child_GeV) -> Float64

Photon momentum `q = (M² - M'²) / 2M` for the radiative transition
`parent → child + γ`, in GeV.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.photon_momentum(3.10, 2.98) > 0
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
photon_momentum(M_parent_GeV::Real, M_child_GeV::Real) =
    (float(M_parent_GeV)^2 - float(M_child_GeV)^2) / (2 * float(M_parent_GeV))

"""
    m1_transition_moment(singlet, triplet, m1_GeV, m2_GeV, terms) -> Float64

The M1 transition moment `μ/μ_N` for `V → P γ` between two mock mesons given by
their momentum-space waves, in Table VI's units of `e/2`.

`terms` is an iterable of `(coefficient, m_i)` pairs — one per contributing
quark line — summed as `Σ c·I_i`. The coefficients follow the transition-moment
rule `μ = e_q I_q - e_q̄ I_q̄` (the antiquark charge enters flipped), so
charmonium is `+4/3 I_c` and bottomonium `-2/3 I_b`.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
mw = momentum_wave(wave, 0)
moment = QMT.m1_transition_moment(mw, mw, 1.628, 1.628, [(4/3, 1.628)])
@assert QMT.m1_radiative_width(moment, 0.1) > 0
```

## Related

- `NUCLEON_MASS_GEV` — mass defining the magnetic-moment unit.
- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
- `m1_radiative_width` — primitive M1 width in GeV.
- `m1_recoil_moment` — M1 moment with the recoil correction.
- `mock_mean_energy` — mean constituent energy.
- `mock_meson_overlap` — Appendix-D magnetic overlap.
- `mock_wave_mass` — mock mass from a momentum wave.
- `neutral_m1_charge` — reference neutral-flavor charge factor.
- `observable_momentum_wave` — momentum representation for observables.
- `photon_recoil_form_factor` — optional Gaussian recoil factor.
"""
function m1_transition_moment(
    singlet::MomentumWave,
    triplet::MomentumWave,
    m1_GeV::Real,
    m2_GeV::Real,
    terms,
    ;
    exponent::Real = ELECTROMAGNETIC_DEFAULTS.magnetic_exponent,
)
    Mx = mock_wave_mass(singlet, m1_GeV, m2_GeV)
    My = mock_wave_mass(triplet, m1_GeV, m2_GeV)
    return sum(
        c * mock_meson_overlap(
            singlet, triplet, m_i; Mx = Mx, My = My, exponent = exponent,
        ) for (c, m_i) in terms
    ) * NUCLEON_MASS_GEV
end

"""
    e1_transition_amplitude(wave_S, mom_S, wave_P, mom_P, m_i, coeff_of_q,
                            M_parent_GeV, M_child_GeV; q = nothing) -> Float64

The E1 (and M2) amplitude `coeff(q) · E₁ⁱ · √(α q)` in `MeV^(1/2)`, with `q` in
MeV inside the square root — Table VI's amplitude convention.

`coeff_of_q` is the row's angular/charge factor as a function of `q` in GeV, so
the paper's printed q-dependence stays visible at the call site. Pass an
explicit `q` (GeV) to override the photon momentum implied by the two masses,
e.g. to use model rather than measured masses.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = OscillatorWave(0, 0.5, [1.0])
mw = momentum_wave(wave, 0)
pw = OscillatorWave(1, 0.5, [1.0])
amp = QMT.e1_transition_amplitude(wave, mw, pw, momentum_wave(pw, 1), 1.628, q -> (2/3)*q, 3.51, 3.10)
@assert isfinite(amp)
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
- `e1_angular_coefficient` — E1 angular factor.
- `mock_meson_radial_moment` — Appendix-D electric radial moment.
- `mock_mean_energy` — mean constituent energy.
"""
function e1_transition_amplitude(
    wave_S::RadialWave,
    mom_S::MomentumWave,
    wave_P::RadialWave,
    mom_P::MomentumWave,
    m_i::Real,
    coeff_of_q,
    M_parent_GeV::Real,
    M_child_GeV::Real;
    q = nothing,
    exponent::Real = ELECTROMAGNETIC_DEFAULTS.electric_exponent,
)
    q_GeV = isnothing(q) ? photon_momentum(M_parent_GeV, M_child_GeV) : float(q)
    E1 = mock_meson_radial_moment(
        wave_S, wave_P,
        mock_mean_energy(mom_S, m_i), mock_mean_energy(mom_P, m_i), m_i;
        n = 1, exponent = exponent,
    )
    return coeff_of_q(q_GeV) * E1 * sqrt(ALPHA_EM * 1000 * q_GeV)
end
