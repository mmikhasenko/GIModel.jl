# Public API (exported from GIModel.jl): StrongDecayModel, decay_momentum,
# reduced_decay_amplitude, strong_decay_amplitude, calibrate_strong_decay_model
#
# Sec. IV / Table IV-V strong-decay amplitudes in the paper's two-parameter
# harmonic-oscillator convention (see docs/observable_ledger.md). A Table V row
# with coefficient `c`, reduced-amplitude class `X`, and orbital power `L` is
#
#   amp = c * X(qbar) * qbar^L * sqrt(q / 2pi) * exp(-q^2 / 16 beta^2)
#
# in MeV^(1/2) with `qbar = q/beta`. The structure-independent classes
# (A, A', A'', A0) all take the single fitted value `A`; the structure-dependent
# classes share `S0 = 3 h beta` with class-specific `qbar^2` subtractions.

"""
    StrongDecayModel(A, S0, beta_GeV)

Two-parameter Table IV/V decay model: `A` is the structure-independent reduced
amplitude (fit to `rho -> pi pi`), `S0 = 3 h beta` the structure-dependent
strength (fit to `B -> [omega pi]_S`), with oscillator scale `beta` (0.40 GeV
in the paper).
"""
struct StrongDecayModel
    A::Float64
    S0::Float64
    beta_GeV::Float64
end

"""
    decay_momentum(M, m1, m2)

Two-body breakup momentum of `M -> m1 + m2` (GeV); zero below threshold.
"""
function decay_momentum(M::Real, m1::Real, m2::Real)
    M <= m1 + m2 && return 0.0
    return sqrt((M^2 - (m1 + m2)^2) * (M^2 - (m1 - m2)^2)) / (2M)
end

"""
    reduced_decay_amplitude(model, class, qbar)

Table IV reduced partial-wave amplitude for `class` at `qbar = q/beta`.
Structure-independent classes (`:A`, `:Aprime`, `:Adoubleprime`, `:A0`) return
the fitted `A` (the paper's `A' ≃ A'' ≃ A0 ≃ A` convention); structure-dependent
classes return `S0 - k A qbar^2` with `k = 1/2, 3/10, 3/4` for `:S`, `:D`, `:P`.

The charm analogues (`:A_c`, `:S_c`, footnote d of Table V) reuse the identical
reduced-amplitude algebra with the SAME fitted strengths `A`, `S0` and
`beta_c = beta`; only the emission form factor and (for the A_c P-wave rows)
an explicit recoil multiplier differ, and those live in
`strong_decay_amplitude` / `_charm_form_factor`, not here.
"""
function reduced_decay_amplitude(model::StrongDecayModel, class::Symbol, qbar::Real)
    class in (:A, :Aprime, :Adoubleprime, :A0) && return model.A
    class === :A_c && return model.A            # charm analogue of :A (footnote d)
    class === :S && return model.S0 - 0.5 * model.A * qbar^2
    class === :S_c && return model.S0 - 0.5 * model.A * qbar^2  # charm analogue of :S
    class === :D && return model.S0 - 0.3 * model.A * qbar^2
    class === :P && return model.S0 - 0.75 * model.A * qbar^2
    throw(ArgumentError("unknown reduced-amplitude class `$class`"))
end

# Default charm-meson quark masses (GeV): m_c and the light spectator m_d, from
# data/parameters.provisional.toml (Table II). Used by the charm form factor.
const CHARM_M_C_GEV = 1.628
const CHARM_M_D_GEV = 0.220

# The Table V caption suppresses (q/2pi)^(1/2) exp(-q^2/16 beta^2) from the
# formula column; q enters the square root in MeV so amplitudes come out in
# MeV^(1/2), matching the numeric column.
function _suppressed_factor(q_GeV::Real, beta_GeV::Real)
    return sqrt(1000q_GeV / (2π)) * exp(-q_GeV^2 / (16 * beta_GeV^2))
end

# Charmed-meson suppressed factor (Table V footnote d): the Gaussian form
# factor is replaced by exp[-(1/4)(m_c/(m_c+m_d))^2 q^2/beta_c^2] with
# beta_c = beta numerically. The sqrt(q/2pi) prefactor is unchanged.
function _charm_suppressed_factor(
    q_GeV::Real, beta_c_GeV::Real, m_c_GeV::Real, m_d_GeV::Real,
)
    r = m_c_GeV / (m_c_GeV + m_d_GeV)
    return sqrt(1000q_GeV / (2π)) * exp(-0.25 * r^2 * q_GeV^2 / beta_c_GeV^2)
end

"""
    charm_decay_amplitude(model, coefficient, class, qbar_power, q_GeV;
                          recoil=false, m_c_GeV=CHARM_M_C_GEV, m_d_GeV=CHARM_M_D_GEV,
                          beta_c_GeV=model.beta_GeV)

Table V charmed-meson amplitude in `MeV^(1/2)` (footnote d). Identical in form
to `strong_decay_amplitude` but with the charmed Gaussian form factor
`exp[-(1/4)(m_c/(m_c+m_d))^2 q^2/beta_c^2]`, `beta_c = beta` numerically, and
the `:A_c` / `:S_c` reduced-amplitude classes (same strengths as `:A` / `:S`).

`recoil=true` applies the unequal-mass multiplier `m_c beta / ((m_c+m_d) beta_c)`
that the paper prints on the A_c P-wave rows (`K*_c`, and the `[D* pi]_D` rows of
`Q1c`/`Q2c`); the S_c S-wave rows and the 1^3S_1 A_c qtilde rows do NOT carry it.
"""
function charm_decay_amplitude(
    model::StrongDecayModel,
    coefficient::Real,
    class::Symbol,
    qbar_power::Integer,
    q_GeV::Real;
    recoil::Bool = false,
    m_c_GeV::Real = CHARM_M_C_GEV,
    m_d_GeV::Real = CHARM_M_D_GEV,
    beta_c_GeV::Real = model.beta_GeV,
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / beta_c_GeV
    recoil_mult = recoil ? m_c_GeV * model.beta_GeV / ((m_c_GeV + m_d_GeV) * beta_c_GeV) : 1.0
    return coefficient *
           reduced_decay_amplitude(model, class, qbar) *
           qbar^qbar_power *
           recoil_mult *
           _charm_suppressed_factor(q_GeV, beta_c_GeV, m_c_GeV, m_d_GeV)
end

"""
    strong_decay_amplitude(model, coefficient, class, qbar_power, q_GeV)

Full Table V amplitude in `MeV^(1/2)`: `coefficient` is the signed flavor/spin
factor from the formula column, `class` the reduced-amplitude class,
`qbar_power` the explicit `qbar^L` power, `q_GeV` the breakup momentum.
"""
function strong_decay_amplitude(
    model::StrongDecayModel,
    coefficient::Real,
    class::Symbol,
    qbar_power::Integer,
    q_GeV::Real,
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / model.beta_GeV
    return coefficient *
           reduced_decay_amplitude(model, class, qbar) *
           qbar^qbar_power *
           _suppressed_factor(q_GeV, model.beta_GeV)
end

"""
    calibrate_strong_decay_model(rho_q_GeV, B_q_GeV; rho_amplitude, B_amplitude, beta_GeV)

Fix `A` from `rho -> pi pi` (`+(4/3)^(1/2) A qbar`, paper `+12.4 MeV^(1/2)`)
and then `S0` from `B -> [omega pi]_S` (`-(2/9)^(1/2) S(qbar)`, paper `-11`),
given the breakup momenta of the two fit decays.
"""
function calibrate_strong_decay_model(
    rho_q_GeV::Real,
    B_q_GeV::Real;
    rho_amplitude::Real = 12.4,
    B_amplitude::Real = -11.0,
    beta_GeV::Real = 0.40,
)
    qbar_rho = rho_q_GeV / beta_GeV
    A = rho_amplitude /
        (sqrt(4 / 3) * qbar_rho * _suppressed_factor(rho_q_GeV, beta_GeV))
    qbar_B = B_q_GeV / beta_GeV
    S_at_B = B_amplitude / (-sqrt(2 / 9) * _suppressed_factor(B_q_GeV, beta_GeV))
    S0 = S_at_B + 0.5 * A * qbar_B^2
    return StrongDecayModel(A, S0, beta_GeV)
end
