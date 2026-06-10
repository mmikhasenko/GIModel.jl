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
"""
function reduced_decay_amplitude(model::StrongDecayModel, class::Symbol, qbar::Real)
    class in (:A, :Aprime, :Adoubleprime, :A0) && return model.A
    class === :S && return model.S0 - 0.5 * model.A * qbar^2
    class === :D && return model.S0 - 0.3 * model.A * qbar^2
    class === :P && return model.S0 - 0.75 * model.A * qbar^2
    throw(ArgumentError("unknown reduced-amplitude class `$class`"))
end

# The Table V caption suppresses (q/2pi)^(1/2) exp(-q^2/16 beta^2) from the
# formula column; q enters the square root in MeV so amplitudes come out in
# MeV^(1/2), matching the numeric column.
function _suppressed_factor(q_GeV::Real, beta_GeV::Real)
    return sqrt(1000q_GeV / (2π)) * exp(-q_GeV^2 / (16 * beta_GeV^2))
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
