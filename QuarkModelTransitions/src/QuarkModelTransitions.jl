module QuarkModelTransitions

using GIModel
using PartialWaveFunctions: CG
using Printf

include("flavor_algebra.jl")
include("spin_algebra.jl")
include("algebraic_decomposition.jl")

export TransitionOperator,
    StrongDecayOperator,
    PhysicalState,
    ReferenceState,
    physical_state,
    TwoMesonChannel,
    PartialWave,
    allowed_partial_waves,
    partial_wave_projection,
    OnShell,
    CMKinematics,
    RelativisticTwoBodyNormalization,
    GITableVNormalization,
    TransitionAmplitude,
    matrix_element,
    partial_waves
include("transition_amplitudes.jl")

export StrongDecayModel,
    TableIVPolynomial,
    LeadingS0,
    TableVReference,
    decay_momentum,
    reduced_decay_amplitude,
    spatial_overlap,
    strong_decay_amplitude,
    calibrate_strong_decay_model,
    DecayChannel,
    StrongDecayAmplitude,
    decay_amplitude,
    reduced_matrix_element,
    decay_width,
    MesonMasses,
    meson_mass
include("strong_decays.jl")

end
