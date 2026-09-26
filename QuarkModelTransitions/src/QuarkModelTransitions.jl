module QuarkModelTransitions

using GIModel
using PartialWaveFunctions: CG
using Printf

export ELECTROMAGNETIC_DEFAULTS, STRONG_DECAY_DEFAULTS
include("observable_inputs.jl")

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

export PseudoscalarEmission
include("pseudoscalar_emission.jl")

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

# =============================================================================
# Table VII annihilation widths: gluonic QQ̄ → gluons
# =============================================================================

export wavefunction_origin_smearing,
    observable_momentum_wave,
    gluonic_annihilation_amplitude,
    gluonic_annihilation_width,
    GLUONIC_CHANNELS,
    mock_meson_mass,
    leptonic_decay_factor,
    LEPTONIC_FACTOR_KINDS,
    two_photon_amplitude,
    charge_radius_squared,
    leptonic_pseudoscalar_width,
    dilepton_vector_width,
    axial_tau_width,
    G_FERMI_GEV,
    HBARC_FM2
include("annihilation_widths.jl")

export ALPHA_EM,
    NUCLEON_MASS_GEV,
    photon_momentum,
    m1_transition_moment,
    e1_transition_amplitude,
    mock_mean_energy,
    mock_wave_mass,
    mock_meson_overlap,
    mock_meson_radial_moment
include("mock_meson_overlaps.jl")
export e1_angular_coefficient, spin_flip_photon_amplitude
export m1_recoil_moment, photon_recoil_form_factor, m1_radiative_width, neutral_m1_charge
include("radiative_decays.jl")

export PhotonEmission, PhotonEmitter, RadiativeAmplitude
include("photon_emission.jl")

end
