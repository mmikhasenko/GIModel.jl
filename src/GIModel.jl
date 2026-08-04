module GIModel

using LinearAlgebra
using Printf
using TOML
using KrylovKit: eigsolve
using SpecialFunctions: erf, gamma

# =============================================================================
# Constants and core layout types (shared by setup bookkeeping and numerics)
# =============================================================================

include("constants.jl")

export ConstituentMasses, reduced_mass, FineStructureMultiplet, RadialWaveOnUniformMesh
include("model_objects.jl")

# =============================================================================
# Setup / bookkeeping: TOML parameters and flavor mass table
# =============================================================================

export GIParameters, GIBasis, FiniteDifferenceBasis, HarmonicOscillatorBasis, with_basis
export ConfinementPotential,
    RelativisticSmearing, RelativisticFactors, FineStructure, AnnihilationAmplitudes
export CentralPotentialMethod,
    PointwiseCentral,
    Coulomb1DSmearing,
    AppendixASmearing3D,
    AppendixADerivativeG,
    AppendixAClosedForm,
    AppendixAMomentumSandwich,
    central_potential_method
export load_parameters
include("parameters.jl")

export QuarkMassTable, load_quark_masses, load_parameters_and_quark_masses
include("quark_mass_table.jl")

export RadialSolver, SpinTerms
include("solver_options.jl")

export AbstractQuark, LightQuark, StrangeQuark, HeavyQuark
export charge, flavor_symbol, mass_GeV
include("quark.jl")

export Meson, is_equal_flavor, flavor_label
include("meson.jl")

# =============================================================================
# Numerics: potentials, Hamiltonian, radial solves, spin-dependent corrections
# =============================================================================

export alpha_s_q, alpha_s_r
include("running_coupling.jl")
include("smearing_appendix_a.jl")
include("radial_1d_coulomb_smear.jl")
include("appendix_a_derivative_potential.jl")

include("radial_grid.jl")

export central_potential_values
include("central_potential_dispatch.jl")

include("hamiltonian.jl")
include("harmonic_oscillator_basis.jl")

export channel_solution
include("channel_solver.jl")

export spin_dot, contact_smearing_sigma, resummed_channel_solution,
    contact_hyperfine_nonperturbative_states
include("contact_hyperfine.jl")

export MixingMechanism,
    AntisymmetricSpinOrbit,
    TensorMixing,
    IsoscalarAnnihilation,
    BasisState,
    MixingBlock,
    MixingResult,
    diagonalize_mixing_block
include("state_mixing.jl")

export CalibratedP1Annihilation,
    PaperP1Annihilation,
    PaperP2Annihilation,
    FDOriginP2Smearing,
    FDMomentumIntegralSmearing,
    pseudoscalar_annihilation_basis_input,
    fix_annihilation_phase!,
    isoscalar_pseudoscalar_annihilation_solution,
    isoscalar_general_annihilation_solution,
    isoscalar_general_s1_solution
include("pseudoscalar_annihilation.jl")

export fine_structure_split,
    fine_structure_components,
    fine_structure_grid_operator,
    ho_first_order_distorted_states,
    ho_full_distorted_states,
    spin_orbit_mixing_components,
    tensor_mixing_components,
    same_j_mixing,
    LdotS,
    tensor_triplet_LJ,
    tensor_triplet_offdiag_sameJ,
    radial_cross_expect_udr
include("spin_fine_structure.jl")

export CentralPotentialPath, central_potential_path
include("appendix_a_status.jl")

export RadialChannelKey, ChannelRadialSolution, SectorComputation, solve_sector
include("sector_solver.jl")

export radial_wave,
    spectrum_levels,
    StateMixing,
    CentralState,
    CorrectedState,
    MixedState,
    Spectrum,
    CentralSpectrum,
    CorrectedSpectrum,
    MixedSpectrum,
    central_spectrum,
    add_spin_corrections,
    add_intra_meson_mixing,
    compute_spectrum,
    spectrum_state,
    parameters
include("spectrum.jl")

export annihilation_basis_input, isoscalar_annihilation_block, pseudoscalar_annihilation_block
include("flavor_mixing.jl")

# =============================================================================
# Strong decays: Table IV/V amplitude model
# =============================================================================

export StrongDecayModel,
    decay_momentum,
    reduced_decay_amplitude,
    spatial_overlap,
    strong_decay_amplitude,
    calibrate_strong_decay_model,
    DecayChannel,
    StrongDecayAmplitude,
    decay_amplitude,
    matrix_element,
    decay_width,
    MesonMasses,
    meson_mass,
    load_table_v
include("strong_decays.jl")

# =============================================================================
# Table VII annihilation widths: gluonic QQ̄ → gluons
# =============================================================================

export wavefunction_origin_smearing,
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

export MockMomentumWave,
    ALPHA_EM,
    NUCLEON_MASS_GEV,
    photon_momentum,
    m1_transition_moment,
    e1_transition_amplitude,
    mock_momentum_wave,
    mock_mean_energy,
    mock_wave_mass,
    mock_meson_overlap,
    mock_meson_radial_moment
include("mock_meson_overlaps.jl")

end
