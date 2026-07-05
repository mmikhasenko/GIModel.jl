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

export Meson, is_equal_flavor, flavor_label
include("meson.jl")

# =============================================================================
# Numerics: potentials, Hamiltonian, radial solves, spin-dependent corrections
# =============================================================================

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

export spin_dot
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
    spin_orbit_mixing_components,
    tensor_mixing_components,
    same_j_mixing,
    LdotS,
    tensor_triplet_LJ,
    tensor_triplet_offdiag_sameJ
include("spin_fine_structure.jl")

export CentralPotentialPath, central_potential_path
include("appendix_a_status.jl")

export RadialChannelKey, ChannelRadialSolution, SectorComputation, solve_sector
include("sector_solver.jl")

export spectrum_levels,
    StateMixing, SpectrumState, Spectrum, compute_spectrum, spectrum_state, parameters
include("spectrum.jl")

export annihilation_basis_input, isoscalar_annihilation_block, pseudoscalar_annihilation_block
include("flavor_mixing.jl")

# =============================================================================
# Strong decays: Table IV/V amplitude model
# =============================================================================

export StrongDecayModel,
    decay_momentum,
    reduced_decay_amplitude,
    strong_decay_amplitude,
    calibrate_strong_decay_model
include("strong_decays.jl")

end
