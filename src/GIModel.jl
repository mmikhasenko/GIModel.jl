module GIModel

using LinearAlgebra
using Printf
using TOML
using CSV
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

export GIParameters, GIBasis, FiniteDifferenceBasis, HarmonicOscillatorBasis, load_parameters
include("parameters.jl")

export QuarkMassTable, load_quark_masses, load_parameters_and_quark_masses
include("quark_mass_table.jl")

# =============================================================================
# Numerics: potentials, Hamiltonian, radial solves, spin-dependent corrections
# =============================================================================

include("running_coupling.jl")
include("smearing_appendix_a.jl")
include("radial_1d_coulomb_smear.jl")
include("appendix_a_derivative_potential.jl")

include("radial_grid.jl")

export central_potential_mode, central_potential_values
include("central_potential_dispatch.jl")

include("hamiltonian.jl")
include("harmonic_oscillator_basis.jl")

export channel_solution
include("channel_solver.jl")

export spin_dot
include("contact_hyperfine.jl")

export BasisState, MixingBlock, MixingResult, diagonalize_mixing_block
include("state_mixing.jl")

export isoscalar_pseudoscalar_annihilation_solution
include("pseudoscalar_annihilation.jl")

export fine_structure_split,
    fine_structure_components,
    spin_orbit_mixing_components,
    same_j_mixing,
    LdotS,
    tensor_triplet_LJ
include("spin_fine_structure.jl")

export CentralPotentialPath, central_potential_path
include("appendix_a_status.jl")

export RadialChannelKey, ChannelRadialSolution, SectorComputation, solve_sector
include("sector_solver.jl")

# =============================================================================
# IO: reference catalog types, CSV loader, string→mass resolution, attach
# Sector batch solves, comparison vs reference rows, residual markdown
# =============================================================================

export ReferenceState, ReferenceStateWithMasses, load_reference_spectrum
include("reference_state.jl")

export compute_sector, compare, mixing_prone_state, nonmixing_deviation_summary, write_residual_report
include("sector_comparison.jl")

export parse_quark_masses, resolve_constituent_masses, attach_constituent_masses
include("masses_from_content.jl")

end
