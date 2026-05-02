module GIModel

using LinearAlgebra
using Printf
using TOML
using CSV
using KrylovKit: eigsolve
using SpecialFunctions: erf

# =============================================================================
# Constants and core layout types (shared by setup bookkeeping and numerics)
# =============================================================================

include("constants.jl")

export ConstituentMasses, reduced_mass, FineStructureMultiplet, RadialWaveOnUniformMesh
include("model_objects.jl")

# =============================================================================
# Setup / bookkeeping: TOML tables, reference CSV rows, attaching masses to rows
# =============================================================================

export GIParameters, load_parameters
include("parameters.jl")

export QuarkMassTable, load_quark_masses, load_parameters_and_quark_masses
include("quark_mass_table.jl")

export ReferenceState, ReferenceStateWithMasses, load_reference_spectrum
include("reference_spectrum.jl")

export parse_quark_masses, resolve_constituent_masses, attach_constituent_masses
include("masses_from_content.jl")

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

export channel_solution
include("channel_solver.jl")

export spin_dot
include("contact_hyperfine.jl")

export fine_structure_split, fine_structure_components, LdotS, tensor_triplet_LJ
include("spin_fine_structure.jl")

export CentralPotentialPath, central_potential_path
include("appendix_a_status.jl")

export RadialChannelKey,
    ChannelRadialSolution,
    SectorComputation,
    solve_sector,
    compute_sector,
    compare,
    write_residual_report
include("sector_workflow.jl")

end
