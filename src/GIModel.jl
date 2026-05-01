module GIModel

using LinearAlgebra
using Printf
using TOML
using CSV
using KrylovKit: eigsolve
using SpecialFunctions: erf

include("constants.jl")

export GIParameters, ReferenceState, load_parameters, load_reference_spectrum
include("parameters.jl")

include("running_coupling.jl")
include("smearing_appendix_a.jl")
include("radial_1d_coulomb_smear.jl")
include("appendix_a_derivative_potential.jl")

export reduced_mass
include("radial_grid.jl")

export central_potential_mode, central_potential_values
include("central_potential_dispatch.jl")

include("hamiltonian.jl")

export channel_solution
include("channel_solver.jl")

export spin_dot
include("contact_hyperfine.jl")

export parse_quark_masses
include("masses_from_content.jl")

export fine_structure_split, LdotS, tensor_triplet_LJ
include("spin_fine_structure.jl")

export CentralPotentialPath, central_potential_path
include("appendix_a_status.jl")

export solve_sector, compare_sector, write_residual_report
include("sector_workflow.jl")

end
