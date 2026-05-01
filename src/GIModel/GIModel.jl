module GIModel

using LinearAlgebra
using Printf
using TOML
using CSV
using KrylovKit: eigsolve
using SpecialFunctions: erf

export GIParameters,
    ReferenceState,
    load_parameters,
    load_reference_spectrum,
    solve_sector,
    compare_sector,
    write_residual_report,
    parse_quark_masses,
    reduced_mass,
    channel_solution,
    fine_structure_split,
    LdotS,
    tensor_triplet_LJ,
    spin_dot,
    central_potential_mode,
    central_potential_values,
    CentralPotentialPath,
    central_potential_path

include("constants.jl")
include("parameters.jl")
include("running_coupling.jl")
include("smearing_appendix_a.jl")
include("radial_1d_coulomb_smear.jl")
include("appendix_a_derivative_potential.jl")
include("radial_grid.jl")
include("central_potential_dispatch.jl")
include("hamiltonian.jl")
include("channel_solver.jl")
include("contact_hyperfine.jl")
include("masses_from_content.jl")
include("spin_fine_structure.jl")
include("appendix_a_status.jl")
include("sector_workflow.jl")

end
