using GIModel
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PARAMETER_FILE = joinpath(ROOT, "data", "parameters.provisional.toml")

function show_return_type(label, f, argument_types)
    return_type = Core.Compiler.return_type(f, argument_types)
    println(
        rpad(label, 34), " => ", return_type,
        "  [concrete: ", isconcretetype(return_type), "]",
    )
end

params, mass_table = load_parameters_and_quark_masses(PARAMETER_FILE)
masses = ConstituentMasses(mass_table["c"], mass_table["c"])
fd_solver = FiniteDifferenceSolver(nlevels_per_channel = 3)
ho_solver = OscillatorSolver(nlevels_per_channel = 3)

println("parameter object")
println("  type:          ", typeof(params))
println("  concrete type: ", isconcretetype(typeof(params)))
println("  field types:   ", fieldtypes(typeof(params)))
println()
println("compiler return-type snapshot")
show_return_type(
    "potential_diagonal",
    GIModel.potential_diagonal,
    Tuple{typeof(params),Float64,Float64,Vector{Float64}},
)
show_return_type(
    "relativistic_hamiltonian",
    GIModel.relativistic_hamiltonian,
    Tuple{typeof(params),typeof(masses),Int},
)
show_return_type(
    "full lowest_eigenpairs",
    GIModel.lowest_eigenpairs,
    Tuple{LinearAlgebra.Symmetric{Float64,Matrix{Float64}},Int,typeof(fd_solver)},
)
show_return_type(
    "HO fixed-beta Hamiltonian",
    GIModel.oscillator_hamiltonian_for_beta,
    Tuple{typeof(params),typeof(masses),Int,Vector{Float64},Float64,Float64},
)
show_return_type(
    "FD channel implementation",
    GIModel._channel_solution,
    Tuple{typeof(fd_solver),typeof(params),typeof(masses),Int,Int},
)
show_return_type(
    "HO channel implementation",
    GIModel._channel_solution,
    Tuple{typeof(ho_solver),typeof(params),typeof(masses),Int,Int},
)

println()
println("Interpretation: the loaded GIParameters and both solver values carry their")
println("method choices in concrete types. The reported boundaries should remain")
println("concretely inferred as later numeric types are generalized.")
