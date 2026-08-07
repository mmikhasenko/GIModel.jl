using GIModel
using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PARAMETER_FILE = joinpath(ROOT, "data", "parameters.provisional.toml")

function measure_warm(label, f; samples = 7)
    f() # compilation/warm-up is deliberately outside the samples
    observations = map(1:samples) do _
        GC.gc()
        @timed f()
    end
    times = sort(getproperty.(observations, :time))
    bytes = sort(getproperty.(observations, :bytes))
    @printf("%-36s median=%8.4f s  min=%8.4f s  min_alloc=%10d bytes\n",
        label, times[cld(samples, 2)], first(times), first(bytes))
    return nothing
end

function main()
    params, mass_table = load_parameters_and_quark_masses(PARAMETER_FILE)
    masses = ConstituentMasses(mass_table["c"], mass_table["c"])
    fd_solver = FiniteDifferenceSolver(
        ngrid = 450, rmax = 24.0, nlevels_per_channel = 3,
    )
    ho_solver = OscillatorSolver(
        ngrid = 450, rmax = 24.0, nlevels_per_channel = 3,
    )

    fd_matrix() = Matrix(first(GIModel.relativistic_hamiltonian(
        params, masses, 0; solver = fd_solver,
    )))
    fd_values() = eigvals(Symmetric(fd_matrix()))[1:3]

    beta = 0.65
    nbasis = 24
    r, h = GIModel.radial_grid(ho_solver.ngrid, ho_solver.rmax)
    ho_matrix() = Matrix(first(GIModel.oscillator_hamiltonian_for_beta(
        params, masses, 0, r, h, beta; nbasis = nbasis,
    )))
    ho_values() = eigvals(Symmetric(ho_matrix()))[1:3]
    ho_selected_values() = channel_solution(
        params, masses, 0; solver = ho_solver,
    ).eigenvalues_GeV

    println("audit baseline: charmonium, L=0, active Appendix-A central path")
    println("Julia: ", VERSION)
    println("FD context: ngrid=", fd_solver.ngrid, ", rmax=", fd_solver.rmax)
    println("HO context: beta=", beta, ", nbasis=", nbasis,
        " (selected route retains native OscillatorWave coefficients)")
    println()
    println("FD central values [GeV]: ",
        join((@sprintf("%.10f", x) for x in fd_values()), ", "))
    println("HO fixed-beta values [GeV]: ",
        join((@sprintf("%.10f", x) for x in ho_values()), ", "))
    println("HO selected-route values [GeV]: ",
        join((@sprintf("%.10f", x) for x in ho_selected_values()), ", "))
    println()
    measure_warm("FD assemble H (450x450)", fd_matrix)
    measure_warm("FD assemble H + eigvals", fd_values)
    measure_warm("HO fixed-beta H (24x24)", ho_matrix)
    measure_warm("HO fixed-beta H + eigvals", ho_values)
    measure_warm("HO selected 22-beta route", ho_selected_values)

    if isdefined(GIModel, :_GAUSS_LAGUERRE_DVR)
        cache = getfield(GIModel, :_GAUSS_LAGUERRE_DVR)
        println()
        println("worktree HO DVR cache entries after measurements: ", length(cache))
        println("The cache is process-global mutable orchestration ",
            "and is not an active AD input.")
    end
    return nothing
end

main()
