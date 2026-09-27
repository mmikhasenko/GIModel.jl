# Normalization of u(r) across every solver and orbital sector (~17 s).

@testset "Every solve returns u with the same normalization" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    # The central solve, with both solvers and every orbital: integral u^2 dr = 1.
    # FD used to return Euclidean eigenvectors here (sum u^2 = 1) while HO
    # returned physical ones, differing by sqrt(h).
    for key in ("q", "s", "c", "b"), L in 0:2
        masses = ConstituentMasses(mq[key], mq[key])
        for slv in (FiniteDifferenceSolver(), OscillatorSolver())
            sol = channel_solution(params, masses, L; solver = slv, nlevels = 2)
            for wave in sol.waves
                @test isapprox(wave_norm(wave), 1.0; rtol = 1e-10)
            end
            if slv isa OscillatorSolver
                @test sol.convergence.status == :converged
                @test sol.convergence.energy_delta_GeV <= slv.energy_tolerance_GeV
                @test sol.convergence.refinements >=
                      GIModel.HO_REQUIRED_CONVERGED_REFINEMENTS
                @test sol.convergence.max_wave_overlap_defect >= 0
            else
                @test isnothing(sol.convergence)
            end
        end
    end

    # The convention has one definition, and it is idempotent.
    raw = reshape(collect(1.0:10.0), 10, 1)
    once = GIModel.physically_normalized_waves(raw, 0.25)
    @test isapprox(sum(abs2, once) * 0.25, 1.0; rtol = 1e-12)
    @test GIModel.physically_normalized_waves(once, 0.25) ≈ once
    # A zero column is left alone rather than producing NaN.
    @test all(iszero, GIModel.physically_normalized_waves(zeros(5, 1), 0.25))
end
