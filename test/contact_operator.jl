@testset "One smeared contact operator across orbital sectors" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["c"], mq["u"])
    beta = 0.7
    r = collect(0.02:0.02:16.0)
    for L in 0:4
        label = orbital_label(L)
        # Independent analytic Gaussian integral on the normalized n=0 HO wave:
        # <delta_tau> = tau^3/pi^(3/2) * (beta²/(beta²+tau²))^(L+3/2).
        sigma = contact_smearing_sigma(params, masses)
        density = sum(zip(GIModel.ALPHA_COEFFS, GIModel.ALPHA_GAMMAS)) do (alpha, gamma)
            tau = inv(sqrt(inv(sigma^2) + inv(gamma^2)))
            alpha * tau^3 / pi^1.5 * (beta^2 / (beta^2 + tau^2))^(L + 1.5)
        end
        expected = (1 + params.factors.epsilon_c) * 32pi / (9masses.m1_GeV * masses.m2_GeV) * density
        ho = OscillatorWave(L, beta, [1.0])
        mesh = sample_wave(ho, r)
        for wave in (ho, mesh)
            singlet = GIModel.contact_hyperfine_shift(params, masses, FineStructureMultiplet(label, 1, L), wave)
            triplet = GIModel.contact_hyperfine_shift(params, masses, FineStructureMultiplet(label, 3, L + 1), wave)
            @test singlet ≈ -3expected / 4 rtol = 2e-7
            @test triplet ≈ expected / 4 rtol = 2e-7
            @test singlet < 0 < triplet
            @test singlet ≈ -3triplet rtol = 1e-12
            # The dressed operator also acts for every L, with the same spin algebra.
            dressed1 = GIModel.contact_hyperfine_shift_momentum_sandwich(params, masses, FineStructureMultiplet(label, 1, L), wave)
            dressed3 = GIModel.contact_hyperfine_shift_momentum_sandwich(params, masses, FineStructureMultiplet(label, 3, L + 1), wave)
            @test dressed1 < 0 < dressed3
            @test dressed1 ≈ -3dressed3 rtol = 1e-12
        end
        # Finite-difference and native HO diagnostics use exactly the production
        # contact-only solve, including centrifugal momentum in each sector.
        for solver in (FiniteDifferenceSolver(ngrid = 120), OscillatorSolver(nbasis = 12, converge = false))
            diagnostic = contact_hyperfine_nonperturbative_states(params, masses, label, 1, 1; solver)
            direct = fixed_channel_solution(params, masses, FineStructureMultiplet(label, 1, L);
                solver, nlevels = 1,
                terms = SpinTerms(contact_hyperfine = true, fine_structure = false,
                                  same_j_spin_orbit = false, tensor = false))
            @test diagnostic.eigenvalues_GeV ≈ direct.eigenvalues_GeV rtol = 1e-12
            wave = radial_wave(direct, 1)
            # Hellmann-Feynman on the actual assembled Hamiltonian checks the
            # reported contribution against energy response, not another wrapper.
            matrices = if wave isa MeshWave
                GIModel._fixed_channel_matrices(solver, params, masses,
                    FineStructureMultiplet(label, 1, L), nothing,
                    SpinTerms(fine_structure = false), 1)
            else
                GIModel._fixed_channel_matrices(solver, params, masses,
                    FineStructureMultiplet(label, 1, L), wave.beta,
                    SpinTerms(fine_structure = false), length(wave.coefficients))
            end
            step = 1e-4
            energy(t) = first(eigvals(Symmetric(Matrix(matrices.total) + t * Matrix(matrices.contact))))
            derivative = (energy(step) - energy(-step)) / (2step)
            contribution = GIModel.contact_hyperfine_shift_active(params, masses, FineStructureMultiplet(label, 1, L), wave)
            @test contribution ≈ derivative atol = 2e-8 rtol = 2e-6
        end
    end
    @test_throws ArgumentError GIModel.contact_hyperfine_shift_active(
        params, masses, FineStructureMultiplet("P", 1, 1), OscillatorWave(0, beta, [1.0]))
    @test_throws ArgumentError GIModel.ContactHyperfine(params, masses, 2)
end

@testset "Explicit local contact prescription and invalid grids" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    local_params = GIParameters(params; factors = RelativisticFactors(
        params.factors; contact_momentum_sandwich = false))
    masses = ConstituentMasses(mq["c"], mq["u"])
    for solver in (FiniteDifferenceSolver(ngrid = 100), OscillatorSolver(nbasis = 12, converge = false))
        sol = contact_hyperfine_nonperturbative_states(local_params, masses, "D", 3, 1; solver)
        wave = radial_wave(sol, 1)
        mult = FineStructureMultiplet("D", 3, 3)
        active = GIModel.contact_hyperfine_shift_active(local_params, masses, mult, wave)
        @test active > 0
        @test active ≈ GIModel.contact_hyperfine_shift(local_params, masses, mult, wave)
    end
    @test_throws ArgumentError contact_hyperfine_nonperturbative_states(
        params, masses, "P", 1, [0.1, 0.21, 0.3], 1)
    @test_throws ArgumentError contact_hyperfine_nonperturbative_states(
        params, masses, "P", 1, [0.2, 0.3, 0.4], 1)
end
