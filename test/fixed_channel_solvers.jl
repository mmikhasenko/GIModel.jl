@testset "Native spin matrices equal wave-interface expectations" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["c"], mq["c"])
    beta, nbasis = 0.75, 12

    cS = collect(range(1.0, 0.1; length = nbasis))
    wS = OscillatorWave(0, beta, cS)
    contact = GIModel.ho_contact_hyperfine_matrix(
        params, masses, 0, 1, beta, nbasis,
    )
    @test Matrix(contact) ≈ Matrix(contact)' atol = 1e-13
    @test dot(wS.coefficients, contact * wS.coefficients) ≈
          GIModel.contact_hyperfine_shift_momentum_sandwich(
              params, masses, FineStructureMultiplet("S", 1, 0), wS,
          ) rtol = 1e-5 atol = 1e-9

    cP = [sin(i) + 0.2cos(2i) for i in 1:nbasis]
    wP = OscillatorWave(1, beta, cP)
    matrices = GIModel.ho_fine_structure_matrices(
        params, masses, 1, 3, 2, beta, nbasis,
    )
    components = fine_structure_components(
        params, masses, FineStructureMultiplet("P", 3, 2), wP,
    )
    @test all(
        M -> isapprox(Matrix(M), Matrix(M)'; atol = 1e-13),
        (matrices.spin_orbit_vector, matrices.spin_orbit_thomas,
         matrices.spin_orbit, matrices.tensor, matrices.total),
    )
    @test dot(wP.coefficients, matrices.spin_orbit_vector * wP.coefficients) ≈
          components.spin_orbit_vector rtol = 1e-5 atol = 1e-9
    @test dot(wP.coefficients, matrices.spin_orbit_thomas * wP.coefficients) ≈
          components.spin_orbit_thomas rtol = 1e-5 atol = 1e-9
    @test dot(wP.coefficients, matrices.tensor * wP.coefficients) ≈
          components.tensor rtol = 1e-5 atol = 1e-9
    @test Matrix(matrices.total) ≈
          Matrix(matrices.spin_orbit) + Matrix(matrices.tensor) atol = 1e-13

    r, h = GIModel.radial_grid(120, 12.0)
    wfd = MeshWave(r .^ 2 .* exp.(-r), r)
    fdmat = GIModel.fine_structure_grid_matrices(params, masses, 2, r, h; L = 1)
    fdcomp = fine_structure_components(
        params, masses, FineStructureMultiplet("P", 3, 2), wfd,
    )
    expect(M) = h * dot(wfd.u, M * wfd.u) / wave_norm(wfd)
    @test expect(fdmat.spin_orbit_vector) ≈ fdcomp.spin_orbit_vector atol = 1e-10
    @test expect(fdmat.spin_orbit_thomas) ≈ fdcomp.spin_orbit_thomas atol = 1e-10
    @test expect(fdmat.tensor) ≈ fdcomp.tensor atol = 1e-10
    @test Matrix(fdmat.total) ≈
          Matrix(fdmat.spin_orbit) + Matrix(fdmat.tensor) atol = 1e-13
end

@testset "Both solvers solve the same native fixed-(L,S,J) Hamiltonian" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r, h = GIModel.radial_grid(450, 24.0)

    for key in ("q", "c", "b")
        masses = ConstituentMasses(mq[key], mq[key])
        multiplet = FineStructureMultiplet("S", 1, 0)
        terms = SpinTerms(fine_structure = false)
        sol_fd = fixed_channel_solution(
            params, masses, multiplet;
            solver = FiniteDifferenceSolver(), terms = terms, nlevels = 2)
        sol_ho = fixed_channel_solution(
            params, masses, multiplet;
            solver = OscillatorSolver(), terms = terms, nlevels = 2)

        # Same normalization, independent of representation.
        # FD used to return Euclidean eigenvectors (sum u^2 = 1), differing from
        # HO by exactly sqrt(h); anything quadratic in u was then off by h.
        for wave in vcat(sol_fd.waves, sol_ho.waves)
            @test isapprox(wave_norm(wave), 1.0; rtol = 1e-10)
        end
        # Same quantity, to the combined accuracy of two now-INDEPENDENT methods.
        # This was 1e-3 while the oscillator path projected the finite-difference
        # p^2 through the mesh: it was then a Galerkin restriction of the FD
        # problem, inherited FD's discretization error, and so agreed artificially
        # well. With Eq. (A17) assembled from exact matrix elements the two share
        # no error, and the gap measures both: +0.18 MeV (charm) and +0.55 MeV
        # (bottom) with the oscillator answer correctly ABOVE — variational in a
        # finite basis — and ~1.6 MeV for light quarks, where the finite-difference
        # mesh is itself least converged and sits above the true value.
        @test abs(sol_ho.eigenvalues_GeV[1] - sol_fd.eigenvalues_GeV[1]) < 3e-3
    end

    # A mesh operator is not an HO input. Unsupported routes fail explicitly.
    masses = ConstituentMasses(mq["c"], mq["c"])
    V = Matrix(GIModel.contact_hyperfine_operator(params, masses, "S", 1, r))
    ho = OscillatorSolver()
    @test_throws ArgumentError resummed_channel_solution(
        params, masses, 0, V; solver = ho, nlevels = 2)

    # V must live on the solver's mesh, whichever solver that is.
    bad = zeros(10, 10)
    @test_throws Exception resummed_channel_solution(
        params, masses, 0, bad; solver = FiniteDifferenceSolver(), nlevels = 1)
    @test_throws Exception resummed_channel_solution(
        params, masses, 0, bad; solver = ho, nlevels = 1)
end

@testset "An unimplemented solve fails instead of degrading" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["q"], mq["q"])
    r, _ = GIModel.radial_grid(450, 24.0)

    # FD has the non-perturbative contact solve, and returns it.
    solution = contact_hyperfine_nonperturbative_states(params, masses, "S", 1, r, 2)
    @test length(solution.eigenvalues_GeV) == 2 && length(solution.waves) == 2

    # All partial waves use the fixed-sector solver; invalid input is an error.
    @test length(contact_hyperfine_nonperturbative_states(
        params, masses, "P", 1, r, 2).waves) == 2
    @test_throws ArgumentError GIModel.contact_hyperfine_nonperturbative_levels(
        params, masses, "S", 2, r, 2)

    # The oscillator basis used to answer "empty" here, which
    # the old spectrum correction stage read as "fall back to first-order PT": the light
    # 1S0 came out 0.2842 GeV against the correctly resummed 0.149 GeV. U1 made that a loud
    # failure, U2 unified the solve, and B1 made this wrapper basis-generic — so
    # the oscillator path now resums in its own space and lands on the same
    # answer. That agreement is what the throw was standing in for.
    ho = OscillatorSolver()
    solution_ho = contact_hyperfine_nonperturbative_states(
        params, masses, "S", 1, 2; solver = ho)
    @test length(solution_ho.eigenvalues_GeV) == 2 && length(solution_ho.waves) == 2
    @test abs(solution_ho.eigenvalues_GeV[1] - solution.eigenvalues_GeV[1]) < 1e-3
    @test 0.12 < solution_ho.eigenvalues_GeV[1] < 0.18
    @test GIModel.contact_hyperfine_nonperturbative_levels(
        params, masses, "S", 1, 2; solver = ho) ≈ solution_ho.eigenvalues_GeV
    @test_throws ArgumentError contact_hyperfine_nonperturbative_states(
        params, masses, "S", 1, r, 2; solver = ho)
    ho_pi = spectrum_state(
        compute_spectrum(params, Meson(mq, :q, :q); solver = ho, levels = spectrum_levels(1)),
        "1^1S_0")
    @test 0.12 < ho_pi.mass_GeV < 0.18

    # The FD front door is untouched and still resums.
    fd_pi = spectrum_state(
        compute_spectrum(params, Meson(mq, :q, :q); levels = spectrum_levels(1)), "1^1S_0")
    @test 0.12 < fd_pi.mass_GeV < 0.18 # resummed; first-order PT lands near 0.28
end

@testset "Resolution walls are detected, not silent" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    # Grid points spanned by the ground state's RMS radius. The FD grid is fixed,
    # so a compact enough state falls between points and comes back silently
    # under-resolved (at 30 GeV the hyperfine splitting collapses to 0.0025).
    function points_across(m; ngrid = 450, rmax = 24.0)
        sol = channel_solution(
            params, ConstituentMasses(m, m), 0;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax),
        )
        wave = radial_wave(sol, 1)
        rms = sqrt(radial_expect(wave, x -> x^2))
        return rms / wave.h
    end
    # Every sector the paper actually uses must sit clear of the threshold.
    for m in (mq["q"], mq["s"], mq["c"], mq["b"])
        @test points_across(m) > GIModel.MIN_POINTS_ACROSS_STATE
    end
    # ...including bottomonium on the deliberately coarser grid the Table III
    # mixing audit uses (ngrid = 220, rmax = 22): 10.3 points is ~2 MeV of
    # error, imprecise but not broken, so it must NOT warn.
    @test points_across(mq["b"]; ngrid = 220, rmax = 22.0) > GIModel.MIN_POINTS_ACROSS_STATE
    # The genuinely broken case must be caught: at 6.9 points the hyperfine
    # splitting collapses.
    @test points_across(30.0) < GIModel.MIN_POINTS_ACROSS_STATE

    # Same story on the oscillator path: the default beta_grid rails only well
    # above bottomonium. The adaptive production solve rejects that bracket.
    default_grid = OscillatorSolver().beta_grid
    function best_beta(m)
        mm = ConstituentMasses(m, m)
        r, h = GIModel.radial_grid(450, 24.0)
        best = nothing
        for b in default_grid
            H, _ = GIModel.oscillator_hamiltonian_for_beta(params, mm, 0, r, h, b; nbasis = 24)
            v, _ = GIModel.lowest_eigenpairs(Matrix(H), 2, Val(:full))
            (isnothing(best) || v[end] < best[2]) && (best = (b, v[end]))
        end
        return best[1]
    end
    @test best_beta(mq["b"]) < last(default_grid)     # bottomonium is fine
    @test best_beta(30.0) == last(default_grid)       # railed

    @test_throws ErrorException channel_solution(
        params, ConstituentMasses(30.0, 30.0), 0;
        solver = OscillatorSolver(), nlevels = 2)
end
