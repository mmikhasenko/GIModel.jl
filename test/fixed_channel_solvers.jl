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

    # Empty is still the right answer where the path is genuinely inactive:
    # not an S wave, or a multiplicity the contact term does not touch.
    @test contact_hyperfine_nonperturbative_states(params, masses, "P", 1, r, 2) === nothing
    @test GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 2, r, 2) == Float64[]

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

@testset "RadialSolver is numerics, SpinTerms is physics" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    meson = Meson(mq, :q, :s)          # unequal mass: every mixing mechanism is live
    levels = spectrum_levels(2)
    base = compute_spectrum(params, meson; levels = levels)
    masses(spec) = [s.mass_GeV for s in spec.states]

    # The struct defaults ARE the old inline defaults: passing them explicitly
    # must reproduce the default call bit for bit.
    @test masses(compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 450, rmax = 24.0, kinetic = :relativistic,
            eigensolver = :full, nlevels_per_channel = 6),
        terms = SpinTerms(contact_hyperfine = true, fine_structure = true,
            same_j_spin_orbit = true, tensor = true))) == masses(base)
    @test masses(compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(), terms = SpinTerms())) == masses(base)

    # RadialSolver is "how well": refining the mesh must not move a mass more
    # than the discretization error it removes (sub-MeV here).
    fine = compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 900, rmax = 32.0))
    @test maximum(abs.(masses(fine) .- masses(base))) < 0.003   # < 3 MeV

    # SpinTerms is "what physics": every switch must move a mass, and turning
    # them all off must land exactly on the central eigenvalues.
    for sw in (:contact_hyperfine, :fine_structure, :same_j_spin_orbit, :tensor)
        off = compute_spectrum(params, meson; levels = levels,
            terms = SpinTerms(; sw => false))
        @test maximum(abs.(masses(off) .- masses(base))) > 1e-4   # > 0.1 MeV
    end
    none = compute_spectrum(params, meson; levels = levels,
        terms = SpinTerms(contact_hyperfine = false, fine_structure = false,
            same_j_spin_orbit = false, tensor = false))
    @test all(s.mass_GeV === s.central_GeV for s in none.states)

    # The solver threads down to the low tier, and its copy constructor keeps
    # the untouched fields.
    sol = channel_solution(params, meson.constituent_masses, 0;
        solver = FiniteDifferenceSolver(nlevels_per_channel = 3))
    @test length(sol.eigenvalues_GeV) == 3
    tuned = FiniteDifferenceSolver(FiniteDifferenceSolver(ngrid = 900); rmax = 32.0)
    @test tuned.ngrid == 900 && tuned.rmax == 32.0 && tuned.kinetic === :relativistic

    # Every path is silent. There is nothing left to deprecate: the loose
    # keywords are gone, so the shim that used to fold them in -- and warn that
    # they silently overrode the objects -- is gone with them. A settings name
    # that no longer exists is now a MethodError at the call, not a mass that
    # quietly came from the wrong grid.
    @test_logs compute_spectrum(params, meson; levels = levels)
    @test_logs compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(), terms = SpinTerms())
    @test_logs central_spectrum(params, meson; levels = levels)
    @test_logs channel_solution(params, meson.constituent_masses, 0;
        solver = FiniteDifferenceSolver())
    # Every retired keyword, by its old name, on the entry point that used to
    # accept it. NOTE: these calls are deliberately written the old way -- that
    # they no longer parse into a method IS the assertion -- so do not rewrite
    # them into the solver/terms form.
    corrected = fixed_spectrum(params, meson; levels = levels)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels, ngrid = 450)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels, rmax = 24.0)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        eigensolver = :krylov)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        nlevels_per_channel = 3)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        tensor_mixing = true)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        contact_hyperfine = true)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        same_j_spin_orbit_mixing = true)
    @test_throws MethodError central_spectrum(params, meson; levels = levels, rmax = 24.0)
    @test_throws MethodError channel_solution(params, meson.constituent_masses, 0;
        kinetic = :relativistic)
    @test_throws MethodError solve_sector(params, 1.628; ngrid = 100)
    @test_throws MethodError fixed_spectrum(
        params, meson; levels = levels, use_fine_structure = false,
    )
    @test_throws MethodError add_intra_meson_mixing(corrected; tensor_mixing = false)

    # Invalid settings are construction errors, not silent fallbacks.
    @test_throws ArgumentError RadialSolver(kinetic = :newtonian)
    @test_throws ArgumentError RadialSolver(eigensolver = :lanczos)
    @test_throws ArgumentError RadialSolver(ngrid = 1)
    @test_throws ArgumentError RadialSolver(rmax = 0)
    @test_throws ArgumentError RadialSolver(nlevels_per_channel = 0)
end

@testset "The solver, not the parameters, chooses the radial method" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    meson = Meson(mq, :c, :c)
    levels = spectrum_levels(2)
    masses(spec) = [s.mass_GeV for s in spec.states]

    # `GIParameters` is the model and nothing else: no basis marker, and one
    # object serves both methods. It used to carry a phantom type parameter that
    # no field used, so switching numerical method meant rebuilding the physics.
    @test isconcretetype(typeof(params))
    @test !any(f -> occursin("asis", String(f)), fieldnames(GIParameters))

    # `RadialSolver` is the interface; calling it builds the default
    # implementation, which is what every pre-split call site meant.
    @test RadialSolver isa Type && !isconcretetype(RadialSolver)
    @test RadialSolver(ngrid = 900) === FiniteDifferenceSolver(ngrid = 900)
    @test FiniteDifferenceSolver() isa RadialSolver
    @test OscillatorSolver() isa RadialSolver

    # Same parameters, two methods, through the front door. Both must answer,
    # and agree to a few MeV -- the residual is each one's own truncation error.
    fd = compute_spectrum(params, meson; levels = levels, solver = FiniteDifferenceSolver())
    ho = compute_spectrum(params, meson; levels = levels, solver = OscillatorSolver())
    @test isconcretetype(typeof(meson))
    @test isconcretetype(typeof(fd.computation))
    @test isconcretetype(typeof(fd))
    @test maximum(abs.(masses(fd) .- masses(ho))) < 0.005

    # The oscillator path is variational in a finite basis, so it can only sit
    # ABOVE the converged answer. Against a fine finite-difference mesh, every
    # charmonium level must therefore come out no lower.
    fine = compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 1800, rmax = 24.0))
    @test all(masses(ho) .>= masses(fine) .- 1e-9)

    # Raising nbasis can only lower an oscillator eigenvalue, for the same reason.
    e24 = channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(nbasis = 24, converge = false),
        nlevels = 3).eigenvalues_GeV
    e40 = channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(nbasis = 40, converge = false),
        nlevels = 3).eigenvalues_GeV
    @test all(e40 .<= e24 .+ 1e-9)

    # The stage-1 solver is recorded and stages 2-3 reuse it, so one spectrum is
    # one calculation. The non-perturbative contact solve in stage 2 resolves its
    # own eigenproblem and would otherwise silently revert to finite differences.
    @test compute_spectrum(params, meson; levels = levels,
        solver = OscillatorSolver()).computation.solver isa OscillatorSolver

    # Settings the oscillator path does not have are absent, not ignored. The
    # `:nonrelativistic` comparator is a finite-difference thing, and asking an
    # OscillatorSolver for it fails at construction rather than silently doing
    # the relativistic calculation you did not ask for.
    @test !hasfield(OscillatorSolver, :kinetic)
    @test !hasfield(OscillatorSolver, :eigensolver)
    @test_throws MethodError OscillatorSolver(kinetic = :nonrelativistic)
    @test_throws MethodError OscillatorSolver(eigensolver = :krylov)
    @test_throws ArgumentError OscillatorSolver(nbasis = 0)
    @test_throws ArgumentError OscillatorSolver(max_nbasis = 1)
    @test_throws ArgumentError OscillatorSolver(basis_step = 0)
    @test_throws ArgumentError OscillatorSolver(energy_tolerance_GeV = 0)
    @test_throws ArgumentError OscillatorSolver(beta_tolerance_GeV = 0)
    @test_throws ArgumentError OscillatorSolver(beta_grid = Float64[])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [0.5, 0.6])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [0.5, -0.5])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [1.5, 0.5])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [0.5, 0.5, 0.6])

    # The paper's complete fixed-sector Hamiltonian is relativistic. The
    # nonrelativistic FD option remains a central-only comparator and must not
    # be silently upgraded to a relativistic production calculation.
    @test_throws ArgumentError compute_spectrum(
        params, meson; levels = levels,
        solver = FiniteDifferenceSolver(kinetic = :nonrelativistic),
    )

    # beta_grid is a solver field, not a module constant to edit in source. A
    # bracket that excludes the minimum now fails rather than returning a
    # warned-but-usable under-resolved result.
    narrow = OscillatorSolver(beta_grid = 0.9:0.1:1.2)
    @test narrow.beta_grid == [0.9, 1.0, 1.1, 1.2]
    @test_throws ErrorException channel_solution(
        params, meson.constituent_masses, 0; solver = narrow, nlevels = 2)

    @testset "PA-12 adaptive HO convergence is certified or fails loudly" begin
        # Explicit fixed-size mode is available only as an unchecked convergence
        # study. Production mode records the continuously refined beta and final
        # basis certificate on the existing solution object.
        unchecked = channel_solution(
            params, meson.constituent_masses, 0;
            solver = OscillatorSolver(nbasis = 24, converge = false),
            nlevels = 2,
        )
        @test unchecked.convergence.status == :unchecked
        @test isnothing(unchecked.convergence.energy_delta_GeV)
        certified = channel_solution(
            params, meson.constituent_masses, 0;
            solver = OscillatorSolver(), nlevels = 2,
        )
        @test certified.convergence.status == :converged
        @test certified.convergence.beta_GeV == radial_wave(certified, 1).beta
        @test certified.convergence.beta_GeV ∉ OscillatorSolver().beta_grid
        @test certified.convergence.energy_delta_GeV <=
              certified.convergence.tolerance_GeV
        high_order_left = OscillatorWave(0, 0.55, [sin(i) for i in 1:64])
        high_order_right = OscillatorWave(0, 0.65, [cos(i / 2) for i in 1:72])
        @test isfinite(radial_overlap(high_order_left, high_order_right, _ -> 1.0))
        high_order_momentum = momentum_wave(high_order_left)
        @test isfinite(momentum_functional(
            high_order_momentum, p -> inv(sqrt(1 + p^2)),
        ))
        @test isfinite(momentum_overlap(
            high_order_momentum, momentum_wave(high_order_right), _ -> 1.0,
        ))
        too_small = OscillatorSolver(
            nbasis = 8,
            max_nbasis = 16,
            basis_step = 4,
            energy_tolerance_GeV = 1e-14,
            beta_grid = [0.8],
            nlevels_per_channel = 1,
        )
        @test_throws ErrorException channel_solution(
            params, meson.constituent_masses, 0; solver = too_small, nlevels = 1)
        @test_throws ArgumentError channel_solution(
            params,
            meson.constituent_masses,
            0;
            solver = OscillatorSolver(nbasis = 1, converge = false),
            nlevels = 2,
        )
    end
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
