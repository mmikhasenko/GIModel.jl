# Solver-convergence sweeps across flavors, solvers and adaptive HO bases.
# ~90 s; run with GI_HEAVY_TESTS=true (scripts/verify_project.sh does).

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
    # The charmonium S-wave optimum is near 1.05 GeV, so a bracket on either side
    # rails at its endpoint. (A 0.9:0.1:1.2 bracket used to rail too, but only
    # because the uncorrected momentum projection made the optimum grid-dependent.)
    narrow = OscillatorSolver(beta_grid = 0.5:0.1:0.8)
    @test narrow.beta_grid == [0.5, 0.6, 0.7, 0.8]
    @test_throws ErrorException channel_solution(
        params, meson.constituent_masses, 0; solver = narrow, nlevels = 2)
    @test_throws ErrorException channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(beta_grid = 1.5:0.1:1.8), nlevels = 2)
    # A bracket that contains the optimum finds the same β as the default grid.
    @test channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(beta_grid = 0.9:0.1:1.2), nlevels = 2).convergence.beta_GeV ≈
          channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(), nlevels = 2).convergence.beta_GeV atol = 5e-3

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
