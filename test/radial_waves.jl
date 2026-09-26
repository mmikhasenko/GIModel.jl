@testset "RadialWave interface: invariants, not fixed numbers" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r, h = GIModel.radial_grid(450, 24.0)

    @test MeshWave <: RadialWave

    # --- closed-form cross-checks: the operations must reproduce integrals we
    # --- can do by hand, independently of any solver.
    # u = r e^{-r}  =>  <r^n> = Gamma(3+n)/(2^n Gamma(3)) exactly.
    u = r .* exp.(-r)
    w = MeshWave(u ./ sqrt(sum(abs2, u) * h), r)
    @test isapprox(wave_norm(w), 1.0; rtol = 1e-10)
    @test isapprox(radial_expect(w, x -> 1.0), 1.0; rtol = 1e-10)
    @test isapprox(radial_expect(w, x -> x), 1.5; rtol = 1e-5)          # 3/2
    @test isapprox(radial_expect(w, x -> x^2), 3.0; rtol = 1e-5)        # 3
    # A normalized wave overlapped with itself is 1, and the overlap is
    # invariant under rescaling either wave (it normalizes internally).
    @test isapprox(radial_overlap(w, w, x -> 1.0), 1.0; rtol = 1e-10)
    w_scaled = MeshWave(7.3 .* w.u, r)
    @test isapprox(radial_overlap(w, w_scaled, x -> 1.0), 1.0; rtol = 1e-10)
    @test radial_overlap(w, w, _ -> 1 + 2im) ≈ 1 + 2im
    @test isapprox(radial_expect(w_scaled, x -> x), radial_expect(w, x -> x); rtol = 1e-12)

    # Integration by parts gives integral r*u*u' dr = -1/2 for every
    # normalized reduced radial wave with vanishing endpoints.
    @test GIModel.radial_derivative_overlap(w, w, identity) ≈ -0.5 atol = 2e-3
    @test GIModel.radial_derivative_overlap(w, w, x -> (1 + im) * x) ≈
          -0.5(1 + im) atol = 3e-3

    ho_ground = OscillatorWave(0, 0.4, [1.0])
    @test GIModel.radial_derivative_overlap(ho_ground, ho_ground, identity) ≈
          -0.5 atol = 1e-10
    @test radial_overlap(ho_ground, ho_ground, _ -> 1 + 2im) ≈
          1 + 2im atol = 1e-10

    # Two different waves: overlap is symmetric and bounded by 1 (Cauchy-Schwarz).
    v = r .* exp.(-0.7 .* r)
    wv = MeshWave(v ./ sqrt(sum(abs2, v) * h), r)
    ov = radial_overlap(w, wv, x -> 1.0)
    @test isapprox(ov, radial_overlap(wv, w, x -> 1.0); rtol = 1e-12)
    @test abs(ov) <= 1.0 + 1e-10

    # --- Parseval: the transform must conserve probability. This ties the
    # --- position-space and momentum-space halves of the interface together.
    for key in ("q", "c", "b")
        masses = ConstituentMasses(mq[key], mq[key])
        sol = channel_solution(params, masses, 0; nlevels = 2)
        for n in 1:2
            wn = radial_wave(sol, n)
            mw = momentum_wave(wn, 0)
            @test isapprox(wave_norm(wn), 1.0; rtol = 1e-10)
            @test isapprox(momentum_expect(mw, p -> 1.0), 1.0; rtol = 1e-6)
            # <E> >= m for a relativistic quark energy, with equality only at p=0.
            @test momentum_expect(mw, p -> sqrt(mq[key]^2 + p^2)) > mq[key]
        end
    end

    # --- orthogonality: distinct radial levels of one channel are orthogonal.
    # --- Nothing in the solve enforces this; it is a property of the operator,
    # --- so it is a real check on the solve rather than on the interface.
    masses = ConstituentMasses(mq["c"], mq["c"])
    sol = channel_solution(params, masses, 0; nlevels = 3)
    for i in 1:3, j in 1:3
        w_i = radial_wave(sol, i)
        w_j = radial_wave(sol, j)
        target = i == j ? 1.0 : 0.0
        @test isapprox(abs(radial_overlap(w_i, w_j, x -> 1.0)), target; atol = 1e-8)
    end

    # --- node counting: the n-th radial level has n-1 interior nodes. Count only
    # --- where the wave is physically present: the exponential tail flips sign
    # --- in floating point (measured: crossings at r = 17.6 and 23.8 with |u|
    # --- at 2e-12 and 2e-16 of the peak), so a naive sign count reports phantom
    # --- nodes in the ground state.
    for n in 1:3
        u_n = radial_wave(sol, n).u
        cut = 1e-8 * maximum(abs, u_n)
        big = [i for i in eachindex(u_n) if abs(u_n[i]) > cut]
        signs = sign.(u_n[big])
        nodes = count(i -> signs[i] != signs[i+1], 1:(length(signs)-1))
        @test nodes == n - 1
    end

    # --- momentum space is its own noun, with a linear and a quadratic operation.
    # --- momentum_expect is quadratic in Phi, momentum_functional is linear;
    # --- conflating them is a whole class of factor-of-Phi bug.
    @test MeshMomentumWave <: MomentumWave
    mc = mq["c"]
    solq = channel_solution(params, ConstituentMasses(mc, mc), 0; nlevels = 1)
    wq = radial_wave(solq, 1)
    mwq = momentum_wave(wq, 0)
    # quadratic with g = 1 is the norm; linear with K = 1 is NOT (different object)
    @test isapprox(momentum_expect(mwq, p -> 1.0), 1.0; rtol = 1e-6)
    @test !isapprox(momentum_functional(mwq, p -> 1.0), 1.0; rtol = 1e-3)
    # both are linear in their kernel argument
    @test isapprox(momentum_expect(mwq, p -> 3.0), 3 * momentum_expect(mwq, p -> 1.0); rtol = 1e-12)
    @test isapprox(momentum_functional(mwq, p -> 3.0), 3 * momentum_functional(mwq, p -> 1.0); rtol = 1e-12)
    # the functional scales linearly in Phi, the expectation quadratically
    scaled = MeshMomentumWave(mwq.p, 2 .* mwq.phi)
    @test isapprox(momentum_functional(scaled, p -> 1.0), 2 * momentum_functional(mwq, p -> 1.0); rtol = 1e-12)
    @test isapprox(momentum_expect(scaled, p -> 1.0), 4 * momentum_expect(mwq, p -> 1.0); rtol = 1e-12)

    # --- scale invariance: every consumer of a wave must give the same answer
    # --- for a rescaled wave, because the amplitude of u carries no physics.
    # --- This is what catches a deleted normalization: migrating
    # --- charge_radius_squared onto `radial_expect` removed the normalization
    # --- that `_rel_momentum_average` was silently relying on, and only a
    # --- rescaled input revealed it (the value moved by a factor of 13.7^2).
    sol3 = channel_solution(params, ConstituentMasses(mq["q"], mq["s"]), 0; nlevels = 1)
    w_one = radial_wave(sol3, 1)
    w_big = MeshWave(13.7 .* w_one.u, w_one.r)
    for probe in (
        w -> radial_expect(w, x -> x^2),
        w -> wave_norm(w) / wave_norm(w),
        w -> momentum_expect(momentum_wave(w, 0), p -> 1.0),
    )
        @test isapprox(probe(w_one), probe(w_big); rtol = 1e-10)
    end

    # Mismatched meshes are an error, not a silently wrong overlap.
    other, _ = GIModel.radial_grid(200, 24.0)
    @test_throws ArgumentError radial_overlap(w, MeshWave(other .* 0 .+ 1.0, other), x -> 1.0)
end

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
