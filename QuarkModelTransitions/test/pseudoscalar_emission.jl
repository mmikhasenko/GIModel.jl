@testset "GI Eq. (19) native S-wave spatial integrals" begin
    mq = QuarkMassTable("u" => 0.22, "d" => 0.22, "c" => 1.7)
    operator = PseudoscalarEmission(0.73, -0.28, mq)
    @test operator.g == 0.73
    @test operator.h == -0.28
    @test operator.quark_masses !== mq
    @test_throws ArgumentError PseudoscalarEmission(Inf, 0.0, mq)
    @test_throws ArgumentError PseudoscalarEmission(0.0, 0.0, QuarkMassTable("u" => 0.0))

    qtop = QuarkModelTransitions._QuarkEmission()
    atop = QuarkModelTransitions._AntiquarkEmission()
    @test QuarkModelTransitions._eq19_momentum_fraction(operator, qtop, (:u, :d)) == 0.5
    @test QuarkModelTransitions._eq19_momentum_fraction(operator, qtop, (:c, :d)) ≈
          mq["d"] / (mq["c"] + mq["d"])
    @test QuarkModelTransitions._eq19_momentum_fraction(operator, atop, (:c, :d)) ≈
          mq["c"] / (mq["c"] + mq["d"])
    @test_throws ArgumentError QuarkModelTransitions._eq19_momentum_fraction(
        operator, qtop, (:s, :u),
    )

    beta = 0.4
    q = 0.31
    wave = OscillatorWave(0, beta, [1.0])
    direct_q = QuarkModelTransitions._Eq19OrbitalLabel(
        QuarkModelTransitions._DirectPseudoscalarPiece(), qtop, 0, -1,
    )
    direct_qbar = QuarkModelTransitions._Eq19OrbitalLabel(
        QuarkModelTransitions._DirectPseudoscalarPiece(), atop, 0, 1,
    )
    recoil_q = QuarkModelTransitions._Eq19OrbitalLabel(
        QuarkModelTransitions._RecoilPseudoscalarPiece(), qtop, 0, -1,
    )
    recoil_qbar = QuarkModelTransitions._Eq19OrbitalLabel(
        QuarkModelTransitions._RecoilPseudoscalarPiece(), atop, 0, 1,
    )
    gaussian = exp(-q^2 / (16beta^2))
    direct_expected = operator.g * q * gaussian
    recoil_magnitude = operator.h * q * gaussian / 4
    spatial(label, momentum=q) = QuarkModelTransitions._eq19_s_wave_spatial_integral(
        operator, label, wave, wave, momentum, (:u, :d),
    )
    @test spatial(direct_q) ≈ direct_expected rtol = 1e-10
    @test spatial(direct_qbar) ≈ direct_expected rtol = 1e-10
    @test spatial(recoil_q) ≈ recoil_magnitude rtol = 1e-10
    @test spatial(recoil_qbar) ≈ -recoil_magnitude rtol = 1e-10

    # With the Phase-III rho -> pi pi coefficients [-1,-1,-1,+1], the
    # numerical columns reduce to -2*q*F*(g+h/4), i.e. Table IV's A relation.
    combined = -spatial(direct_q) - spatial(direct_qbar) - spatial(recoil_q) +
               spatial(recoil_qbar)
    @test combined ≈ -2q * gaussian * (operator.g + operator.h / 4) rtol = 1e-10

    rho = BasisState(1, "S", 3, 1; label = "rho", flavors = (:u, :u))
    pi = BasisState(1, "S", 1, 0; label = "pi", flavors = (:d, :u))
    decomposition = QuarkModelTransitions._eq19_angular_decomposition(
        pi,
        rho,
        QuarkModelTransitions._appendix_b_flavor(:pi_minus),
        QuarkModelTransitions._appendix_b_flavor(:pi_plus),
        QuarkModelTransitions._appendix_b_flavor(:pi_zero),
    )
    evaluated = QuarkModelTransitions._eq19_wave_amplitude(
        operator, decomposition, wave, wave, q, (:u, :d),
    )
    @test length(evaluated.spatial_integrals) == 4
    @test only(evaluated.helicity).second ≈ combined rtol = 1e-10
    @test only(evaluated.partial_waves).second ≈ combined rtol = 1e-10

    complex_q = 0.17 + 0.08im
    complex_gaussian = exp(-complex_q^2 / (16beta^2))
    @test spatial(direct_q, complex_q) ≈
          operator.g * complex_q * complex_gaussian rtol = 1e-9

    # A sampled native mesh follows the analytic SHO result and exercises the
    # independent finite-difference derivative implementation.
    r = collect(range(0.01, 30.0; step = 0.01))
    mesh = sample_wave(wave, r)
    mesh_spatial(label) = QuarkModelTransitions._eq19_s_wave_spatial_integral(
        operator, label, mesh, mesh, q, (:u, :d),
    )
    @test mesh_spatial(direct_q) ≈ direct_expected rtol = 2e-5
    @test mesh_spatial(recoil_q) ≈ recoil_magnitude rtol = 2e-5

    transverse = QuarkModelTransitions._Eq19OrbitalLabel(
        QuarkModelTransitions._RecoilPseudoscalarPiece(), qtop, 1, -1,
    )
    @test_throws ArgumentError spatial(transverse)
end


@testset "Native solved HO/FD S-wave transition agreement" begin
    params, mq = load_parameters_and_quark_masses(
        joinpath(REPOSITORY_ROOT, "data", "parameters.provisional.toml"),
    )
    masses = ConstituentMasses(mq["q"], mq["q"])
    ho = channel_solution(
        params,
        masses,
        0;
        solver = OscillatorSolver(nlevels_per_channel = 2),
        nlevels = 2,
    )
    fd = channel_solution(
        params,
        masses,
        0;
        solver = FiniteDifferenceSolver(
            ngrid = 1200, rmax = 28.0, nlevels_per_channel = 2,
        ),
        nlevels = 2,
    )
    operator = PseudoscalarEmission(0.7, 0.3, mq)
    topology = QuarkModelTransitions._QuarkEmission()
    labels = (
        QuarkModelTransitions._Eq19OrbitalLabel(
            QuarkModelTransitions._DirectPseudoscalarPiece(), topology, 0, -1,
        ),
        QuarkModelTransitions._Eq19OrbitalLabel(
            QuarkModelTransitions._RecoilPseudoscalarPiece(), topology, 0, -1,
        ),
    )
    value(solution, label, final_level, initial_level) =
        QuarkModelTransitions._eq19_s_wave_spatial_integral(
            operator,
            label,
            radial_wave(solution, final_level),
            radial_wave(solution, initial_level),
            0.45,
            (:q, :q),
        )

    # The 0.3% gate is measured, not guessed: the companion checkpoint records
    # an FD 450--1800 refinement. The worst HO/FD case is the cancellation-prone
    # 1S--2S direct overlap (0.233%); all non-node-sensitive columns are <0.1%.
    for label in labels, (final_level, initial_level) in
        ((1, 1), (1, 2), (2, 1), (2, 2))
        ho_value = value(ho, label, final_level, initial_level)
        fd_value = value(fd, label, final_level, initial_level)
        @test abs(ho_value - fd_value) / max(abs(ho_value), 1e-12) < 3e-3
    end
end
