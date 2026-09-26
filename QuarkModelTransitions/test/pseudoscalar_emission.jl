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
    @test spatial(transverse) == 0
end

@testset "General Eq. (19) orbital integrals" begin
    mq = QuarkMassTable("u" => 0.22, "d" => 0.22)
    operator = PseudoscalarEmission(0.73, -0.28, mq)
    qtop = QuarkModelTransitions._QuarkEmission()
    atop = QuarkModelTransitions._AntiquarkEmission()
    direct = QuarkModelTransitions._DirectPseudoscalarPiece()
    recoil = QuarkModelTransitions._RecoilPseudoscalarPiece()
    beta = 0.4
    momentum = 0.31
    alpha = 0.5
    s_wave = OscillatorWave(0, beta, [1.0])
    p_wave = OscillatorWave(1, beta, [1.0])
    d_wave = OscillatorWave(2, beta, [1.0])

    orbital(piece, topology, component, sign, Li, mi, Lf, mf) =
        QuarkModelTransitions._SpectroscopicOrbitalLabel(
            QuarkModelTransitions._Eq19OrbitalLabel(
                piece, topology, component, sign,
            ),
            Li,
            mi,
            Lf,
            mf,
        )
    spatial(label, final_wave, initial_wave, q=momentum) =
        QuarkModelTransitions._eq19_spatial_integral(
            operator, label, final_wave, initial_wave, q, (:u, :d),
        )

    # Generic multipole machinery must reduce exactly to the old S-to-S slice.
    for piece in (direct, recoil), topology in (qtop, atop)
        sign = topology isa QuarkModelTransitions._QuarkEmission ? -1 : 1
        elementary = QuarkModelTransitions._Eq19OrbitalLabel(
            piece, topology, 0, sign,
        )
        general = orbital(piece, topology, 0, sign, 0, 0, 0, 0)
        @test spatial(general, s_wave, s_wave) ≈
              QuarkModelTransitions._eq19_s_wave_spatial_integral(
                  operator, elementary, s_wave, s_wave, momentum, (:u, :d),
              ) rtol = 2e-13
    end

    # Independent Cartesian Gaussian result:
    # <P_0|exp(s*i*k*z)|S> = s*i*k/(sqrt(2)*beta) exp[-k^2/(4beta^2)].
    gaussian = exp(-(alpha * momentum)^2 / (4beta^2))
    for (topology, sign) in ((qtop, -1), (atop, 1))
        s_to_p = orbital(direct, topology, 0, sign, 0, 0, 1, 0)
        expected = operator.g * momentum * sign * im * alpha * momentum /
                   (sqrt(2) * beta) * gaussian
        @test spatial(s_to_p, p_wave, s_wave) ≈ expected rtol = 2e-10

        # The axial plane wave preserves orbital m; forbidden projections are
        # exact zeros rather than targeted special cases.
        forbidden = orbital(direct, topology, 0, sign, 0, 0, 1, 1)
        @test spatial(forbidden, p_wave, s_wave) == 0
    end

    # At q=0, the two gradient branches reproduce the analytic oscillator
    # momentum matrix element and Hermitian conjugation between S and P.
    p_to_s = orbital(recoil, qtop, 0, -1, 1, 0, 0, 0)
    s_to_p = orbital(recoil, qtop, 0, -1, 0, 0, 1, 0)
    @test spatial(p_to_s, s_wave, p_wave, 0.0) ≈
          im * operator.h * beta / sqrt(2) atol = 2e-10
    @test spatial(s_to_p, p_wave, s_wave, 0.0) ≈
          -im * operator.h * beta / sqrt(2) atol = 2e-10
    @test spatial(p_to_s, s_wave, p_wave, 0.0) ≈
          conj(spatial(s_to_p, p_wave, s_wave, 0.0)) atol = 2e-10

    # Rotational covariance: all three spherical components of S<->P have the
    # same magnitude when their magnetic projections satisfy m_i=m_f-nu.
    transverse_values = ComplexF64[]
    for component in -1:1
        label = orbital(recoil, qtop, component, -1, 1, -component, 0, 0)
        push!(transverse_values, spatial(label, s_wave, p_wave, 0.0))
    end
    @test all(
        value -> isapprox(
            abs(value), abs(first(transverse_values)); atol = 2e-10,
        ),
        transverse_values,
    )

    # Nothing in the evaluator is specialized to L<=1 or to real momentum.
    generic_d = orbital(recoil, atop, -1, 1, 2, 1, 1, 0)
    @test isfinite(spatial(generic_d, p_wave, d_wave, 0.19 + 0.07im))
    @test_throws ArgumentError spatial(generic_d, s_wave, d_wave)

    # One P-wave parent produces two helicities and correlated S/D partial
    # waves from the same general integral columns; no decay-class formula is
    # selected by the evaluator.
    rho = BasisState(1, "S", 3, 1; label = "rho", flavors = (:u, :u))
    a1 = BasisState(1, "P", 3, 1; label = "A1", flavors = (:u, :u))
    decomposition = QuarkModelTransitions._eq19_angular_decomposition(
        rho,
        a1,
        QuarkModelTransitions._appendix_b_flavor(:pi_zero),
        QuarkModelTransitions._appendix_b_flavor(:pi_plus),
        QuarkModelTransitions._appendix_b_flavor(:pi_plus),
    )
    amplitude = QuarkModelTransitions._eq19_wave_amplitude(
        operator, decomposition, s_wave, p_wave, momentum, (:u, :d),
    )
    @test length(amplitude.helicity) == 2
    @test length(amplitude.partial_waves) == 2
    @test all(isfinite(last(value)) for value in amplitude.helicity)
    @test all(isfinite(last(value)) for value in amplitude.partial_waves)
    @test Set(first(value).relative_L for value in amplitude.partial_waves) == Set((0, 2))
    table_iv_A = (operator.g + operator.h / 4) * beta
    table_iv_S = (
        3operator.h -
        0.5 * (operator.g + operator.h / 4) * momentum^2 / beta^2
    ) * beta
    expected_D_over_S = -table_iv_A * (momentum / beta)^2 /
                        (sqrt(8) * table_iv_S)
    @test last(amplitude.partial_waves[2]) / last(amplitude.partial_waves[1]) ≈
          expected_D_over_S rtol = 2e-12

    for order in 0:6
        x = 0.37 + 0.11im
        value = QuarkModelTransitions._spherical_bessel_j(order, x)
        @test isfinite(value)
        if order >= 2
            @test value ≈
                  (2order - 1) / x *
                  QuarkModelTransitions._spherical_bessel_j(order - 1, x) -
                  QuarkModelTransitions._spherical_bessel_j(order - 2, x) rtol = 1e-10
        end
    end
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
    ho_p = channel_solution(
        params,
        masses,
        1;
        solver = OscillatorSolver(nlevels_per_channel = 1),
        nlevels = 1,
    )
    fd_p = channel_solution(
        params,
        masses,
        1;
        solver = FiniteDifferenceSolver(
            ngrid = 1200, rmax = 28.0, nlevels_per_channel = 1,
        ),
        nlevels = 1,
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

    # The first genuinely general-orbital certificate: P->S direct and all
    # three recoil components agree across independent native solvers below
    # 0.1%. Magnetic labels follow m_i=m_f-nu.
    general_labels = (
        QuarkModelTransitions._SpectroscopicOrbitalLabel(
            labels[1], 1, 0, 0, 0,
        ),
        QuarkModelTransitions._SpectroscopicOrbitalLabel(
            labels[2], 1, 0, 0, 0,
        ),
        QuarkModelTransitions._SpectroscopicOrbitalLabel(
            QuarkModelTransitions._Eq19OrbitalLabel(
                QuarkModelTransitions._RecoilPseudoscalarPiece(), topology, 1, -1,
            ),
            1, -1, 0, 0,
        ),
        QuarkModelTransitions._SpectroscopicOrbitalLabel(
            QuarkModelTransitions._Eq19OrbitalLabel(
                QuarkModelTransitions._RecoilPseudoscalarPiece(), topology, -1, -1,
            ),
            1, 1, 0, 0,
        ),
    )
    for label in general_labels
        ho_value = QuarkModelTransitions._eq19_spatial_integral(
            operator,
            label,
            radial_wave(ho, 1),
            radial_wave(ho_p, 1),
            0.45,
            (:q, :q),
        )
        fd_value = QuarkModelTransitions._eq19_spatial_integral(
            operator,
            label,
            radial_wave(fd, 1),
            radial_wave(fd_p, 1),
            0.45,
            (:q, :q),
        )
        @test abs(ho_value - fd_value) / max(abs(ho_value), 1e-12) < 1e-3
    end
end
