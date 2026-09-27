@testset "Table V strong-decay model (light 1S+1P)" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)
    # the two fit rows are exact by construction
    @test strong_decay_amplitude(model, sqrt(4 / 3), :A, 1, q_rho) ≈ 12.4 atol = 1e-9
    @test strong_decay_amplitude(model, -sqrt(2 / 9), :S, 0, q_B) ≈ -11.0 atol = 1e-9
    # predictions against the paper's numeric column (2 significant figures)
    q_kstar = decay_momentum(0.8921, 0.4957, 0.138)
    @test isapprox(strong_decay_amplitude(model, 1.0, :A, 1, q_kstar), 7.9; rtol = 0.05)
    q_a2 = decay_momentum(1.318, 0.5488, 0.138)
    @test isapprox(strong_decay_amplitude(model, sqrt(1 / 30), :A, 2, q_a2), 4.5; rtol = 0.05)
    q_f = decay_momentum(1.273, 0.138, 0.138)
    @test isapprox(strong_decay_amplitude(model, -sqrt(1 / 10), :A, 2, q_f), -11.0; rtol = 0.05)
    # below threshold: zero momentum and zero amplitude
    @test decay_momentum(1.0, 0.6, 0.6) == 0.0
    @test strong_decay_amplitude(model, 1.0, :A, 1, 0.0) == 0.0
    # reduced-amplitude classes
    @test reduced_decay_amplitude(model, :A0, 1.0) == model.A
    @test reduced_decay_amplitude(model, :Aprime, 2.0) == model.A
    @test reduced_decay_amplitude(model, :S, 1.0) ≈ model.S0 - 0.5 * model.A
    @test reduced_decay_amplitude(model, :D, 1.0) ≈ model.S0 - 0.3 * model.A
    @test reduced_decay_amplitude(model, :P, 1.0) ≈ model.S0 - 0.75 * model.A
    @test_throws ArgumentError reduced_decay_amplitude(model, :bogus, 1.0)
end

@testset "Table V charmed decays (A_c/S_c, footnote d)" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)

    # A_c is structure-independent, so it is the light A unchanged.
    @test reduced_decay_amplitude(model, :A_c, 3.0) == model.A

    # The charm mass ratio comes from the parameters TOML, not from constants
    # frozen inside src/ (which used to shadow it).
    _, mq_charm = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m_c, m_d = mq_charm["c"], mq_charm["d"]
    r = m_c / (m_c + m_d)
    @test r ≈ 1.628 / (1.628 + 0.220)   # unfreezing changed no value

    # GI (1985) Table IV (PDF page 13) prints
    #   S   = [3h - (1/2)      (g + h/4) q^2/beta^2  ] beta
    #   S_c = [3h - m_c/(m_d+m_c) A_c/beta q^2/beta_c^2] beta_c
    # so the S-family polynomial coefficient is the heavy fraction r, and the
    # light 1/2 is simply r at equal constituent masses. S_c must therefore be
    # MORE suppressed than S (r = 0.881 > 1/2); it is not the light formula.
    @test reduced_decay_amplitude(model, :S_c, 1.0; heavy_fraction = 0.5) ≈
          reduced_decay_amplitude(model, :S, 1.0)
    @test reduced_decay_amplitude(model, :S_c, 1.0; heavy_fraction = r) ≈
          model.S0 - r * model.A
    @test reduced_decay_amplitude(model, :S_c, 1.0; heavy_fraction = r) <
          reduced_decay_amplitude(model, :S, 1.0)

    # The charmed form factor exp[-(1/4)(m_c/(m_c+m_d))^2 q^2/beta^2] replaces
    # the light exp(-q^2/16 beta^2); at fixed q, coefficient, and class the
    # charm amplitude equals the light one rescaled by the Gaussian ratio.
    let q = 0.4, c = 0.9, L = 2, beta = model.beta_GeV
        light = strong_decay_amplitude(model, c, :A, L, q)
        # rebuild the light amplitude under the :A_c class = same reduced value,
        # so the only difference is the form factor and (optionally) recoil.
        charm = strong_decay_amplitude(model, c, :A_c, L, q; heavy_fraction = r)
        ratio = exp(-0.25 * r^2 * q^2 / beta^2) / exp(-q^2 / (16 * beta^2))
        @test isapprox(charm / light, ratio; rtol = 1e-10)
    end

    # The recoil multiplier m_c*beta/((m_c+m_d)*beta_c) = m_c/(m_c+m_d) with
    # beta_c=beta; A_c P-wave rows carry it, S-wave/1^3S_1 rows do not.
    let q = 0.5, c = -sqrt(1 / 5), L = 2
        no_rec = strong_decay_amplitude(model, c, :A_c, L, q; heavy_fraction = r, recoil = false)
        with_rec = strong_decay_amplitude(model, c, :A_c, L, q; heavy_fraction = r, recoil = true)
        # The multiplier itself is exactly r when beta_c = beta, but this ratio
        # is a quotient of two products, so it carries its own rounding and is
        # not bit-exact (it lands 1 ulp off r).
        @test isapprox(with_rec / no_rec, r; rtol = 1e-10)
    end

    # below threshold -> zero
    @test strong_decay_amplitude(model, 1.0, :A_c, 1, 0.0; heavy_fraction = r) == 0.0

    # Clean 1^3S_1 D* -> D pi rows reproduce the paper's column with NO refit
    # (masses: D*+=2.010, D0=1.865, pi+=0.1396; D*0=2.007, pi0=0.135).
    q_dstarp = decay_momentum(2.010, 1.865, 0.1396)
    @test isapprox(strong_decay_amplitude(model, -sqrt(2 / 3), :A_c, 1, q_dstarp; heavy_fraction = r),
        -0.34; atol = 0.05)
    q_dstar0 = decay_momentum(2.007, 1.865, 0.135)
    @test isapprox(strong_decay_amplitude(model, -sqrt(1 / 3), :A_c, 1, q_dstar0; heavy_fraction = r),
        -0.27; atol = 0.05)
end

@testset "Row-oriented decay API (DecayChannel decomposition)" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B) # TableIVPolynomial default
    model_lead = calibrate_strong_decay_model(q_rho, q_B; convention = LeadingS0())
    @test calibrate_strong_decay_model(q_rho, q_B; convention = TableIVPolynomial()) == model
    @test_throws TypeError calibrate_strong_decay_model(q_rho, q_B; convention = "leading")

    # spatial_overlap is the pure [DERIVED] SHO momentum factor.
    @test spatial_overlap(0.0, 1, 0.40) == 0.0
    @test spatial_overlap(0.359, 0, 0.40) ≈
          QuarkModelTransitions._suppressed_factor(0.359, 0.40)

    # The 3-way decomposition multiplies back to the scalar amplitude.
    ch_rho = DecayChannel("rho", "pi", "pi", sqrt(4 / 3), :A, 1)
    a = decay_amplitude(model, ch_rho, q_rho; convention = TableIVPolynomial())
    @test a.total ≈ a.coefficient * a.reduced * a.spatial_overlap
    @test a.total ≈ strong_decay_amplitude(model, sqrt(4 / 3), :A, 1, q_rho)
    @test a.total ≈ 12.4 atol = 1e-9
    @test reduced_matrix_element(a) ≈ a.coefficient * a.reduced
    @test reduced_matrix_element(a) isa Float64
    @test decay_width(a) ≈ a.total^2

    # leading vs table_iv: identical for structure-independent A, differ for D.
    ch_D = DecayChannel("rho2", "omega", "pi", -sqrt(1 / 36), :D, 1)
    q = 0.66
    @test decay_amplitude(model, ch_rho, q; convention = LeadingS0()).total ≈
          decay_amplitude(model, ch_rho, q; convention = TableIVPolynomial()).total
    @test decay_amplitude(model_lead, ch_D, q; convention = LeadingS0()).reduced ≈ model_lead.S0
    @test decay_amplitude(model, ch_D, q; convention = TableIVPolynomial()).reduced <
          decay_amplitude(model, ch_D, q; convention = LeadingS0()).reduced

    # MesonMasses resolves the momentum internally.
    masses = MesonMasses(Dict("rho" => 0.769, "pi" => 0.138))
    @test meson_mass(masses, "rho") == 0.769
    @test masses["pi"] == 0.138
    @test decay_amplitude(model, ch_rho, masses; convention = TableIVPolynomial()).total ≈ 12.4 atol = 1e-9
    @test_throws Exception meson_mass(masses, "unregistered")

    # A charmed channel now carries its OWN mass ratio r = m_c/(m_c+m_d) rather
    # than having it inferred from the class; the A_c P-wave recoil multiplier is
    # still keyed on the class.
    _, mq_row = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r_c = mq_row["c"] / (mq_row["c"] + mq_row["d"])
    ch_c = DecayChannel("Kstar_c", "D", "pi", -sqrt(1 / 5), :A_c, 2; heavy_fraction = r_c)
    q_c = 0.4
    got = decay_amplitude(model, ch_c, q_c; convention = LeadingS0()).total
    want = strong_decay_amplitude(model, -sqrt(1 / 5), :A_c, 2, q_c;
        heavy_fraction = r_c, recoil = true, convention = LeadingS0())
    @test got ≈ want

    # Omitting it on an unequal-mass class is a construction error, never a
    # silent fallback to the light (r = 1/2) form factor.
    @test_throws ArgumentError DecayChannel("Kstar_c", "D", "pi", 1.0, :A_c, 2)
    @test_throws ArgumentError DecayChannel("x", "y", "z", 1.0, :A, 1; heavy_fraction = 1.5)
    @test DecayChannel("rho", "pi", "pi", 1.0, :A, 1).heavy_fraction == 0.5
end
