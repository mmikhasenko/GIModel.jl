@testset "GI Eq. (19) pure coefficient decomposition" begin
    direct = QuarkModelTransitions._DirectPseudoscalarPiece()
    recoil = QuarkModelTransitions._RecoilPseudoscalarPiece()
    q = QuarkModelTransitions._QuarkEmission
    qbar = QuarkModelTransitions._AntiquarkEmission

    singlet = QuarkModelTransitions._coupled_spin(0, 0)
    triplet0 = QuarkModelTransitions._coupled_spin(1, 0)
    rho0 = QuarkModelTransitions._appendix_b_flavor(:pi_zero)
    piminus = QuarkModelTransitions._appendix_b_flavor(:pi_minus)
    piplus = QuarkModelTransitions._appendix_b_flavor(:pi_plus)

    direct_terms = QuarkModelTransitions._eq19_coefficient_decomposition(
        direct, singlet, triplet0, piminus, piplus, rho0,
    )
    recoil_terms = QuarkModelTransitions._eq19_coefficient_decomposition(
        recoil, singlet, triplet0, piminus, piplus, rho0,
    )

    @test length(direct_terms) == 6
    @test all(term -> term.provenance.source == :GI1985_Eq19, direct_terms)
    @test all(term -> term.provenance.coupling == :g, direct_terms)
    @test all(term -> term.provenance.coupling == :h, recoil_terms)
    @test all(term -> term.orbital_label.vector_component == -term.spin_component,
              direct_terms)

    direct_q0 = only(term for term in direct_terms
                     if term.topology isa q && term.spin_component == 0)
    direct_qbar0 = only(term for term in direct_terms
                        if term.topology isa qbar && term.spin_component == 0)
    @test direct_q0.flavor_factor == -1
    @test direct_qbar0.flavor_factor == -1
    @test direct_q0.spin_factor ≈ 1
    @test direct_qbar0.spin_factor ≈ -1
    @test direct_q0.topology_phase == 1
    @test direct_qbar0.topology_phase == -1
    @test direct_q0.coefficient ≈ -1
    @test direct_qbar0.coefficient ≈ -1
    @test direct_q0.orbital_label.plane_wave_sign == -1
    @test direct_qbar0.orbital_label.plane_wave_sign == 1

    recoil_q0 = only(term for term in recoil_terms
                     if term.topology isa q && term.spin_component == 0)
    recoil_qbar0 = only(term for term in recoil_terms
                        if term.topology isa qbar && term.spin_component == 0)
    @test recoil_q0.topology_phase == 1
    @test recoil_qbar0.topology_phase == 1
    @test recoil_q0.coefficient ≈ -1
    @test recoil_qbar0.coefficient ≈ 1

    # In the helicity frame q points along +z, so g sigma.q has only mu=0.
    @test all(iszero(term.coefficient) for term in direct_terms
              if term.spin_component != 0)
    @test all(term.component_selection == (term.spin_component == 0 ? 1 : 0)
              for term in direct_terms)
    @test all(term.component_selection == 1 for term in recoil_terms)

    # The spherical scalar product is sum_mu (-1)^mu sigma_mu v_-mu.
    triplet_minus = QuarkModelTransitions._coupled_spin(1, -1)
    transverse = QuarkModelTransitions._eq19_coefficient_decomposition(
        recoil, singlet, triplet_minus, piminus, piplus, rho0,
    )
    qplus = only(term for term in transverse
                 if term.topology isa q && term.spin_component == 1)
    @test qplus.spherical_phase == -1
    @test qplus.component_selection == 1
    @test qplus.orbital_label.vector_component == -1
    @test qplus.coefficient ≈ -1

    # A spectator/OZI-forbidden channel remains exactly zero at the term level.
    pure_ss = QuarkModelTransitions._appendix_b_flavor(:M_s)
    forbidden = QuarkModelTransitions._eq19_coefficient_decomposition(
        direct, singlet, triplet0, piminus, piplus, pure_ss,
    )
    @test all(term -> iszero(term.flavor_factor), forbidden)
    @test all(term -> iszero(term.coefficient), forbidden)

    @test_throws ArgumentError QuarkModelTransitions.orbital_decomposition(
        direct, q(), 2,
    )
end

@testset "Spectroscopic helicity and shared-integral decomposition" begin
    rho = BasisState(1, "S", 3, 1; label = "rho", flavors = (:u, :u))
    pi = BasisState(1, "S", 1, 0; label = "pi", flavors = (:d, :u))
    a1 = BasisState(1, "P", 3, 1; label = "A1", flavors = (:u, :u))
    flavors = (
        daughter = QuarkModelTransitions._appendix_b_flavor(:pi_minus),
        emitted = QuarkModelTransitions._appendix_b_flavor(:pi_plus),
        parent = QuarkModelTransitions._appendix_b_flavor(:pi_zero),
    )

    rho_to_pipi = QuarkModelTransitions._eq19_angular_decomposition(
        pi, rho, flavors.daughter, flavors.emitted, flavors.parent,
    )
    @test length(rho_to_pipi.helicities) == 1
    @test rho_to_pipi.partial_waves == (PartialWave(1, 0),)
    @test size(rho_to_pipi.helicity_coefficients, 1) == 1
    @test rho_to_pipi.partial_wave_coefficients ==
          rho_to_pipi.helicity_coefficients
    @test rho_to_pipi.provenance.numerical_spatial_integrals == false
    # Independent hand reduction of B1--B3, B25--B29, and Eq. (19), ordered
    # as direct-q, direct-qbar, recoil-q, recoil-qbar.
    @test rho_to_pipi.helicity_coefficients ≈
          ComplexF64[-1 -1 -1 1]

    # A1 -> rho pi has h0/h1 and S/D waves. Both waves are rows over the exact
    # same spatial-integral columns; only the Appendix-C linear combination differs.
    a1_to_rhopi = QuarkModelTransitions._eq19_angular_decomposition(
        rho,
        a1,
        QuarkModelTransitions._appendix_b_flavor(:pi_zero),
        QuarkModelTransitions._appendix_b_flavor(:pi_plus),
        QuarkModelTransitions._appendix_b_flavor(:pi_plus),
    )
    @test length(a1_to_rhopi.helicities) == 2
    @test a1_to_rhopi.partial_waves == (PartialWave(0, 1), PartialWave(2, 1))
    @test size(a1_to_rhopi.helicity_coefficients) ==
          (2, length(a1_to_rhopi.integral_labels))
    @test size(a1_to_rhopi.partial_wave_coefficients) ==
          (2, length(a1_to_rhopi.integral_labels))
    @test !isempty(a1_to_rhopi.integral_labels)
    expected_projection = [sqrt(1 / 3) sqrt(2 / 3);
                           -sqrt(2 / 3) sqrt(1 / 3)]
    @test a1_to_rhopi.partial_wave_coefficients ≈
          expected_projection * a1_to_rhopi.helicity_coefficients
    @test a1_to_rhopi.helicity_coefficients ≈ ComplexF64[
        -inv(sqrt(2)) -inv(sqrt(2)) -inv(sqrt(2)) -inv(sqrt(2)) 0 0 0 0;
        0 0 -1 -1 1 -1 1 1
    ]

    # Exact flavor zeros survive the entire recoupling/projection pipeline.
    forbidden = QuarkModelTransitions._eq19_angular_decomposition(
        pi,
        rho,
        flavors.daughter,
        flavors.emitted,
        QuarkModelTransitions._appendix_b_flavor(:M_s),
    )
    @test isempty(forbidden.integral_labels)
    @test size(forbidden.helicity_coefficients) == (1, 0)
    @test size(forbidden.partial_wave_coefficients) == (1, 0)

    # The same machinery accepts unequal-mass/heavy flavor without a special
    # angular route; mass-dependent q' belongs to the Phase-IV integral value.
    heavy_parent = QuarkModelTransitions._heavy_flavor_state(:c, :d)
    heavy_daughter = QuarkModelTransitions._heavy_flavor_state(:c, :u)
    heavy = QuarkModelTransitions._eq19_angular_decomposition(
        pi,
        rho,
        heavy_daughter,
        QuarkModelTransitions._appendix_b_flavor(:pi_plus),
        heavy_parent,
    )
    @test !isempty(heavy.integral_labels)

    # A pure ss tensor cannot emit pi+ into pi-, but its nonstrange physical
    # admixture produces a nonzero symbolic helicity vector coherently.
    tensor = BasisState(1, "P", 3, 2; label = "fprime", flavors = (:s, :s))
    pure_ss = QuarkModelTransitions._appendix_b_flavor(:M_s)
    angle = 0.12
    mixed_fprime = QuarkModelTransitions._FlavorState("fprime", [
        (:s, :s) => cos(angle),
        (:u, :u) => sin(angle) / sqrt(2),
        (:d, :d) => sin(angle) / sqrt(2),
    ])
    pure_tensor = QuarkModelTransitions._eq19_angular_decomposition(
        pi, tensor, flavors.daughter, flavors.emitted, pure_ss,
    )
    mixed_tensor = QuarkModelTransitions._eq19_angular_decomposition(
        pi, tensor, flavors.daughter, flavors.emitted, mixed_fprime,
    )
    @test isempty(pure_tensor.integral_labels)
    @test !isempty(mixed_tensor.integral_labels)
    @test any(!iszero, mixed_tensor.partial_wave_coefficients)

    @test_throws ArgumentError QuarkModelTransitions._spectroscopic_decomposition(
        QuarkModelTransitions._DirectPseudoscalarPiece(),
        pi,
        rho,
        flavors.daughter,
        flavors.emitted,
        flavors.parent,
        2,
    )
end
