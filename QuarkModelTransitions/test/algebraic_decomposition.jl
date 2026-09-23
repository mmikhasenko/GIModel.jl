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
