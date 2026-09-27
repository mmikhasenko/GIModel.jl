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
    @test rho_to_pipi.partial_waves == (PartialWave(1, 0),)
    @test rho_to_pipi.partial_wave_coefficients ==
          rho_to_pipi.helicity_coefficients
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
    @test a1_to_rhopi.partial_waves == (PartialWave(0, 1), PartialWave(2, 1))
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
    @test any(!iszero, mixed_tensor.partial_wave_coefficients)

end
