@testset "Isoscalar electromagnetic currents resolve coherent light flavors" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    levels = spectrum_levels(1; L_labels = ("S",))
    solver = FiniteDifferenceSolver(ngrid = 90, nlevels_per_channel = 1)
    iso = compute_isoscalar_spectrum(params, Meson(mq, :q, :q), Meson(mq, :s, :s);
        levels, solver, amplitudes = Dict(("S", 3, 1) => 2.5))
    phi = physical_state(iso, BasisState(1, "S", 3, 1; flavors = (:s, :s)))
    omega = physical_state(iso, BasisState(1, "S", 3, 1; flavors = (:q, :q)))
    eta = physical_state(iso, BasisState(1, "S", 1, 0; flavors = (:q, :q)))
    em = LeptonicCurrent(:electromagnetic, mq)
    effective = LeptonicCurrent(:V_V, mq, [
        AnnihilationTerm((:q, :q), 2sqrt(3) / (3sqrt(2))),
        AnnihilationTerm((:s, :s), -2sqrt(3) / 3),
    ])
    for state in (omega, phi)
        amplitude = matrix_element(Vacuum(), em, state; npoints = 120)
        @test amplitude.value ≈ matrix_element(Vacuum(), effective, state; npoints = 120).value
        @test Set(t.basis.flavors for t in amplitude.terms) == Set(((:u,:u), (:d,:d), (:s,:s)))
        @test decay_width(MasslessLeptonPair(), em, state; npoints = 120) > 0
        @test isfinite(mass_correction_factor(Vacuum(), em, state;
                                             target_mass = 1.1state.mass_GeV, npoints = 120))
    end
    @test decay_width(eta, PhotonEmission(mq), phi) > 0
    # The coarse source is unchanged for explicit effective-charge operators.
    @test any(c.basis.flavors == (:q, :q) for c in phi.components)
    ambiguous = physical_state(compute_spectrum(params, Meson(mq, :q, :q);
        levels, solver), "1^3S_1")
    @test_throws ArgumentError matrix_element(Vacuum(), em, ambiguous; npoints = 120)
end
