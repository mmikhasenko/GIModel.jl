struct DummyNativeDecayOperator <: StrongDecayOperator end

@testset "Transition-amplitude domain objects" begin
    wave = OscillatorWave(0, 0.4, [1.0])
    rho_basis = BasisState(1, "S", 3, 1; label = "rho", flavors = (:u, :d))
    pip_basis = BasisState(1, "S", 1, 0; label = "pi+", flavors = (:u, :d))
    pim_basis = BasisState(1, "S", 1, 0; label = "pi-", flavors = (:d, :u))

    rho = PhysicalState("rho", 0.769,
        [(basis = rho_basis, coefficient = 1.0, wave = wave)])
    pip = PhysicalState("pi+", 0.138,
        [(basis = pip_basis, coefficient = 1.0, wave = wave)])
    pim = PhysicalState("pi-", 0.138,
        [(basis = pim_basis, coefficient = 1.0, wave = wave)])

    @test rho.J == 1
    @test rho.parity == -1
    @test length(rho.components) == 1
    @test_throws ArgumentError PhysicalState("empty", 1.0, [])
    @test_throws ArgumentError PhysicalState("bad mass", -1.0,
        [(basis = rho_basis, coefficient = 1.0, wave = wave)])

    final = TwoMesonChannel(pip, pim)
    reversed = TwoMesonChannel(pim, pip)
    @test final.first.label == reversed.first.label
    @test final.second.label == reversed.second.label
    @test allowed_partial_waves(final, rho) == [PartialWave(1, 0)]

    projection = partial_wave_projection(final, rho)
    @test projection.partial_waves == (PartialWave(1, 0),)
    @test projection.matrix == ones(1, 1)
    vector = ReferenceState("vector", 0.8; J = 1, parity = -1)
    @test_throws ArgumentError partial_wave_projection(
        TwoMesonChannel(vector, ReferenceState("p", 0.1; J = 0, parity = -1)),
        ReferenceState("parent", 1.2; J = 1, parity = 1),
    )

    identical = TwoMesonChannel(pip, pip)
    @test isempty(allowed_partial_waves(identical, rho))
    @test_throws ArgumentError PartialWave(-1, 0)
    @test_throws ArgumentError CMKinematics(-0.1)
    @test CMKinematics(0.2im).momentum_GeV == 0.2im

    incomplete = ReferenceState("epsilon", 0.7)
    @test_throws ArgumentError allowed_partial_waves(
        TwoMesonChannel(incomplete, ReferenceState("pi", 0.138)),
        ReferenceState("A1", 1.2),
    )
end

@testset "Typed Table V reference adapter" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)
    row = DecayChannel("rho", "pi+", "pi-", sqrt(4 / 3), :A, 1;
        label = "rho -> pi pi")
    operator = TableVReference(
        model,
        row,
        PartialWave(1, 0);
        convention = TableIVPolynomial(),
    )
    parent = ReferenceState("rho", 0.769; J = 1, parity = -1)
    final = TwoMesonChannel(
        ReferenceState("pi+", 0.138; J = 0, parity = -1),
        ReferenceState("pi-", 0.138; J = 0, parity = -1),
    )

    amplitude = matrix_element(final, operator, parent; kinematics = OnShell())
    @test @inferred(matrix_element(final, operator, parent; kinematics = OnShell())) isa
          TransitionAmplitude
    legacy = decay_amplitude(model, row, q_rho; convention = TableIVPolynomial())
    @test isempty(amplitude.helicity)
    @test partial_waves(amplitude) == [PartialWave(1, 0)]
    @test amplitude[PartialWave(1, 0)] == legacy.total
    @test decay_width(amplitude) == decay_width(legacy)
    @test amplitude.kinematics.momentum_GeV ≈ q_rho
    @test amplitude.provenance.helicity_available == false
    @test all(type -> type !== Any, fieldtypes(typeof(amplitude)))
    @test occursin("rho -> pi+ + pi-", sprint(show, MIME"text/plain"(), amplitude))
    @test_throws KeyError amplitude[PartialWave(3, 0)]
    @test_throws ArgumentError matrix_element(
        final, operator, parent; kinematics = CMKinematics(0.1im),
    )
    @test_throws ArgumentError matrix_element(
        final, DummyNativeDecayOperator(), parent; kinematics = CMKinematics(0.1),
    )

    closed_parent = ReferenceState("x", 0.2; J = 0, parity = 1)
    closed_final = TwoMesonChannel(
        ReferenceState("a", 0.15; J = 0, parity = -1),
        ReferenceState("b", 0.15; J = 0, parity = -1),
    )
    closed_row = DecayChannel("x", "a", "b", 1.0, :A, 0)
    closed_operator = TableVReference(model, closed_row, PartialWave(0, 0))
    @test decay_width(closed_final, closed_operator, closed_parent) == 0.0
    @test_throws GIModel.ClosedChannelError matrix_element(
        closed_final, closed_operator, closed_parent; kinematics = OnShell(),
    )

    quasi_row = DecayChannel("A1", "epsilon", "pi", 1.0, :A0, 0)
    quasi_operator = TableVReference(model, quasi_row, PartialWave(0, 0))
    quasi_parent = ReferenceState("A1", 1.4)
    quasi_final = TwoMesonChannel(
        ReferenceState("epsilon", 0.7),
        ReferenceState("pi", 0.138),
    )
    quasi = matrix_element(
        quasi_final,
        quasi_operator,
        quasi_parent;
        kinematics = CMKinematics(0.2),
    )
    @test quasi[PartialWave(0, 0)] ==
          decay_amplitude(model, quasi_row, 0.2).total
end
