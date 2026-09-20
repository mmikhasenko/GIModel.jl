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

@testset "Coherent three-state composition" begin
    wave = OscillatorWave(0, 0.4, [1.0])
    parent_basis = [
        BasisState(1, "P", 3, 0; label = "A1", flavors = (:u, :u)),
        BasisState(2, "P", 3, 0; label = "A2", flavors = (:u, :u)),
    ]
    first_basis = [
        BasisState(1, "S", 1, 0; label = "B1", flavors = (:u, :d)),
        BasisState(2, "S", 1, 0; label = "B2", flavors = (:u, :d)),
    ]
    second_basis = [
        BasisState(1, "S", 1, 0; label = "C1", flavors = (:d, :u)),
        BasisState(2, "S", 1, 0; label = "C2", flavors = (:d, :u)),
    ]
    parent_coefficients = (sqrt(0.6), im * sqrt(0.4))
    first_coefficients = (inv(sqrt(2)), im / sqrt(2))
    second_coefficients = (sqrt(0.7), -im * sqrt(0.3))
    parent = PhysicalState("A", 1.8, [
        (basis = parent_basis[i], coefficient = parent_coefficients[i], wave = wave)
        for i in eachindex(parent_basis)
    ])
    first = PhysicalState("B", 0.4, [
        (basis = first_basis[i], coefficient = first_coefficients[i], wave = wave)
        for i in eachindex(first_basis)
    ])
    second = PhysicalState("C", 0.5, [
        (basis = second_basis[i], coefficient = second_coefficients[i], wave = wave)
        for i in eachindex(second_basis)
    ])
    final = TwoMesonChannel(first, second)
    wave_S = PartialWave(0, 0)
    kernel(b, c, a) = complex(
        10 * a.basis.n + 2 * b.basis.n + c.basis.n,
        b.basis.n - c.basis.n,
    )

    composed = GIModel._compose_physical_decay(kernel, final, parent, wave_S)
    expected = sum(
        a.coefficient * conj(b.coefficient) * conj(c.coefficient) * kernel(b, c, a)
        for a in parent.components for b in final.first.components for c in final.second.components
    )
    @test composed.value ≈ expected
    @test length(composed.terms) == 8
    @test composed.external_masses_GeV ==
          (parent.mass_GeV, final.first.mass_GeV, final.second.mass_GeV)
    @test all(type -> type !== Any, fieldtypes(typeof(composed)))
    @test all(type -> type !== Any, fieldtypes(typeof(composed.terms[1])))
    @test all(term.identical_normalization == 1.0 for term in composed.terms)
    @test all(
        term.mixing_coefficient == term.initial_component.coefficient *
        conj(term.first_component.coefficient) * conj(term.second_component.coefficient)
        for term in composed.terms
    )

    phases = (cis(0.31), cis(-0.47), cis(0.83))
    rephased_parent = PhysicalState("A", 1.8, [
        (basis = component.basis, coefficient = phases[1] * component.coefficient,
         wave = component.wave) for component in parent.components
    ])
    rephased_first = PhysicalState("B", 0.4, [
        (basis = component.basis, coefficient = phases[2] * component.coefficient,
         wave = component.wave) for component in first.components
    ])
    rephased_second = PhysicalState("C", 0.5, [
        (basis = component.basis, coefficient = phases[3] * component.coefficient,
         wave = component.wave) for component in second.components
    ])
    rephased = GIModel._compose_physical_decay(
        kernel,
        TwoMesonChannel(rephased_first, rephased_second),
        rephased_parent,
        wave_S,
    )
    covariance = phases[1] * conj(phases[2]) * conj(phases[3])
    @test rephased.value ≈ covariance * composed.value
    @test abs2(rephased.value) ≈ abs2(composed.value)

    pure_first = PhysicalState("A", 1.8,
        [(basis = parent_basis[1], coefficient = 1.0, wave = wave)])
    mixed_plus = PhysicalState("A", 1.8, [
        (basis = parent_basis[1], coefficient = inv(sqrt(2)), wave = wave),
        (basis = parent_basis[2], coefficient = inv(sqrt(2)), wave = wave),
    ])
    mixed_minus = PhysicalState("A", 1.8, [
        (basis = parent_basis[1], coefficient = inv(sqrt(2)), wave = wave),
        (basis = parent_basis[2], coefficient = -inv(sqrt(2)), wave = wave),
    ])
    pure_first_daughter = PhysicalState("B", 0.4,
        [(basis = first_basis[1], coefficient = 1.0, wave = wave)])
    pure_second_daughter = PhysicalState("C", 0.5,
        [(basis = second_basis[1], coefficient = 1.0, wave = wave)])
    pure_final = TwoMesonChannel(pure_first_daughter, pure_second_daughter)
    equal_kernel(b, c, a) = 1.0
    mixing_only_kernel(b, c, a) = a.basis.n == 2 ? 2.0 : 0.0
    @test GIModel._compose_physical_decay(equal_kernel, pure_final, mixed_plus, wave_S).value ≈ sqrt(2)
    @test GIModel._compose_physical_decay(equal_kernel, pure_final, mixed_minus, wave_S).value ≈ 0.0
    @test GIModel._compose_physical_decay(mixing_only_kernel, pure_final, pure_first, wave_S).value == 0.0
    @test GIModel._compose_physical_decay(mixing_only_kernel, pure_final, mixed_plus, wave_S).value ≈ sqrt(2)

    identical = PhysicalState("B", 0.4, [
        (basis = first_basis[1], coefficient = inv(sqrt(2)), wave = wave),
        (basis = first_basis[2], coefficient = im / sqrt(2), wave = wave),
    ])
    identical_reordered_rephased = PhysicalState("B", 0.4, [
        (basis = first_basis[2], coefficient = -inv(sqrt(2)), wave = wave),
        (basis = first_basis[1], coefficient = im / sqrt(2), wave = wave),
    ])
    identical_final = TwoMesonChannel(identical, identical_reordered_rephased)
    @test GIModel._same_external_state(identical_final.first, identical_final.second)
    relative_rephased = PhysicalState("B", 0.4, [
        (basis = first_basis[1], coefficient = inv(sqrt(2)), wave = wave),
        (basis = first_basis[2], coefficient = -im / sqrt(2), wave = wave),
    ])
    @test !GIModel._same_external_state(identical, relative_rephased)
    @test allowed_partial_waves(identical_final, pure_first) == [wave_S]
    identical_composed = GIModel._compose_physical_decay(
        equal_kernel, identical_final, pure_first, wave_S,
    )
    raw_identical = sum(
        a.coefficient * conj(b.coefficient) * conj(c.coefficient)
        for a in pure_first.components for b in identical_final.first.components for
        c in identical_final.second.components
    )
    @test identical_composed.symmetry.identical
    @test identical_composed.symmetry.exchange_phase == 1
    @test identical_composed.value ≈ raw_identical / sqrt(2)
    @test all(
        term.identical_normalization == inv(sqrt(2)) for term in identical_composed.terms
    )

    vector_parent_basis = BasisState(1, "S", 3, 1; label = "V", flavors = (:u, :u))
    vector_parent = PhysicalState("V", 1.8,
        [(basis = vector_parent_basis, coefficient = 1.0, wave = wave)])
    @test isempty(allowed_partial_waves(identical_final, vector_parent))
    @test_throws ArgumentError GIModel._compose_physical_decay(
        equal_kernel, identical_final, vector_parent, PartialWave(1, 0),
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
