@testset "GI Appendix C / Table XI projection" begin
    vector = ReferenceState("V", 0.7; J = 1, parity = -1)
    pseudoscalar = ReferenceState("P", 0.1; J = 0, parity = -1)
    final = TwoMesonChannel(vector, pseudoscalar)

    expected = Dict(
        (0, -1) => [(1, [-1.0, 0.0])],
        (1, 1) => [
            (0, [sqrt(1 / 3), sqrt(2 / 3)]),
            (2, [-sqrt(2 / 3), sqrt(1 / 3)]),
        ],
        (1, -1) => [(1, [0.0, -1.0])],
        (2, 1) => [(2, [0.0, -1.0])],
        (2, -1) => [
            (1, [sqrt(2 / 5), sqrt(3 / 5)]),
            (3, [-sqrt(3 / 5), sqrt(2 / 5)]),
        ],
        (3, 1) => [
            (2, [sqrt(3 / 7), sqrt(4 / 7)]),
            (4, [-sqrt(4 / 7), sqrt(3 / 7)]),
        ],
        (3, -1) => [(3, [0.0, -1.0])],
        (4, 1) => [(4, [0.0, -1.0])],
        (4, -1) => [
            (3, [sqrt(4 / 9), sqrt(5 / 9)]),
            (5, [-sqrt(5 / 9), sqrt(4 / 9)]),
        ],
        (5, 1) => [
            (4, [sqrt(5 / 11), sqrt(6 / 11)]),
            (6, [-sqrt(6 / 11), sqrt(5 / 11)]),
        ],
        (5, -1) => [(5, [0.0, -1.0])],
    )

    empty_projection = partial_wave_projection(
        final, ReferenceState("J0+", 2.0; J = 0, parity = 1),
    )
    @test isempty(empty_projection.partial_waves)
    @test size(empty_projection.matrix) == (0, 2)

    for ((J, parity), rows) in expected
        projection = partial_wave_projection(
            final,
            ReferenceState("parent", 2.0; J, parity),
        )
        @test [wave.relative_L for wave in projection.partial_waves] == first.(rows)
        @test projection.matrix ≈ reduce(vcat, permutedims(values) for (_, values) in rows)
        identity = [i == j ? 1.0 : 0.0 for i in eachindex(rows), j in eachindex(rows)]
        @test projection.matrix * transpose(projection.matrix) ≈
              identity
    end
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

    composed = QuarkModelTransitions._compose_physical_decay(kernel, final, parent, wave_S)
    expected = sum(
        a.coefficient * conj(b.coefficient) * conj(c.coefficient) * kernel(b, c, a)
        for a in parent.components for b in final.first.components for c in final.second.components
    )
    @test composed.value ≈ expected
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
    rephased = QuarkModelTransitions._compose_physical_decay(
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
    @test QuarkModelTransitions._compose_physical_decay(
        equal_kernel, pure_final, mixed_plus, wave_S,
    ).value ≈ sqrt(2)
    @test QuarkModelTransitions._compose_physical_decay(
        equal_kernel, pure_final, mixed_minus, wave_S,
    ).value ≈ 0.0
    @test QuarkModelTransitions._compose_physical_decay(
        mixing_only_kernel, pure_final, pure_first, wave_S,
    ).value == 0.0
    @test QuarkModelTransitions._compose_physical_decay(
        mixing_only_kernel, pure_final, mixed_plus, wave_S,
    ).value ≈ sqrt(2)

    identical = PhysicalState("B", 0.4, [
        (basis = first_basis[1], coefficient = inv(sqrt(2)), wave = wave),
        (basis = first_basis[2], coefficient = im / sqrt(2), wave = wave),
    ])
    identical_reordered_rephased = PhysicalState("B", 0.4, [
        (basis = first_basis[2], coefficient = -inv(sqrt(2)), wave = wave),
        (basis = first_basis[1], coefficient = im / sqrt(2), wave = wave),
    ])
    identical_final = TwoMesonChannel(identical, identical_reordered_rephased)
    @test QuarkModelTransitions._same_external_state(
        identical_final.first, identical_final.second,
    )
    relative_rephased = PhysicalState("B", 0.4, [
        (basis = first_basis[1], coefficient = inv(sqrt(2)), wave = wave),
        (basis = first_basis[2], coefficient = -im / sqrt(2), wave = wave),
    ])
    @test !QuarkModelTransitions._same_external_state(identical, relative_rephased)
    @test allowed_partial_waves(identical_final, pure_first) == [wave_S]
    identical_composed = QuarkModelTransitions._compose_physical_decay(
        equal_kernel, identical_final, pure_first, wave_S,
    )
    raw_identical = sum(
        a.coefficient * conj(b.coefficient) * conj(c.coefficient)
        for a in pure_first.components for b in identical_final.first.components for
        c in identical_final.second.components
    )
    @test identical_composed.value ≈ raw_identical / sqrt(2)
    vector_parent_basis = BasisState(1, "S", 3, 1; label = "V", flavors = (:u, :u))
    vector_parent = PhysicalState("V", 1.8,
        [(basis = vector_parent_basis, coefficient = 1.0, wave = wave)])
    @test isempty(allowed_partial_waves(identical_final, vector_parent))
end
