@testset "GI Appendix-B spin algebra" begin
    q = QuarkModelTransitions._QuarkEmission()
    qbar = QuarkModelTransitions._AntiquarkEmission()
    singlet = QuarkModelTransitions._coupled_spin(0, 0)
    triplet = Dict(m => QuarkModelTransitions._coupled_spin(1, m) for m in -1:1)

    @test sum(abs2, values(singlet.components)) ≈ 1
    @test all(state -> sum(abs2, values(state.components)) ≈ 1, values(triplet))
    @test QuarkModelTransitions._spin_factor(q, singlet, triplet[0], 0) ≈ 1
    @test QuarkModelTransitions._spin_factor(qbar, singlet, triplet[0], 0) ≈ -1
    @test QuarkModelTransitions._spin_factor(q, singlet, triplet[-1], 1) ≈ -1
    @test QuarkModelTransitions._spin_factor(qbar, singlet, triplet[-1], 1) ≈ 1
    @test QuarkModelTransitions._spin_factor(q, singlet, triplet[1], -1) ≈ -1
    @test QuarkModelTransitions._spin_factor(qbar, singlet, triplet[1], -1) ≈ 1

    # A rank-one component changes the total projection by its spherical index.
    for topology in (q, qbar), mf in -1:1, mi in -1:1, mu in -1:1
        value = QuarkModelTransitions._spin_factor(
            topology, triplet[mf], triplet[mi], mu,
        )
        mf == mi + mu || @test value == 0
    end

    # Spherical-tensor Hermiticity: sigma_mu^dagger = (-1)^mu sigma_-mu.
    states = (singlet, triplet[-1], triplet[0], triplet[1])
    for topology in (q, qbar), left in states, right in states, mu in -1:1
        lhs = QuarkModelTransitions._spin_factor(topology, left, right, mu)
        rhs = (-1)^mu * conj(QuarkModelTransitions._spin_factor(
            topology, right, left, -mu,
        ))
        @test lhs ≈ rhs
    end

    @test_throws ArgumentError QuarkModelTransitions._coupled_spin(2, 0)
    @test_throws ArgumentError QuarkModelTransitions._pauli_spherical(1, 2)
end
