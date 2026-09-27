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

    # A rank-one component changes projection by its spherical index.
    # Compare the operators as matrices rather than testing every scalar entry
    # repeatedly (most entries are the same selection-rule zero).
    states = (singlet, triplet[-1], triplet[0], triplet[1])
    projections = (0, -1, 0, 1)
    for topology in (q, qbar)
        matrix(mu) = [QuarkModelTransitions._spin_factor(topology, l, r, mu)
                      for l in states, r in states]
        z, plus, minus = matrix(0), matrix(1), matrix(-1)
        @test z ≈ z'
        @test plus ≈ -minus'
        @test all(iszero(plus[i,j]) for i in 1:4, j in 1:4
                  if projections[i] != projections[j] + 1)
    end
end
