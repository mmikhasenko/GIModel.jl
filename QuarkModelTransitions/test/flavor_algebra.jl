@testset "GI Appendix-B sparse flavor algebra" begin
    names = (
        :pi_plus, :pi_zero, :pi_minus, :K_plus, :K_zero, :Kbar_zero,
        :K_minus, :eta8, :eta1, :M_ns, :M_s, :eta, :eta_prime,
    )
    states = Dict(name => QuarkModelTransitions._appendix_b_flavor(name) for name in names)
    @test all(
        state -> sum(abs2, values(state.components)) ≈ 1.0,
        values(states),
    )

    q = QuarkModelTransitions._QuarkEmission()
    qbar = QuarkModelTransitions._AntiquarkEmission()
    @test QuarkModelTransitions._flavor_transfer(q, states[:pi_plus], :d, :u) ≈ -sqrt(2)
    @test QuarkModelTransitions._flavor_transfer(q, states[:pi_zero], :u, :u) ≈ 1
    @test QuarkModelTransitions._flavor_transfer(q, states[:pi_zero], :d, :d) ≈ -1
    @test QuarkModelTransitions._flavor_transfer(qbar, states[:pi_plus], :u, :d) ≈ sqrt(2)
    @test QuarkModelTransitions._flavor_transfer(q, states[:eta], :u, :u) ≈ inv(sqrt(2))
    @test QuarkModelTransitions._flavor_transfer(q, states[:eta], :s, :s) ≈ -1

    # Both topologies obey the de Swart phases for rho0 -> pi+ pi-.
    @test QuarkModelTransitions._flavor_factor(
        q, states[:pi_minus], states[:pi_plus], states[:pi_zero],
    ) ≈ -1
    @test QuarkModelTransitions._flavor_factor(
        qbar, states[:pi_minus], states[:pi_plus], states[:pi_zero],
    ) ≈ -1

    # Spectator mismatch and OZI-forbidden pure-ss -> pi pi transitions are exact zeros.
    @test QuarkModelTransitions._flavor_factor(
        q, states[:pi_minus], states[:pi_plus], states[:M_s],
    ) == 0
    @test QuarkModelTransitions._flavor_factor(
        qbar, states[:pi_minus], states[:pi_plus], states[:M_s],
    ) == 0

    # The physical f' fixture is zero for pure ss and nonzero through an nn admixture.
    angle = 0.12
    fprime = QuarkModelTransitions._FlavorState("f'", [
        (:s, :s) => cos(angle),
        (:u, :u) => sin(angle) / sqrt(2),
        (:d, :d) => sin(angle) / sqrt(2),
    ])
    pure = QuarkModelTransitions._flavor_factor(
        q, states[:pi_minus], states[:pi_plus], states[:M_s],
    )
    mixed = QuarkModelTransitions._flavor_factor(
        q, states[:pi_minus], states[:pi_plus], fprime,
    )
    @test pure == 0
    @test mixed ≈ -sin(angle)

    # Overall flavor-state phases transform a bra/operator/ket matrix element covariantly.
    phases = (cis(0.2), cis(-0.4), cis(0.7))
    rephase(state, phase) = QuarkModelTransitions._FlavorState(
        state.label,
        [key => phase * value for (key, value) in state.components],
    )
    base = QuarkModelTransitions._flavor_factor(
        q, states[:pi_minus], states[:pi_plus], states[:pi_zero],
    )
    changed = QuarkModelTransitions._flavor_factor(
        q,
        rephase(states[:pi_minus], phases[1]),
        rephase(states[:pi_plus], phases[2]),
        rephase(states[:pi_zero], phases[3]),
    )
    @test changed ≈ conj(phases[1]) * conj(phases[2]) * phases[3] * base

    @test QuarkModelTransitions._heavy_flavor_state(:c, :d).components ==
          Dict((:c, :d) => -1.0)
end
