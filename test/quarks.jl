@testset "Quark types (mass is dynamics, charge is the only discrete datum)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    c = HeavyQuark{:up}(mq["c"], :c)
    b = HeavyQuark{:down}(mq["b"], :b)
    s = StrangeQuark(mq["s"])
    l = LightQuark(mq["q"])

    # Charge is the discrete flavor datum, carried by the weak-isospin tag.
    @test charge(c) == 2 // 3
    @test charge(b) == -1 // 3
    @test charge(s) == -1 // 3

    # LightQuark has NO charge: the model sets m_u = m_d, so it cannot resolve
    # up from down. Charge lives on the meson's flavor wavefunction instead
    # (pi+ = u dbar is charged, pi0 = (u ubar - d dbar)/sqrt(2) is not, at the
    # same constituent masses). This must raise, never guess.
    @test_throws MethodError charge(l)

    # The tag is weak-isospin class, validated at construction.
    @test_throws ArgumentError HeavyQuark{:sideways}(1.5, :x)
    @test_throws ArgumentError HeavyQuark{:up}(-1.0, :x)
    @test_throws ArgumentError LightQuark(0.0)
    @test_throws ArgumentError StrangeQuark(-0.4)

    @test flavor_symbol(l) === :q
    @test flavor_symbol(s) === :s
    @test flavor_symbol(c) === :c
    @test mass_GeV(c) == mq["c"]

    # Quark-built mesons agree with the flavor-table form.
    @test Meson(c, c) == Meson(mq, :c, :c)
    @test Meson(l, s) == Meson(mq, :q, :s)

    # Mass is the only dynamical input: same mass + different charge class ->
    # identical constituent masses, different charge.
    up_at_mb = HeavyQuark{:up}(mq["b"], :hypothetical)
    @test Meson(up_at_mb, up_at_mb).constituent_masses ==
          Meson(b, b).constituent_masses
    @test charge(up_at_mb) != charge(b)
end
