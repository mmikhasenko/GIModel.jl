@testset "reference loading" begin
    ccbar = load_ref(reference_spectrum_path("charmonium"))
    bbbar = load_ref(reference_spectrum_path("bottomonium"))
    @test length(ccbar) == 28
    @test length(bbbar) == 30
    @test ccbar[1].quark_content == "c cbar"
    @test ccbar[1].composition == "1^1S_0"
    @test_throws ArgumentError reference_spectrum_path("nonexistent_sector")
end

@testset "reference_meson resolution (loud, no fallback)" begin
    row(sector, content) = GIPaper.ReferenceState(sector, content, "1^1S_0", 1, 1, "S", 0, 1.0, "high")
    m = reference_meson(mq, row("ccbar", "ignore"))
    @test m.constituent_masses.m1_GeV ≈ m.constituent_masses.m2_GeV ≈ mq["c"]
    m = reference_meson(mq, row("charmonium", "c cbar"))
    @test (m.flavor1, m.flavor2) == (:c, :c)
    # u/d/n/q all resolve to the one LightQuark the model can represent, so the
    # light leg reads `:q` regardless of how the CSV spells it. Masses are
    # unchanged (m_u = m_d = m_q in the parameter set).
    m = reference_meson(mq, row("charmed", "-c dbar; c ubar"))
    @test (m.flavor1, m.flavor2) == (:c, :q)
    @test m.constituent_masses.m1_GeV ≈ mq["c"]
    @test m.constituent_masses.m2_GeV ≈ mq["d"]
    m = reference_meson(mq, row("charmed_strange", "c sbar"))
    @test (m.flavor1, m.flavor2) == (:c, :s)
    m = reference_meson(mq, row("bottom_light", "b ubar; -b dbar"))
    @test (m.flavor1, m.flavor2) == (:b, :q)
    m = reference_meson(mq, row("strange", "-u sbar; -d sbar"))
    @test (m.flavor1, m.flavor2) == (:q, :s)
    m = reference_meson(mq, row("isoscalar", "n nbar / s sbar mixed"))
    @test (m.flavor1, m.flavor2) == (:q, :q)
    @test m.constituent_masses.m1_GeV ≈ 0.5 * (mq["u"] + mq["d"])
    # failures are loud and name the offending row
    @test_throws ArgumentError reference_meson(mq, row("unknown_sector", "c cbar"))
    @test_throws ArgumentError reference_meson(mq, row("charmed", "gibberish"))
    @test_throws ArgumentError reference_meson(mq, row("charmed", "x ybar"))
    @test_throws ArgumentError reference_meson(mq, row("charmed", "c d"))  # antiquark must end in bar
end
