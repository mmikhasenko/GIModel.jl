@testset "Canonical Table VI identities and predictions" begin
    rows = load_table_vi()
    @test length(rows) == 79
    refs = Dict(r.id => r for r in rows)
    # Adjacent hindered rows have distinct final radial levels and signs.
    @test refs[(:M1, "Upsilondoubleprime", "etaprime_b")].predicted == 0.007
    @test refs[(:M1, "Upsilondoubleprime", "eta_b")].predicted == -0.004
    @test refs[(:M2, "A2", "pi")].formula == "q^2/(sqrt(60) m_u) E1^u(pi,A2)"
    @test refs[(:M1, "Upsilon", "pi0")].approximate
    @test refs[(:M1, "Upsilon", "pi0")].parenthetical
    @test refs[(:M1, "psi", "eta")].predicted == 0.001
    @test count(r -> first(r.id) == :M1, rows) == 42
    @test count(r -> first(r.id) == :E1, rows) == 35
    @test count(r -> first(r.id) == :M2, rows) == 2
    # Reordering CSV rows cannot change transition/value association.
    mktemp() do path, io
        source = readlines(joinpath(GIPaper.paper_data_dir(), "raw",
            "digitized_tables", "table_vi_photon_decays.csv"))
        println(io, first(source))
        foreach(line -> println(io, line), reverse(source[2:end]))
        close(io)
        @test Dict(r.id => r.predicted for r in load_table_vi(path)) ==
              Dict(r.id => r.predicted for r in rows)
    end
    mktemp() do path, io
        source = readlines(joinpath(GIPaper.paper_data_dir(), "raw",
            "digitized_tables", "table_vi_photon_decays.csv"))
        foreach(line -> println(io, line), source)
        println(io, source[2])
        close(io)
        @test_throws ArgumentError load_table_vi(path)
    end
end

@testset "Complete Table VI state mapping" begin
    rows = load_table_vi()
    states = load_table_vi_states()
    @test Set(keys(states)) == Set(id for r in rows for id in r.id[2:3])
    @test states["Bstarminus"].flavors == (:u, :b) # printed historical label
    @test states["etaprime_b"].basis.n == 2
    @test states["chiprime_2b"].basis.n == 2
    @test states["eta_r"].mixed
    @test only(r for r in rows if r.decay == "A1 -> pi gamma").id[1] == :E1
    for row in rows
        p, d = states[row.id[2]], states[row.id[3]]
        if first(row.id) != :M1 || any(f in split(row.footnotes, ",") for f in ("c", "g"))
            if !isnothing(p.mass_GeV) && !isnothing(d.mass_GeV)
                @test p.mass_GeV > d.mass_GeV
            else
                @test p.mass_source == "unavailable" || d.mass_source == "unavailable"
            end
        end
    end
end
