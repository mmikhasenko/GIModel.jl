@testset "Experimental mass registry and decay kinematics" begin
    inputs=load_mass_inputs()
    @test all(isnothing(r.mass_GeV) || (isfinite(r.mass_GeV) && r.mass_GeV>0 && r.edition=="2026") for r in values(inputs))
    @test_throws KeyError experimental_mass("VII","unknown state")
    @test isnothing(experimental_mass("VII","1^3D_1(bb) -> e+ e-"))
    @test isnothing(experimental_mass("VI","etadoubleprime_b"))
    @test experimental_mass("VII","pi -> gamma gamma") ≈ 0.1349768277676847
    @test experimental_mass("VII","pi -> mu nu") ≈ 0.1395703909836813
    @test experimental_mass("V","pi") ≈ (2experimental_mass("VII","pi -> mu nu") + experimental_mass("VII","pi -> gamma gamma"))/3
    @test experimental_mass("V","eta_r") == experimental_mass("VI","eta_r") == experimental_mass("VII","eta_r -> gamma gamma")
    @test experimental_mass("VI","psi") == experimental_mass("VII","psi -> e+ e-") == experimental_mass("VII","psi -> 3g")
    states=load_table_vi_states()
    for (key,r) in states
        @test r.mass_GeV === experimental_mass("VI",key)
    end
    # Saved audit: verify ALL photon kinematics against the independent two-body formula.
    report=joinpath(GIPaper.paper_data_dir(),"..","docs","residual_reports","table_vi_photon_decays.csv")
    rows=CSV.File(report)
    @test length(rows)==79
    for r in rows
        M=experimental_mass("VI",String(r.parent)); m=experimental_mass("VI",String(r.daughter))
        if isnothing(M) || isnothing(m)
            @test ismissing(r.q_GeV)
            @test isnan(r.q_display_GeV)
        else
            @test r.q_GeV ≈ (M*M-m*m)/(2M)
            @test r.q_GeV == r.q_display_GeV
        end
    end
end

@testset "Mass registry rejects ambiguous or invalid data" begin
    for (file, line) in (("pdg-2026.csv", "duplicate"), ("averages.csv", "bad,unknown_component"), ("assignments.csv", "TEST,unknown,not_registered,test"))
        mktempdir() do dir
            source = joinpath(GIPaper.paper_data_dir(), "mass_inputs")
            for name in ("pdg-2026.csv", "averages.csv", "assignments.csv")
                cp(joinpath(source, name), joinpath(dir, name))
            end
            addition = line == "duplicate" ? readlines(joinpath(dir, file))[2] : line
            open(joinpath(dir, file), "a") do io
                println(io, addition)
            end
            @test_throws ArgumentError load_mass_inputs(dir)
        end
    end
end

@testset "Paper policy provenance and referential integrity" begin
    policy = load_table_policy()
    canonical_v = CSV.File(joinpath(GIPaper.paper_data_dir(), "raw", "digitized_tables", "table_v_strong_decays.csv"))
    for row in policy["table_v"]["mixing"]
        partners = Set(String(r.parent) for r in canonical_v if r.section == row["section"])
        @test row["singlet"] in partners
        @test row["triplet"] in partners
    end
    states = load_table_vi_states()
    @test all(haskey(states, s) for s in policy["table_vi"]["isovectors"])
    @test policy["table_vii"]["validation"]["psi_dilepton_width_GeV"] ==
        only(r["width_GeV"] for r in policy["table_vii"]["dilepton_validation"] if r["label"] == "psi -> e+ e-")
    @test_throws ArgumentError load_parameters(joinpath(GIPaper.paper_data_dir(), "clean", "parameters.toml"))
    @test load_quark_masses(joinpath(GIPaper.paper_data_dir(), "clean", "parameters.toml")) == mq
    mktempdir() do dir
        push!(policy["table_v"]["mixing"], first(policy["table_v"]["mixing"]))
        path = joinpath(dir, "policy.toml")
        open(io -> TOML.print(io, policy), path, "w")
        @test_throws ArgumentError load_table_policy(path)
    end
end
