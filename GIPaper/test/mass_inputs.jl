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
