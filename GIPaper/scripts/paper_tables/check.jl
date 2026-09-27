#!/usr/bin/env julia
# Consistency checks on the generated paper-layout objects and the canonical
# Table III/V/VI/VII reports they reference. Run after generate.jl:
#
#   julia GIPaper/scripts/paper_tables/check.jl
using Test, TOML
const GIPAPER_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const DATA = joinpath(GIPAPER_DIR, "docs", "paper_tables")
const REPORTS = joinpath(GIPAPER_DIR, "docs", "residual_reports")
function tsv(name)
    lines = readlines(joinpath(DATA, name * ".tsv"))
    header = split(first(lines), '\t')
    [Dict(header .=> split(line, '\t'; keepempty=true)) for line in lines[2:end]]
end
cells(line) = strip.(split(line, '|')[2:end-1])
@testset "Requested numerical quantities" begin
    confinement = tsv("table_i")
    @test length(confinement) == 5
    @test all(parse(Float64, row["computed_splitting_GeV"]) > 0 for row in confinement)
    @test all(0 < parse(Float64, row["computed_confinement_percent"]) < 100
        for row in confinement)
    @test issorted(parse(Float64, row["computed_confinement_percent"])
        for row in confinement; rev=true)
    parameters = tsv("table_ii")
    @test length(parameters) == 14
    critical = only(row for row in parameters if row["parameter"] == "alpha_s_critical")
    @test parse(Float64, critical["active_value"]) ≈ 0.60
    @test parse(Float64, critical["gi_value"]) == 0.60
    coupling = tsv("fig_i")
    for row in coupling
        M = parse(Float64,row["M_GeV"])
        expected = .25exp(-M^2) + .15exp(-M^2/10) + .20exp(-M^2/1000)
        @test parse(Float64,row["alpha_s_M"]) ≈ expected
    end
    mix = tsv("caption_mixing")
    @test Set(r["figure"] for r in mix) == Set(("iii", "iv", "vi"))
    for figure in ("iii", "iv", "vi")
        rows = filter(r -> r["figure"] == figure, mix)
        @test length(rows) == 6
        @test sum(parse(Float64,r["probability"]) for r in rows) ≈ 1
        for row in rows
            @test parse(Float64,row["probability"]) ≈ abs2(parse(Float64,row["amplitude"]))
        end
    end
    report = readlines(joinpath(REPORTS,"table_iii_mixing_audit.md"))
    headers = String[]
    ncomponents = 0
    for line in report
        startswith(line,"|") || continue
        cs = cells(line)
        if first(cs) == "state"
            headers = cs
        elseif startswith(first(cs),"`") && any(==("1 ns model"), headers)
            indices = findall(h -> occursin(r"^[12] (ns|ss|cc|bb) model$",h),headers)
            vector = parse.(Float64,cs[indices])
            @test isapprox(sum(abs2,vector),1; atol=3e-6)
            ncomponents += length(vector)
        end
    end
    @test ncomponents == 90
    columns = TOML.parsefile(joinpath(DATA,"quantity_columns.toml"))["column"]
    @test only(q["role"] for q in columns if q["column"] == "M_P (GeV)") == "physical mass input"
    @test only(q["role"] for q in columns if q["column"] == "q_eff") == "algebraic input"
    @test all(isfile(joinpath(GIPAPER_DIR,q["source"])) for q in columns)
    report7 = readlines(joinpath(REPORTS,"table_vii_annihilation_em.md"))
    report7text = join(report7, '\n')
    for label in ("zeta -> e+ e-", "eta_t -> gamma gamma", "eta_t -> 2g", "zeta -> 3g")
        @test occursin("`$label`", report7text)
    end
    @test occursin("All 61 canonical", report7text)
    widths = false
    registry=TOML.parsefile(joinpath(DATA,"mass_inputs.toml"))["inputs"]
    masses=Dict(r["label"]=>r["mass_GeV"] for r in registry if r["context"]=="VII" && haskey(r,"mass_GeV"))
    checked = 0
    for line in report7
        startswith(line,"| decay | f_model |") && (widths=true; continue)
        widths && startswith(line,"| `") || continue
        cs = cells(line)
        label = replace(cs[1],"`"=>"")
        haskey(masses,label) || continue
        f = parse(Float64,cs[2]); width = parse(Float64,first(split(cs[3])))
        expected = 4pi/3 * (1/137.036)^2 * masses[label] * f^2
        # Reports print f to four decimals; propagate its rounding error.
        @test isapprox(width,expected; rtol=0.006)
        checked += 1
    end
    @test checked == 4
    report5 = read(joinpath(REPORTS,"table_v_reproduction.md"), String)
    comparison, diagnoses = split(report5, "\n## Investigated misses"; limit=2)
    @test count(line -> startswith(line, "| `"), split(comparison, '\n')) == 220
    @test count(line -> startswith(line, "| `"), split(diagnoses, '\n')) == 13
end

@testset "Complete displayed photon kinematics" begin
    rows = [cells(l) for l in readlines(joinpath(REPORTS,"table_vi_photon_decays.md"))
        if startswith(l,"|") && occursin("gamma",l)]
    @test length(rows) == 79
    @test all(isnan(parse(Float64,r[4])) || parse(Float64,r[4]) >= 0 for r in rows)
    @test count(r -> r[5] == "—", rows) == 9
end
