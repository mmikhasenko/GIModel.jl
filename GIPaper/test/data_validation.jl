function csv_rows(path)
    isfile(path) || return NamedTuple[]
    return collect(CSV.File(path; normalizenames = false))
end

function cell(row, name)
    value = getproperty(row, Symbol(name))
    return ismissing(value) ? "" : strip(string(value))
end

function has_columns(rows, required)
    isempty(rows) && return false
    names = Set(String.(propertynames(first(rows))))
    return issubset(Set(required), names)
end

function provenance_exists(relative)
    isempty(relative) && return false
    return ispath(joinpath(GIPAPER_ROOT, relative)) ||
           ispath(joinpath(REPOSITORY_ROOT, relative))
end

@testset "paper data invariants" begin
    seed_path = joinpath(GIPAPER_ROOT, "data", "seed", "godfrey_isgur_seed_masses.csv")
    seed = csv_rows(seed_path)
    required = (
        "sector", "state_label", "assignment", "mass_MeV",
        "source_status", "source_short", "source_url",
    )
    @test has_columns(seed, required)
    @test !isempty(seed)
    keys = Tuple{String,String,String,String}[]
    for row in seed
        @test !isempty(cell(row, "sector"))
        @test !isempty(cell(row, "assignment"))
        @test !isempty(cell(row, "source_status"))
        @test tryparse(Float64, cell(row, "mass_MeV")) !== nothing
        push!(keys, (
            cell(row, "sector"),
            cell(row, "state_label"),
            cell(row, "assignment"),
            cell(row, "mass_MeV"),
        ))
    end
    @test length(keys) == length(unique(keys))

    reference_required = (
        "sector", "quark_content", "composition_raw", "n",
        "multiplicity", "L", "J", "mass_GeV", "confidence",
    )
    reference_files = filter(
        path -> startswith(basename(path), "reference_spectrum_") &&
                endswith(path, ".csv"),
        readdir(joinpath(GIPAPER_ROOT, "data"); join = true),
    )
    @test !isempty(reference_files)
    for path in reference_files
        rows = csv_rows(path)
        @test !isempty(rows)
        @test has_columns(rows, reference_required)
    end
end

@testset "Table II agrees with model parameters" begin
    rows = csv_rows(joinpath(GIPAPER_ROOT, "data", "table_ii_parameters.csv"))
    table_values = Dict(cell(row, "parameter_key") => parse(Float64, cell(row, "value"))
                        for row in rows)
    model = TOML.parsefile(default_parameters_path())
    mappings = (
        "m_ud_avg" => ("masses", "m_ud_avg_MeV"),
        "m_s" => ("masses", "m_s_MeV"),
        "m_c" => ("masses", "m_c_MeV"),
        "m_b" => ("masses", "m_b_MeV"),
        "b" => ("potential", "b_GeV2"),
        "Lambda" => ("potential", "Lambda_MeV"),
        "c" => ("potential", "c_MeV"),
        "sigma0" => ("relativistic_smearing", "sigma0_GeV"),
        "s" => ("relativistic_smearing", "s"),
        "epsilon_c" => ("relativistic_factors", "epsilon_c"),
        "epsilon_t" => ("relativistic_factors", "epsilon_t"),
        "epsilon_so_vector" => ("relativistic_factors", "epsilon_so_vector"),
        "epsilon_so_scalar" => ("relativistic_factors", "epsilon_so_scalar"),
    )
    for (key, path) in mappings
        @test haskey(table_values, key)
        @test table_values[key] ≈ Float64(model[path[1]][path[2]]) rtol = 1e-9
    end
end

@testset "promoted data is internally consistent" begin
    valid_confidence = Set(("high", "medium", "low"))

    masses = csv_rows(joinpath(GIPAPER_ROOT, "data", "clean", "masses.csv"))
    @test !isempty(masses)
    mass_ids = String[]
    for row in masses
        clean_id = cell(row, "clean_id")
        mass = tryparse(Float64, cell(row, "mass_model_MeV"))
        @test !isempty(clean_id)
        @test mass !== nothing
        @test mass === nothing || mass > 0
        @test provenance_exists(cell(row, "provenance_file"))
        @test cell(row, "confidence") in valid_confidence
        push!(mass_ids, clean_id)
    end
    @test length(mass_ids) == length(unique(mass_ids))

    mixings = csv_rows(joinpath(GIPAPER_ROOT, "data", "clean", "mixings.csv"))
    @test !isempty(mixings)
    mixing_ids = String[]
    for row in mixings
        clean_id = cell(row, "clean_id")
        @test !isempty(clean_id)
        @test tryparse(Float64, cell(row, "amplitude")) !== nothing
        @test provenance_exists(cell(row, "provenance_file"))
        @test cell(row, "confidence") in valid_confidence
        push!(mixing_ids, clean_id)
    end
    @test length(mixing_ids) == length(unique(mixing_ids))
end
