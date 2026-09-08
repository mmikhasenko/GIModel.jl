@testset "GIModel package boundary" begin
    @test model_parameters_path() == default_parameters_path()

    source_files = filter(
        path -> endswith(path, ".jl"),
        readdir(joinpath(GIPAPER_ROOT, "src"); join = true),
    )
    qualified_names = Set{Symbol}()
    for path in source_files
        for match in eachmatch(
            r"GIModel\.([A-Za-z_][A-Za-z_0-9!]*)",
            read(path, String),
        )
            push!(qualified_names, Symbol(only(match.captures)))
        end
    end
    @test all(name -> Base.isexported(GIModel, name), qualified_names)
end
