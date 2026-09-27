include(joinpath(@__DIR__, "common.jl"))
using TOML
using CSV

function compute_mass_inputs()
    historical = Dict{Tuple{String,String},Vector{String}}()
    for r in CSV.File(joinpath(GIPAPER_DIR,"data","mass_inputs","historical.csv"); types=String)
        push!(get!(historical,(String(r.context),String(r.label)),String[]),
            "$(r.mass_GeV) GeV — $(r.source); $(r.classification)")
    end
    rows = Dict[]
    for ((context,label), r) in sort!(collect(load_mass_inputs()); by=first)
        d = Dict{String,Any}(string(k)=>v for (k,v) in pairs(r) if !isnothing(v))
        d["historical"] = get(historical,(context,label),String[])
        push!(rows,d)
    end
    path=joinpath(PAPER_TABLES_DIR,"mass_inputs.toml")
    open(path,"w") do io
        TOML.print(io,Dict("edition"=>"PDG 2026", "inputs"=>rows);sorted=true)
    end
    return path
end
abspath(PROGRAM_FILE) == (@__FILE__) && compute_mass_inputs()
