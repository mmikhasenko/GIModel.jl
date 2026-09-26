#!/usr/bin/env julia
# Read-only freshness and canonical coverage gate. No spectrum computation.
using Pkg
Pkg.activate(@__DIR__; io=devnull)
using TOML, SHA, CSV
const ROOT = dirname(dirname(@__DIR__))
const DIR = joinpath(dirname(@__DIR__), "docs", "input_traces")
common = TOML.parsefile(joinpath(DIR, "common.toml"))
common["schema_version"] == 1 || error("unsupported trace schema")
for (path, expected) in common["source_sha256"]
    file = joinpath(ROOT, path)
    isfile(file) || error("trace source missing: $path")
    bytes2hex(sha256(read(file))) == expected || error("stale input trace: $path changed; regenerate trace_rate_inputs.jl")
end
for dir in ("src", "QuarkModelTransitions/src", "GIPaper/src", "GIPaper/data/mass_inputs")
    for (base, _, names) in walkdir(joinpath(ROOT,dir)), name in names
        path = relpath(joinpath(base,name),ROOT)
        haskey(common["source_sha256"], path) || error("unrecorded source: $path")
    end
end
for (name, expected) in common["trace_sha256"]
    bytes2hex(sha256(read(joinpath(DIR,name)))) == expected || error("trace result changed: $name")
end
haskey(common, "paper_table_policy") || error("missing paper-policy ownership record")
for (table, csv, rowkey) in (("v", "table_v_strong_decays.csv", "canonical_rows"),
                            ("vi", "table_vi_photon_decays.csv", "rows"),
                            ("vii", "table_vii_annihilation_em.csv", "canonical_rows"))
    trace = TOML.parsefile(joinpath(DIR, "table_$table.toml"))
    trace["common_inputs"] == "common.toml" || error("unresolved common inputs")
    canonical = CSV.File(joinpath(dirname(@__DIR__), "data", "raw", "digitized_tables", csv))
    rows = trace[rowkey]
    length(rows) == length(canonical) || error("Table $table row coverage mismatch")
    sort([r["decay"] for r in rows]) == sort([String(r.decay) for r in canonical]) || error("Table $table row identities mismatch")
    if table == "vii"
        sum(r["implemented"] for r in rows) == sum(length(trace[k]) for k in ("gluonic", "leptonic", "two_photon", "mixed_two_photon", "charge_radius")) || error("Table VII result coverage mismatch")
    end
    println("Table ", uppercase(table), ": ", length(rows), " canonical rows traced")
end
println("Input traces: source fingerprints, ownership record and row coverage verified")
