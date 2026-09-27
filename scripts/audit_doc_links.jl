#!/usr/bin/env julia
# Audit explicit @ref links in core API docstrings (not Julia call dependencies)
# and regenerate docs/discoverability_graph.md. Stdlib only.
#
# Run from any directory: julia scripts/audit_doc_links.jl

const ROOT = dirname(@__DIR__)
const CORE = split("""load_parameters load_quark_masses load_parameters_and_quark_masses
GIParameters QuarkMassTable ConstituentMasses Meson reduced_mass flavor_label
is_equal_flavor BasisState spectrum_levels
compute_spectrum MixedSpectrum MixedState StateMixing spectrum_state
radial_wave physical_components RadialWave sample_wave wave_norm
physical_state_amplitude compute_isoscalar_spectrum""")
const DOC = r"\"\"\"(.*?)\"\"\""s
# Anchored at the end of the docstring: skip comments, then read the documented name.
const DECL = r"\G\s*(?:#[^\n]*\n\s*)*(?:(?:mutable\s+)?struct\s+|abstract\s+type\s+|function\s+|const\s+)?([A-Za-z_]\w*)"
const LINK = r"\[`([^`]+)`\]\(@ref\)"

function collect_links()
    docs = Dict{String,Vector{String}}()
    srcdir = joinpath(ROOT, "src")
    for file in sort(filter(endswith(".jl"), readdir(srcdir)))
        text = read(joinpath(srcdir, file), String)
        for m in eachmatch(DOC, text)
            declaration = match(DECL, text, m.offset + ncodeunits(m.match))
            isnothing(declaration) && continue
            push!(get!(docs, declaration[1], String[]), m[1])
        end
    end
    return Dict(name => setdiff(Set(String(l[1]) for l in eachmatch(LINK, join(bodies, '\n'))), [name])
                for (name, bodies) in docs)
end

function metrics(graph)
    edges = Set((a, b) for a in CORE for b in get(graph, a, Set{String}()) if b in CORE)
    incoming = Dict(n => count(e -> e[2] == n, edges) for n in CORE)
    outgoing = Dict(n => count(e -> e[1] == n, edges) for n in CORE)
    return edges, incoming, outgoing
end

function reachable(graph, start)
    seen, todo = Set([start]), [start]
    while !isempty(todo)
        for target in setdiff(intersect(get(graph, pop!(todo), Set{String}()), CORE), seen)
            push!(seen, target)
            push!(todo, target)
        end
    end
    return seen
end

function main()
    graph = collect_links()
    edges, incoming, outgoing = metrics(graph)
    lines = ["# Public API documentation graph", "",
        "Explore $(length(CORE)) public entries along the input → spectrum → wavefunction workflow.",
        "Arrows lead from a help entry to a related entry linked in its docstring.",
        "Start with the suggested reading paths, or expand the complete graph.", ""]
    route = [
        ("QuarkMassTable", "load_quark_masses"),
        ("load_quark_masses", "load_parameters_and_quark_masses"),
        ("load_parameters", "load_parameters_and_quark_masses"),
        ("load_parameters_and_quark_masses", "Meson"),
        ("ConstituentMasses", "Meson"), ("Meson", "reduced_mass"),
        ("Meson", "flavor_label"), ("Meson", "is_equal_flavor"),
        ("Meson", "spectrum_levels"),
        ("spectrum_levels", "BasisState"), ("Meson", "compute_spectrum"),
        ("compute_spectrum", "MixedSpectrum"), ("MixedSpectrum", "MixedState"),
        ("MixedSpectrum", "spectrum_state"), ("spectrum_state", "physical_components"),
        ("MixedState", "StateMixing"), ("physical_components", "RadialWave"),
        ("physical_components", "physical_state_amplitude"),
        ("spectrum_state", "radial_wave"), ("radial_wave", "sample_wave"),
        ("radial_wave", "wave_norm"),
        ("BasisState", "compute_isoscalar_spectrum"),
    ]
    append!(lines, ["", "## Suggested reading paths", "",
        "A selected subset of actual links makes the main workflow easier to read.",
        "", "```mermaid", "flowchart LR"])
    append!(lines, ["    $a --> $b" for (a, b) in route if (a, b) in edges])
    append!(lines, ["```", "", "## Complete graph", "",
        "<details>", "<summary>Show all core documentation links</summary>", "",
        "```mermaid", "flowchart LR"])
    append!(lines, ["    $n[\"$n\"]" for n in CORE])
    append!(lines, ["    $a --> $b" for (a, b) in sort!(collect(edges))])
    append!(lines, ["```", "", "</details>", "", "## Link coverage", "",
        "The selected entries contain **$(length(edges))** links. From `Meson`,",
        "**$(length(reachable(graph, "Meson"))) of $(length(CORE))** entries are reachable (including itself).", "",
        "| Entry | Incoming links | Outgoing links |", "|---|---:|---:|"])
    append!(lines, ["| `$n` | $(incoming[n]) | $(outgoing[n]) |" for n in CORE])
    append!(lines, ["",
        "Counts combine constructor and method docstrings and exclude self-links.",
        "This graph covers explicit docstring links among the selected entries;",
        "it does not represent function calls or every exported method.", "",
        "See [Discoverability](discoverability.md) for help conventions and regeneration.", ""])
    output = joinpath(ROOT, "docs", "discoverability_graph.md")
    write(output, join(lines, '\n'))
    println("$output: $(length(edges)) core links")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
