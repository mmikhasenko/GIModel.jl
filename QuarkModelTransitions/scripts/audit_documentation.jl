# Execute every public docstring example and check help links. Run from the
# repository root; as a script it also prints the documentation link graph:
# julia --project=. QuarkModelTransitions/scripts/audit_documentation.jl
using GIModel, GIModel.QuarkModelTransitions
using Test

const PACKAGE_ROOT = dirname(@__DIR__)
const QMT = QuarkModelTransitions
const PUBLIC_NAMES = Set(string(n) for n in names(QMT) if n != :QuarkModelTransitions)
const DOC_PATTERN = r"\"\"\"(.*?)\"\"\""s
const DECL_PATTERN = r"^\s*(?:#[^\n]*\n\s*)*(?:(?:mutable\s+)?struct\s+|abstract\s+type\s+|function\s+|const\s+)?([A-Za-z_]\w*)"
const LINK_PATTERN = r"\[`([^`]+)`\]\(@ref\)"
const EXAMPLE_PATTERN = r"```julia\s*\n(.*?)```"s

function public_docstrings()
    entries = Tuple{String,String,String}[]
    for file in sort(readdir(joinpath(PACKAGE_ROOT, "src"); join=true))
        endswith(file, ".jl") || continue
        source = read(file, String)
        for doc in eachmatch(DOC_PATTERN, source)
            tail = SubString(source, doc.offset + ncodeunits(doc.match))
            declaration = match(DECL_PATTERN, tail)
            isnothing(declaration) && continue
            name = declaration[1]
            name in PUBLIC_NAMES || continue
            push!(entries, (name, String(doc[1]), file))
        end
    end
    entries
end

function reachable(graph, start)
    seen, todo = Set([start]), [start]
    while !isempty(todo)
        for neighbor in setdiff(graph[pop!(todo)], seen)
            push!(seen, neighbor)
            push!(todo, neighbor)
        end
    end
    seen
end

function audit_documentation(; write_graph=false)
    entries = public_docstrings()
    graph = Dict(n => Set{String}() for n in PUBLIC_NAMES)
    example_count = 0
    @testset "Public help and runnable examples" begin
        @test Set(first.(entries)) == PUBLIC_NAMES
        for name in ("matrix_element", "decay_width", "mass_correction_factor")
            @test count(entry -> first(entry) == name, entries) == 1
        end
        for name in sort(collect(PUBLIC_NAMES))
            @test Base.Docs.hasdoc(QMT, Symbol(name))
        end
        for (name, body, file) in entries
            @testset "$name ($(basename(file)))" begin
                examples = collect(eachmatch(EXAMPLE_PATTERN, body))
                @test !isempty(examples)
                @test occursin("## Related", body)
                for link in eachmatch(LINK_PATTERN, body)
                    target = String(link[1])
                    # Public package links or actual documented GIModel producers.
                    @test target in PUBLIC_NAMES ||
                        (isdefined(GIModel, Symbol(target)) && Base.Docs.hasdoc(GIModel, Symbol(target)))
                    target in PUBLIC_NAMES && target != name && push!(graph[name], target)
                end
                for example in examples
                    # No shared globals: a snippet must work after copying it alone.
                    sandbox = Module(gensym(:DocumentationExample))
                    @test begin
                        Base.include_string(sandbox, example[1], "$file:$name")
                        true
                    end
                    example_count += 1
                end
            end
        end
    end
    @testset "Public documentation graph" begin
        for name in PUBLIC_NAMES
            @test !isempty(graph[name])
            @test any(name in targets for targets in values(graph))
            # Strong connectivity rejects isolated clusters as well as isolated nodes.
            @test reachable(graph, name) == PUBLIC_NAMES
        end
    end
    @testset "README examples" begin
        for example in eachmatch(EXAMPLE_PATTERN, read(joinpath(PACKAGE_ROOT, "README.md"), String))
            @test begin
                Base.include_string(Module(gensym(:ReadmeExample)), example[1], "README.md")
                true
            end
        end
    end
    if write_graph
        edges = sort([(a,b) for (a,targets) in graph for b in targets])
        lines = ["# Public documentation graph", "",
            "Generated from all $(length(PUBLIC_NAMES)) supported API names (all exported).",
            "All $(length(entries)) public docstrings have independently executable examples ($example_count blocks).",
            "Every entry reaches every other entry through explicit help links. GIModel links are validated separately.",
            "", "Regenerate and execute the examples from the repository root:", "",
            "```sh", "julia --project=. QuarkModelTransitions/scripts/audit_documentation.jl", "```", "",
            "The graph describes help navigation, not function calls or numerical dependencies.", "",
            "| Entry | Access | Incoming | Outgoing |", "|---|---|---:|---:|"]
        for name in sort(collect(PUBLIC_NAMES))
            access = Base.isexported(QMT, Symbol(name)) ? "exported" : "qualified"
            push!(lines, "| `$name` | $access | $(count(e -> e[2] == name, edges)) | $(length(graph[name])) |")
        end
        append!(lines, ["", "<details>", "<summary>Complete graph ($(length(edges)) links)</summary>", "", "```mermaid", "flowchart LR"])
        append!(lines, ["    $a --> $b" for (a,b) in edges])
        append!(lines, ["```", "", "</details>", ""])
        println(join(lines, '\n'))
    end
    println("Documentation: $(length(PUBLIC_NAMES)) public names, $(length(entries)) docstrings, $example_count executable examples.")
end

abspath(PROGRAM_FILE) == (@__FILE__) && audit_documentation(; write_graph=true)
