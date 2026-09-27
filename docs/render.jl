#!/usr/bin/env julia
# Execute the documentation pages that run code and write them as Markdown
# into docs/src/, where Documenter picks them up. Run from the repository root
# before a release, or whenever the physics or the examples change:
#
#   julia docs/render.jl            # every page in docs/quarto/
#   julia docs/render.jl tutorials/charmonium.qmd   # selected pages
#
# Requires Quarto (https://quarto.org) and Julia; the notebook environment is
# docs/quarto/Project.toml. The generated files are committed, so the CI
# documentation build only assembles the site and runs no calculations.

using SHA

const DOCS = @__DIR__
const QUARTO = joinpath(DOCS, "quarto")
const OUTPUT = joinpath(QUARTO, "_output")
const SRC = joinpath(DOCS, "src")

run(`julia --project=$QUARTO -e "using Pkg; Pkg.instantiate()"`)

sources = isempty(ARGS) ?
    [relpath(joinpath(root, f), QUARTO) for (root, _, files) in walkdir(QUARTO)
     for f in files if endswith(f, ".qmd") && !occursin("_output", root)] :
    ARGS
# `quarto render` takes a single input: render the whole project, or each
# requested page in turn.
cd(QUARTO) do
    if isempty(ARGS)
        run(`quarto render`)
    else
        foreach(source -> run(`quarto render $source`), sources)
    end
end

# Quarto writes cell outputs as indented code blocks; turn them into fenced
# blocks so Documenter and VitePress render them consistently.
function fence_outputs(text)
    lines = split(text, '\n')
    out = String[]
    i = 1
    in_fence = false
    while i <= length(lines)
        line = lines[i]
        if startswith(line, "```")
            in_fence = !in_fence
            push!(out, replace(line, r"^``` (\w)" => s"```\1"))
            i += 1
        elseif !in_fence && startswith(line, "    ") && (i == 1 || isempty(strip(lines[i-1])))
            block = String[]
            while i <= length(lines) && (startswith(lines[i], "    ") ||
                    (isempty(strip(lines[i])) && i < length(lines) && startswith(lines[i+1], "    ")))
                push!(block, isempty(strip(lines[i])) ? "" : lines[i][5:end])
                i += 1
            end
            append!(out, ["```", block..., "```"])
        else
            push!(out, line)
            i += 1
        end
    end
    text = join(out, '\n')
    # Quarto emits figures as HTML <img> tags; Documenter handles Markdown images.
    return replace(text, r"<img src=\"([^\"]+)\"[^>]*/>" => s"![](\1)")
end

for source in sources
    name = splitext(source)[1]
    rendered = joinpath(OUTPUT, name * ".md")
    target = joinpath(SRC, name * ".md")
    mkpath(dirname(target))
    # The source hash lets the site build detect a .qmd edited without re-rendering.
    digest = bytes2hex(sha256(read(joinpath(QUARTO, source))))
    header = "<!-- Generated from docs/quarto/$source by docs/render.jl. Edit the .qmd file. source-sha256: $digest -->\n\n"
    write(target, header * fence_outputs(read(rendered, String)))
    # Figures live next to the page in `<name>_files/`.
    figures = joinpath(OUTPUT, name * "_files")
    destination = joinpath(SRC, name * "_files")
    rm(destination; force = true, recursive = true)
    isdir(figures) && cp(figures, destination)
    println("wrote ", relpath(target, dirname(DOCS)))
end
