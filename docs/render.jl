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
using TOML

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

# Landing-page figures (docs/src/public/home/) come from a plain script.
isempty(ARGS) && run(`julia --project=$QUARTO $(joinpath(QUARTO, "home_figures.jl"))`)

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

# Hash of each rendered source, so the site build can detect a .qmd edited
# without re-rendering (checked in docs/make.jl).
const MANIFEST = joinpath(QUARTO, "rendered.toml")
manifest = isfile(MANIFEST) ? TOML.parsefile(MANIFEST) : Dict{String,Any}()

for source in sources
    name = splitext(source)[1]
    rendered = joinpath(OUTPUT, name * ".md")
    target = joinpath(SRC, name * ".md")
    mkpath(dirname(target))
    manifest[source] = bytes2hex(sha256(read(joinpath(QUARTO, source))))
    # Invisible on the page: the "Edit this page" link opens the .qmd source.
    edit_url = relpath(joinpath(QUARTO, source), dirname(target))
    header = "```@meta\nEditURL = \"$edit_url\"\n```\n\n"
    write(target, header * fence_outputs(read(rendered, String)))
    # Figures live next to the page in `<name>_files/`.
    figures = joinpath(OUTPUT, name * "_files")
    destination = joinpath(SRC, name * "_files")
    rm(destination; force = true, recursive = true)
    isdir(figures) && cp(figures, destination)
    println("wrote ", relpath(target, dirname(DOCS)))
end

open(MANIFEST, "w") do io
    println(io, "# Written by docs/render.jl: SHA-256 of each .qmd when it was last rendered.")
    TOML.print(io, manifest; sorted = true)
end
