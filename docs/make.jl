# Assemble the documentation site from the repository root:
#
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#   julia --project=docs docs/make.jl
#
# This step runs no calculations. Pages with executed examples are written in
# docs/quarto/*.qmd and rendered to docs/src/ by `julia docs/render.jl`, whose
# output is committed; re-render before a release. Documenter resolves the
# cross-references and API docstrings, and VitePress builds the static site
# (it installs its own Node.js). For a live preview after a build, run from
# docs/:
#
#   npm run docs:dev

using Documenter
using SHA
using DocumenterVitepress
using GIModel
using GIModel.QuarkModelTransitions

const DOCS = @__DIR__
const ROOT = dirname(DOCS)

# GIPaper writes its comparison reports as Markdown. They are copied into the
# site as they are, so the published numbers are exactly the recorded ones.
const REPORTS = joinpath(ROOT, "GIPaper", "reports")
const REPORT_DEST = joinpath(DOCS, "src", "paper", "reports")
const REPORT_PAGES = [
    "Spectrum scorecard" => "scorecard.md",
    "Mixing angles" => "mixing_angles.md",
    "Table III isoscalar mixing" => "table_iii_mixing_audit.md",
    "Table V strong decays" => "table_v_reproduction.md",
    "Table VI photon decays" => "table_vi_photon_decays.md",
    "Table VII annihilation" => "table_vii_annihilation_em.md",
]
const SECTOR_REPORTS = [
    "isovector_residuals.md", "isoscalar_residuals.md", "strange_residuals.md",
    "charmonium_residuals.md", "charmed_residuals.md",
    "bottomonium_residuals.md", "b_flavored_residuals.md",
]

function copy_reports()

# Executed pages are rendered locally (docs/render.jl) and committed. Refuse to
# publish a page whose .qmd source changed after it was rendered.
function check_rendered_pages()
    quarto = joinpath(DOCS, "quarto")
    stale = String[]
    for (root, _, files) in walkdir(quarto), file in files
        endswith(file, ".qmd") && !occursin("_output", root) || continue
        source = relpath(joinpath(root, file), quarto)
        page = joinpath(DOCS, "src", splitext(source)[1] * ".md")
        digest = bytes2hex(sha256(read(joinpath(root, file))))
        isfile(page) && occursin("source-sha256: $digest", read(page, String)) ||
            push!(stale, source)
    end
    isempty(stale) || error("rendered pages are out of date: $(join(stale, ", ")); " *
                            "run `julia docs/render.jl` and commit docs/src")
end

check_rendered_pages()
    rm(REPORT_DEST; force = true, recursive = true)
    mkpath(REPORT_DEST)
    for file in vcat(last.(REPORT_PAGES), SECTOR_REPORTS)
        text = read(joinpath(REPORTS, file), String)
        # Links between copied reports stay relative; every other relative
        # link (data files, CSV companions) points at the file on GitHub.
        copied = Set(vcat(last.(REPORT_PAGES), SECTOR_REPORTS))
        text = replace(text, r"\]\((?!https?://|#)([^)\s]+)\)" => function (link)
            target = match(r"\]\((.*)\)", link)[1]
            target in copied && return link
            path = normpath(joinpath("GIPaper", "reports", target))
            return "](https://github.com/mmikhasenko/GIModel.jl/blob/main/$path)"
        end)
        write(joinpath(REPORT_DEST, file), text)
    end
end

copy_reports()

# Executed pages are rendered locally (docs/render.jl) and committed. Refuse to
# publish a page whose .qmd source changed after it was rendered.
function check_rendered_pages()
    quarto = joinpath(DOCS, "quarto")
    stale = String[]
    for (root, _, files) in walkdir(quarto), file in files
        endswith(file, ".qmd") && !occursin("_output", root) || continue
        source = relpath(joinpath(root, file), quarto)
        page = joinpath(DOCS, "src", splitext(source)[1] * ".md")
        digest = bytes2hex(sha256(read(joinpath(root, file))))
        isfile(page) && occursin("source-sha256: $digest", read(page, String)) ||
            push!(stale, source)
    end
    isempty(stale) || error("rendered pages are out of date: $(join(stale, ", ")); " *
                            "run `julia docs/render.jl` and commit docs/src")
end

check_rendered_pages()

DocMeta.setdocmeta!(GIModel, :DocTestSetup, :(using GIModel); recursive = true)

makedocs(;
    root = DOCS,
    sitename = "GIModel.jl",
    authors = "Mikhail Mikhasenko",
    modules = [GIModel, GIModel.QuarkModelTransitions],
    repo = Remotes.GitHub("mmikhasenko", "GIModel.jl"),
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/mmikhasenko/GIModel.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "mmikhasenko.github.io/GIModel.jl",
        # `DOCS_MD_ONLY=true` stops after the Markdown stage (no Node.js needed).
        build_vitepress = get(ENV, "DOCS_MD_ONLY", "false") != "true",
    ),
    checkdocs = :exports,
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Manual" => [
            "The model" => "manual/physics.md",
            "Inputs and parameters" => "manual/inputs.md",
            "Computing a spectrum" => "manual/spectra.md",
            "Solvers and convergence" => "manual/solvers.md",
            "Wavefunctions" => "manual/wavefunctions.md",
            "Isoscalar flavor mixing" => "manual/isoscalar.md",
            "Transitions and decays" => "manual/transitions.md",
            "Conventions and units" => "manual/conventions.md",
        ],
        "Tutorials" => [
            "Charmonium, start to finish" => "tutorials/charmonium.md",
            "Heavy-light mesons and mixing" => "tutorials/heavy_light.md",
            "From charm to bottom" => "tutorials/heavy_quark_sweep.md",
            "η and η′" => "tutorials/eta_etaprime.md",
            "Strong decays beyond the paper" => "tutorials/strong_decays.md",
        ],
        "The 1985 paper" => [
            "Reproducing Godfrey–Isgur" => "paper/overview.md",
            "Results at a glance" => "paper/results.md",
            "Mixing angles: a paper erratum" => "paper/mixing_angles.md",
            "Equations to code" => "paper/formula_map.md",
            "Reading the paper" => "paper/navigation.md",
            "Scope and limitations" => "paper/scope.md",
            "Recorded reports" => [
                (title => "paper/reports/$file" for (title, file) in REPORT_PAGES)...,
                "Residuals by sector" =>
                    ["paper/reports/$file" for file in SECTOR_REPORTS],
            ],
        ],
        "API reference" => [
            "GIModel" => "api/gimodel.md",
            "QuarkModelTransitions" => "api/transitions.md",
        ],
        "Developer notes" => [
            "Architecture" => "developer/architecture.md",
            "Testing and CI" => "developer/testing.md",
            "Writing docstrings" => "developer/docstrings.md",
        ],
        "Learning track" => "course.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/mmikhasenko/GIModel.jl",
    target = joinpath(DOCS, "build"),
    devbranch = "main",
    push_preview = true,
)
