# Shared setup for the paper-table drivers. Every driver writes one object of
# Godfrey & Isgur, Phys. Rev. D 32, 189 (1985), in the paper's own row/column
# layout, under GIPaper/docs/paper_tables/.
if !isdefined(@__MODULE__, :PAPER_TABLES_DIR)
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

    using Printf
    using GIModel
    using GIPaper

    const GIPAPER_DIR = pkgdir(GIPaper)
    const PAPER_TABLES_DIR = joinpath(GIPAPER_DIR, "docs", "paper_tables")
    const RESIDUAL_REPORTS_DIR = joinpath(GIPAPER_DIR, "docs", "residual_reports")

    function write_paper_table(name::AbstractString, body::AbstractString)
        mkpath(PAPER_TABLES_DIR)
        path = joinpath(PAPER_TABLES_DIR, name)
        open(path, "w") do io
            write(io, body)
            endswith(body, '\n') || write(io, '\n')
        end
        return path
    end
end
