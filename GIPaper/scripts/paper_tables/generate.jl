#!/usr/bin/env julia
# Regenerate every paper-layout object under GIPaper/docs/paper_tables/.
#
#   julia GIPaper/scripts/paper_tables/generate.jl
#
# Tables III, V, VI, VII and the mixing-angle captions are the canonical
# residual reports (docs/residual_reports/); Table IV and the column registry
# read them, so regenerate those reports first (scripts/verify_project.sh does).
include(joinpath(@__DIR__, "common.jl"))
for driver in ("mass_inputs", "table_i", "table_ii", "table_iv", "fig_i", "fig_ii",
               "spectra", "caption_mixing", "quantity_columns")
    include(joinpath(@__DIR__, "$driver.jl"))
end

artifacts = String[
    compute_mass_inputs(), compute_table_i(), compute_table_ii(), compute_table_iv(),
    compute_fig_i(), compute_fig_ii(), values(compute_spectra())...,
    compute_caption_mixing(), compute_quantity_columns(),
]
println("paper tables: wrote $(length(artifacts)) objects")
foreach(path -> println("  " * relpath(path, GIPAPER_DIR)), sort(artifacts))
