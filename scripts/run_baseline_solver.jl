#!/usr/bin/env julia

using Printf

root = dirname(@__DIR__)
include(joinpath(root, "src", "GIModel", "GIModel.jl"))
using .GIModel

params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))

sectors = [
    (
        "ccbar",
        "c",
        joinpath(root, "data", "reference_spectrum_charmonium.csv"),
        joinpath(root, "docs", "residual_reports", "ccbar_baseline.md"),
        "Charmonium Baseline Residuals",
    ),
    (
        "bbbar",
        "b",
        joinpath(root, "data", "reference_spectrum_bottomonium.csv"),
        joinpath(root, "docs", "residual_reports", "bbbar_baseline.md"),
        "Bottomonium Baseline Residuals",
    ),
]

for (name, flavor, reference_path, report_path, title) in sectors
    reference = load_reference_spectrum(reference_path)
    rows = compare_sector(params, reference, flavor; kinetic = :relativistic, contact_hyperfine = true)
    write_residual_report(report_path, title, rows; kinetic = :relativistic, contact_hyperfine = true)
    absres = abs.([row.residual_MeV for row in rows])
    @printf("%s: rows=%d mean_abs=%.1f MeV max_abs=%.1f MeV -> %s\n", name, length(rows), sum(absres) / length(absres), maximum(absres), report_path)
end
