#!/usr/bin/env julia
# Recompute residual reports for all `data/reference_spectrum_*.csv` files.

using Printf

root = dirname(@__DIR__)
include(joinpath(root, "src", "GIModel", "GIModel.jl"))
using .GIModel

params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
data_dir = joinpath(root, "data")
report_dir = joinpath(root, "docs", "residual_reports")
mkpath(report_dir)

const FLAVOR = Dict(
    "reference_spectrum_charmonium.csv" => "c",
    "reference_spectrum_bottomonium.csv" => "b",
    "reference_spectrum_charmed.csv" => "c",
    "reference_spectrum_b_flavored.csv" => "b",
    "reference_spectrum_strange.csv" => "q",
    "reference_spectrum_isovector.csv" => "q",
    "reference_spectrum_isoscalar.csv" => "q",
)

for fn in readdir(data_dir)
    startswith(fn, "reference_spectrum_") && endswith(fn, ".csv") || continue
    flavor = get(FLAVOR, fn, "q")
    ref_path = joinpath(data_dir, fn)
    base = fn[length("reference_spectrum_")+1:end-length(".csv")]
    report_path = joinpath(report_dir, string(base, "_residuals.md"))
    title = @sprintf("Residuals: %s (GI-style)", base)
    reference = load_reference_spectrum(ref_path)
    rows = compare_sector(
        params, reference, flavor;
        kinetic = :relativistic,
        contact_hyperfine = true,
        use_fine_structure = params.fine_structure,
    )
    write_residual_report(
        report_path, title, rows;
        kinetic = :relativistic,
        contact_hyperfine = true,
        appendix_a_central = params.appendix_a_central,
        use_fine_structure = params.fine_structure,
    )
    absres = [abs(r.residual_MeV) for r in rows]
    isempty(absres) && @printf("%s: no rows\n", fn)
    !isempty(absres) &&
        @printf("%s: n=%d mean_abs=%.1f max_abs=%.1f\n", fn, length(rows), sum(absres) / length(absres), maximum(absres))
end
