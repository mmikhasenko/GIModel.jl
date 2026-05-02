#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf

root = dirname(@__DIR__)
using GIModel

# Parameters + quark table → reference CSV → attach_constituent_masses → compute_sector / compare → markdown reports.

params_path = joinpath(root, "data", "parameters.provisional.toml")
params, mq = load_parameters_and_quark_masses(params_path)

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
    annotated = attach_constituent_masses(mq, reference, mq[flavor])
    computed = compute_sector(params, annotated; kinetic = :relativistic)
    rows = compare(
        computed,
        annotated;
        contact_hyperfine = true,
        use_fine_structure = params.fine_structure,
    )
    write_residual_report(
        report_path,
        title,
        rows;
        kinetic = :relativistic,
        contact_hyperfine = true,
        appendix_a_smearing = params.appendix_a_smearing,
        appendix_a_derivative_g = params.appendix_a_derivative_g,
        appendix_a_closed_form = params.appendix_a_closed_form,
        appendix_a_momentum_sandwich = params.appendix_a_momentum_sandwich,
        contact_momentum_sandwich = params.contact_momentum_sandwich,
        fine_structure_momentum_sandwich = params.fine_structure_momentum_sandwich,
        fine_structure_smeared_kernels = params.fine_structure_smeared_kernels,
        coulomb_1d_smear = params.coulomb_1d_smear,
        use_fine_structure = params.fine_structure,
    )
    absres = abs.([row.residual_MeV for row in rows])
    @printf(
        "%s: rows=%d mean_abs=%.1f MeV max_abs=%.1f MeV -> %s\n",
        name,
        length(rows),
        sum(absres) / length(absres),
        maximum(absres),
        report_path
    )
end
