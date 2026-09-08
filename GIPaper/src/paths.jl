# Canonical locations of the paper reference data and the model parameter TOML.
#
# Public API (exported from GIPaper.jl): paper_data_dir, reference_spectrum_path,
#   model_parameters_path

"""Directory holding the digitized paper data (reference CSVs, seed/raw/clean extractions)."""
paper_data_dir() = normpath(joinpath(@__DIR__, "..", "data"))

"""
    reference_spectrum_path(sector) -> String

Path of `reference_spectrum_<sector>.csv` under [`paper_data_dir`](@ref), e.g.
`reference_spectrum_path("charmonium")`. Throws `ArgumentError` when absent.
"""
function reference_spectrum_path(sector::AbstractString)
    path = joinpath(paper_data_dir(), "reference_spectrum_$(sector).csv")
    isfile(path) || throw(ArgumentError(
        "no reference spectrum for sector `$sector` (looked at $path)",
    ))
    return path
end

"""Path to the parameter file supplied by the installed GIModel package."""
model_parameters_path() = GIModel.default_parameters_path()
