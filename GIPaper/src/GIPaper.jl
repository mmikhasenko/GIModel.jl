module GIPaper

# Comparison layer for the Godfrey-Isgur paper reproduction: reference-spectrum
# catalogs (digitized paper data), row -> Meson mapping, model-vs-paper
# comparison, and residual markdown reports. The physics computation lives in
# GIModel; this package only consumes its public API.

using GIModel
using CSV
using Printf
using LinearAlgebra

export paper_data_dir, reference_spectrum_path, model_parameters_path
include("paths.jl")

export ReferenceState, load_reference_spectrum
include("reference_state.jl")

export reference_meson, quark_for
include("reference_meson.jl")

export GI_PSEUDOSCALAR_FIG5_TARGETS_GEV, table_iii_amplitude
include("table_iii_annihilation.jl")

export load_table_vi, load_table_vi_states
include("table_vi.jl")

export load_table_v
include("table_v.jl")

export compare_reference
include("comparison.jl")

export mixing_prone_state, nonmixing_deviation_summary, write_residual_report
include("residual_report.jl")

end
