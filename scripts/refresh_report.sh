#!/usr/bin/env bash
# Run from any directory. Uses the current checkout, including local changes.
# Keep the experimental PDG snapshot fixed; this refresh concerns model outputs.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export GI_REUSE_SPECTRUM=false
PYTHON="${PYTHON:-python3}"

bash scripts/verify_project.sh
GI_TABLE_VI_SOLVER=fd julia GIPaper/scripts/audit_table_vi_photon_decays.jl
julia GIPaper/scripts/compare_table_vi_solvers.jl
julia --project=GIPaper/scripts GIPaper/scripts/investigate_excited_eta_moments.jl
julia --project=GIPaper/scripts GIPaper/scripts/investigate_tensor_photons.jl
julia GIPaper/scripts/plot_ten_meson_spectra.jl
julia GIPaper/scripts/plot_ten_meson_spectra_simplified.jl
for stem in ten_meson_spectra_pdg2026 ten_meson_spectra_simplified_pdg2026; do
    cp "GIPaper/scripts/spectrum_plots/$stem.png" "report/figures/$stem.png"
    cp "GIPaper/scripts/spectrum_plots/$stem.pdf" "report/figures/$stem.pdf"
done
julia report/sources/pwave_residuals.jl
julia --project=examples examples/density_candidates.jl precompute
julia --project=examples report/sources/vector_meson_densities.jl

refresh_census_masses() {
    julia --project=GIPaper/scripts GIPaper/scripts/charmed_threshold_masses.jl "$@"
    julia --project=GIPaper/scripts GIPaper/scripts/five_sector_threshold_masses.jl "$@"
}
refresh_census_masses
refresh_census_masses --refined
"$PYTHON" GIPaper/scripts/filter_charmed_thresholds.py
julia --project=GIPaper/scripts GIPaper/scripts/count_five_sector_thresholds.jl
"$PYTHON" GIPaper/scripts/verify_five_sector_thresholds.py

# This deliberately uses central P/D waves and resummed-contact S waves.
# The all-L contact correction leaves that prescription unchanged. Recheck
# archived traces and rerender; a full new study is a separate explicit command.
"$PYTHON" GIPaper/scripts/verify_a_class_artifacts.py
julia report/sources/a_class_threshold.jl
"$PYTHON" report/scripts/update_computed_text.py
(cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex)
