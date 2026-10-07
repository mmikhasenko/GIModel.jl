#!/usr/bin/env bash
# Model tests and direct paper comparisons.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
GI_HEAVY_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
for check in run_all_spectrum_checks audit_table_iii_mixings audit_mixing_angles reproduce_table_v audit_table_vi_photon_decays audit_heavy_light_e1 audit_table_vii; do
    julia "GIPaper/checks/$check.jl"
done
echo "verify_project.sh: all checks passed."
