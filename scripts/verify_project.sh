#!/usr/bin/env bash
# Full project verification (same spirit as `scripts/autonomous_loop.py` gate).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
python3 scripts/validate_seed.py
python3 scripts/verify_table_ii_toml.py
python3 scripts/validate_reference_spectra.py
julia --project=. test/runtests.jl
julia --project=. scripts/compare_central_pointwise_vs_a7a8.jl
julia --project=. scripts/analyze_heavy_quarkonium.jl
julia --project=. scripts/run_all_spectrum_checks.jl
echo "verify_project.sh: all checks passed."
