#!/usr/bin/env bash
# Full project verification.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

python3 scripts/data_checks.py validate
python3 scripts/data_checks.py score-annihilation
# `test/runtests.jl` expects the GIModel environment.
julia --project=. test/runtests.jl
julia --project=. scripts/analyze_heavy_quarkonium.jl
julia --project=. scripts/run_all_spectrum_checks.jl
echo "verify_project.sh: all checks passed."
