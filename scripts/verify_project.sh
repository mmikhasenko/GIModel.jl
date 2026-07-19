#!/usr/bin/env bash
# Full project verification: GIModel (pure computation, repo root) and
# GIPaper (paper comparison sub-package, GIPaper/).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

python3 GIPaper/scripts/data_checks.py validate
python3 GIPaper/scripts/data_checks.py score-annihilation
# FiniteDifferences lives in the test target, so run the suites through Pkg.test.
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
julia GIPaper/scripts/audit_table_iii_mixings.jl
julia GIPaper/scripts/reproduce_table_v.jl
julia GIPaper/scripts/audit_table_vii.jl   # Table VII gluonic annihilation (zero-parameter)
julia GIPaper/scripts/analyze_heavy_quarkonium.jl
julia GIPaper/scripts/run_all_spectrum_checks.jl
julia GIPaper/scripts/check_manifest.jl   # anti-drift gate: manifest links must resolve
echo "verify_project.sh: all checks passed."
