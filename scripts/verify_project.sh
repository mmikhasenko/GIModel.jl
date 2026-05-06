#!/usr/bin/env bash
# Full project verification (same sequence as `scripts/autonomous_loop.py` DEFAULT_VERIFY).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# Julia is installed via `juliaup` on many systems, and the launcher creates a lockfile
# next to its configuration. In sandboxed environments (like Codex), `$HOME` may be
# read-only, causing `julia` invocations to fail unless `JULIAUP_DEPOT_PATH` points to
# a writable location.
export JULIAUP_DEPOT_PATH="${JULIAUP_DEPOT_PATH:-${TMPDIR:-/tmp}/juliaup_depot}"
mkdir -p "$JULIAUP_DEPOT_PATH"

python3 scripts/validate_seed.py
python3 scripts/verify_table_ii_toml.py
python3 scripts/validate_reference_spectra.py
python3 scripts/validate_clean_data.py
# `test/runtests.jl` expects the GIModel environment (unlike driver scripts that call `Pkg.activate`).
julia --project="$ROOT" test/runtests.jl
julia scripts/compare_central_pointwise_vs_a7a8.jl
julia scripts/analyze_heavy_quarkonium.jl
julia scripts/run_all_spectrum_checks.jl
echo "verify_project.sh: all checks passed."
