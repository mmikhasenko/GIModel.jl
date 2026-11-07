#!/usr/bin/env bash
# Run from any directory. JULIA/PYTHON may select preconfigured executables.
set -euo pipefail
cd "$(dirname "$0")/../.."
JULIA="${JULIA:-julia}"
PYTHON="${PYTHON:-python3}"
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/gi-weak-mpl}"
mkdir -p WeakDecays/results
"$JULIA" --project=. WeakDecays/test/runtests.jl > WeakDecays/results/tests.log 2>&1
for method in ho32 ho48 fd600; do
  "$JULIA" --project=. WeakDecays/scripts/solve_waves.jl "$method" > "WeakDecays/results/solve_$method.log" 2>&1
  "$JULIA" --project=. WeakDecays/scripts/benchmark.jl "$method" > "WeakDecays/results/benchmark_$method.log" 2>&1
done
for experiment in lightfront heavy_limit rescattering strong_vertex extended_checks; do
  "$JULIA" --project=. "WeakDecays/scripts/$experiment.jl" > "WeakDecays/results/$experiment.log" 2>&1
done
"$PYTHON" WeakDecays/scripts/analyze.py
"$PYTHON" WeakDecays/scripts/make_tables.py
"$PYTHON" WeakDecays/scripts/validate_outputs.py
