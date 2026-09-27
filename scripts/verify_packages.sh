#!/usr/bin/env bash
# Fast GIModel/GIPaper gate; full convergence and reproduction gate: verify_project.sh.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
