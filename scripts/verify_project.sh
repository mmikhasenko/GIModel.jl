#!/usr/bin/env bash
# Full project verification: GIModel (pure computation, repo root) and
# GIPaper (paper comparison sub-package, GIPaper/).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

python3 GIPaper/scripts/data_checks.py validate
# FiniteDifferences lives in the test target, so run the suites through Pkg.test.
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
julia GIPaper/scripts/audit_table_iii_mixings.jl
julia GIPaper/scripts/reproduce_table_v.jl
julia GIPaper/scripts/audit_table_vii.jl   # Table VII gluonic annihilation (zero-parameter)
julia GIPaper/scripts/audit_w6_ho_order.jl # W6: paper-order (HO first-order) spin-distorted waves
# The two FD-vs-oscillator audits. Ungated until now, and it showed: the
# GIModel/GIPaper split (4dfd4c5) broke audit_nonmixing_contact.jl outright, and
# nobody noticed for a month because nothing ran it -- so both reports were
# frozen at that commit and silently stale. ~90 s for the pair.
julia GIPaper/scripts/audit_appendix_a_ho_comparison.jl
julia GIPaper/scripts/audit_nonmixing_contact.jl
julia GIPaper/scripts/analyze_heavy_quarkonium.jl
julia GIPaper/scripts/run_all_spectrum_checks.jl
# Scores by READING two reports -- table_iii_mixing_audit.md and
# isoscalar_residuals.md -- so it must run after both are regenerated. It used to
# run first, which made the scorecard show the PREVIOUS run's numbers: a change
# surfaced as report drift one gate later, attributed to whatever was in flight
# then. Keep this line last but one.
python3 GIPaper/scripts/data_checks.py score-annihilation
julia GIPaper/scripts/check_manifest.jl   # anti-drift gate: manifest links must resolve
echo "verify_project.sh: all checks passed."
