#!/usr/bin/env bash
# Full project verification: GIModel (including transitions) and GIPaper.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# FiniteDifferences lives in the test target, so run the suites through Pkg.test.
# GIPaper's Julia test suite owns the CSV, TOML, and provenance invariants.
GI_HEAVY_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
julia GIPaper/scripts/audit_ho_convergence.jl # PA-12: adaptive native-HO certificate
julia GIPaper/scripts/audit_fd_convergence.jl # FD-COMP: independent grid/domain certificate
julia GIPaper/scripts/audit_table_iii_mixings.jl
julia --project=GIPaper/scripts GIPaper/scripts/trace_rate_inputs.jl # regenerates Tables V/VI/VII and input traces
julia --project=GIPaper/scripts GIPaper/scripts/check_input_traces.jl
julia GIPaper/scripts/audit_w6_ho_order.jl # W6: paper-order HO full diagonalization
# The two FD-vs-oscillator audits. Ungated until now, and it showed: the
# GIModel/GIPaper split (4dfd4c5) broke audit_nonmixing_contact.jl outright, and
# nobody noticed for a month because nothing ran it -- so both reports were
# frozen at that commit and silently stale. ~90 s for the pair.
julia GIPaper/scripts/audit_appendix_a_ho_comparison.jl
julia GIPaper/scripts/audit_nonmixing_contact.jl
# The rest of the report writers, likewise ungated until C2 went looking. The
# Table VI run found the last place A17 had not reached: the hindered
# Upsilon' -> eta_b gamma amplitude, a near-cancelling overlap of distorted
# bottomonium S waves, moved 0.0084 -> 0.0083 (paper +0.007, verdict unchanged).
# mixing_studies.jl also rewrites two PNGs; they are byte-reproducible run to run.
julia GIPaper/scripts/audit_mixing_angles.jl
julia GIPaper/scripts/audit_realistic_factors.jl
julia --project=GIPaper/scripts GIPaper/docs/mixing_studies.jl
julia GIPaper/scripts/analyze_heavy_quarkonium.jl
julia GIPaper/scripts/run_all_spectrum_checks.jl
# Scores by READING two reports -- table_iii_mixing_audit.md and
# isoscalar_residuals.md -- so it must run after both are regenerated. It used to
# run first, which made the scorecard show the PREVIOUS run's numbers: a change
# surfaced as report drift one gate later, attributed to whatever was in flight
# then. Keep it after every report writer.
julia GIPaper/scripts/score_annihilation.jl
# Paper-layout Tables I/II/IV and Figs. 1-9 read the reports regenerated above.
julia GIPaper/scripts/paper_tables/generate.jl
julia GIPaper/scripts/paper_tables/check.jl
julia GIPaper/scripts/check_manifest.jl   # anti-drift gate: manifest links must resolve
echo "verify_project.sh: all checks passed."
