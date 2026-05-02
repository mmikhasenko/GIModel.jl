# Environment Notes

This file tracks environment issues observed while running nested autonomous
agents. These are not necessarily project blockers: the outer autonomous runner
still accepts an iteration only after the standard project gate passes.

## Julia launcher in nested Codex sandboxes

- Symptom inside some nested Codex CLI executions: `julia --version` or
  `julia --project=...` can fail through the `juliaup` launcher with a lockfile
  or channel-selection error.
- Root cause: the nested sandbox may not be able to write to the configured
  juliaup depot/config path. A Julia runtime is installed locally; this is a
  launcher/depot access issue, not evidence that the project lacks Julia.
- Impact for nested agents: direct `julia` commands may fail unless the agent
  uses a writable depot or the resolved Julia binary.
- Impact for accepted project progress: the outer loop has successfully run:
  - `julia --project=. test/runtests.jl`
  - `julia scripts/analyze_heavy_quarkonium.jl`
  - `julia scripts/run_all_spectrum_checks.jl`

### Suggested remediation

- Prefer the standard verification commands from a normal shell (from the repo
  root). Driver scripts call `Pkg.activate` on the repo; **`test/runtests.jl` does not**—run it as
  `julia --project=. test/runtests.jl` (as in `scripts/verify_project.sh`) or via `julia --project=. -e 'using Pkg; Pkg.test()'`.
- If running inside a restricted nested agent sandbox, set writable
  `JULIAUP_DEPOT_PATH`/`JULIA_DEPOT_PATH`, or invoke the resolved Julia binary
  from the local juliaup installation.
- Do not weaken the repository verification gate because of this nested-tooling
  issue.
