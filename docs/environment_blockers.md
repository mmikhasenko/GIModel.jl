# Environment Blockers

This file tracks environment issues that prevent running the repository’s
standard verification gates in sandboxed/CI-like contexts.

## Julia runtime missing (juliaup installed, no channel configured)

- Symptom: `julia --version` (and any `julia --project=...`) fails with a
  juliaup/launcher error about not being able to create a lockfile or not being
  able to determine a juliaup channel.
- Root cause in this sandbox: only `juliaup` is present; no Julia toolchain
  (channel) is installed/configured, and the sandbox disallows writing to the
  default home directory locations juliaup tries to use.
- Impact: cannot run:
  - `julia --project=. test/runtests.jl`
  - `julia --project=. scripts/analyze_heavy_quarkonium.jl`
  - `julia --project=. scripts/run_all_spectrum_checks.jl`

### Suggested remediation (outside strict sandboxes)

- Install a Julia channel via juliaup (example): `juliaup add release` and set a
  default channel (example): `juliaup default release`.
- If the runtime must run in restricted environments, ensure the launcher and
  depot/config paths are writable (for example by setting `JULIAUP_DEPOT_PATH`
  and `JULIA_DEPOT_PATH` to a writable directory).

