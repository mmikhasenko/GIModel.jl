# Clean data (Phase 2 target)

`data/clean/` holds **audited, promoted** tables once raw extraction is verified.
Nothing here supersedes the 1985 paper: each row should trace to raw lines and
page-level evidence (see `docs/orchestrator_task.md` Phase 2).

Planned files from the plan:

- `masses.csv` — normalized masses with `source`, `page`, and confidence.
- `mixings.csv` — isospin/strange—nonstrange (Table III) when promoted.
- `parameters.toml` — long-term: single canonical TOML read by the solver (today
  the project uses `data/parameters.provisional.toml` at the top level for the
  diagnostic build).

The repository currently uses top-level `data/reference_spectrum_*.csv` and
`data/parameters.provisional.toml` for the working solver. Migration into
`data/clean/` is optional bookkeeping until a full extraction audit is done.
