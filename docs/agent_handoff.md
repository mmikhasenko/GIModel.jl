# Agent Handoff

This repository is ready for a new agent or a small set of specialized agents
to continue the Godfrey-Isgur reproduction.

## Current Status

- Git repository is initialized on `main`.
- Phase 0 scaffold is in place.
- Primary paper is available at `paper/Godfrey-Isgur-1985.pdf`.
- Preferred prose reference is `paper/text/godfrey_isgur_1985_prose.md`.
- Seed data is available under `data/seed/`.
- Seed schema validation passes with `python3 scripts/validate_seed.py`.
- Spectrum pages 6-10 are rendered under `paper/screenshots/spectrum_pages/`.
- First-pass table digitizations live under `data/raw/digitized_tables/`.
- First-pass figure-label digitizations live under
  `data/raw/digitized_figures/`.
- A diagnostic Godfrey-Isgur solver lives under `src/` as Julia package
  **GIModel** (`Project.toml`; module entry `src/GIModel.jl` plus included
  sources beside it): semirelativistic kinetic + central + running Coulomb,
  contact hyperfine, first-order fine-structure, heavy-quarkonium comparisons.
  Provisional parameters are in `data/parameters.provisional.toml`. Formula
  intent and flags are described in `docs/formula_map.md`; sector residuals under
  `docs/residual_reports/`.
- `test/runtests.jl` encodes several convention checks (Coulomb derivative,
  fine-structure angular factors, reduced-radial expectations, smearing
  constant-preservation).
- No cleaned *promoted* reference dataset in `data/clean/` has been finished
  yet; extraction remains partly first-pass and needs audit.
- **Autonomous coding loop:** the guarded outer loop and program file are
  `docs/autonomous_loop.md` and `docs/autonomous_program.md` (runner:
  `scripts/autonomous_loop.py`). Logs are written to `docs/autonomous_runs/`
  (gitignored). If the agent CLI hits a usage limit, the runner stops the batch
  early; resume after credits reset or by doing a manual bounded iteration
  (same program, same full verification gate as in the docs).

## Start Here

1. Read `README.md`.
2. Read `docs/orchestrator_task.md`.
3. Read `docs/paper_navigation.md`.
4. Use `paper/text/godfrey_isgur_1985_prose.md` for prose search.
5. Use `paper/Godfrey-Isgur-1985.pdf` as the authority for every equation,
   symbol, table value, and state label.

## Recommended Subagent Split

### Orchestrator

Owns phase gates, status updates, and final integration.

Start files:

- `docs/orchestrator_task.md`
- `docs/agent_handoff.md`
- `docs/extraction_notes.md`

First task:

- Convert the next milestone into issues or checklist items and enforce that
  extracted values keep provenance.

### Extractor

Owns raw extraction from the original PDF.

Start files:

- `paper/Godfrey-Isgur-1985.pdf`
- `paper/text/pdftotext_layout.txt`
- `paper/text/pdftotext_bbox.html`
- `paper/text/pdftotext_words.tsv`
- `docs/paper_navigation.md`

First task:

- Render page images, identify all numerical tables/figure spectra, and begin
  `data/raw/tables_from_pdf_text.csv`.

### Physics Reviewer

Owns conventions, formula interpretation, and state labels.

Start files:

- `paper/Godfrey-Isgur-1985.pdf`
- `paper/text/godfrey_isgur_1985_prose.md`
- `docs/conventions.md`
- `docs/formula_map.md`
- `docs/paper_navigation.md`

First task:

- Draft conventions for sectors, spectroscopic notation, unequal-mass mixing,
  and which paper pages support each convention.

### Data Cleaner

Owns promotion from raw extraction to clean physics tables.

Start files:

- `data/seed/godfrey_isgur_seed_masses.csv`
- `data/seed/godfrey_isgur_sources.csv`
- `data/raw/`
- `docs/conventions.md`

First task:

- Wait for raw extracted rows, then define `data/clean/masses.csv`,
  `data/clean/mixings.csv`, and provenance links back to raw rows.

### Numerics Implementer

Owns solver code after parameter/formula extraction is stable.

Start files:

- `docs/formula_map.md`
- `docs/conventions.md`
- `data/parameters.provisional.toml` (provisional Table II–style input)
- `src/`

First task:

- Move the diagnostic implementation toward a paper-faithful Appendix A
  central/smeared operator and keep `docs/formula_map.md` synchronized; avoid
  silent refits and prefer tests plus residual classification over ad hoc tuning
  (see `docs/autonomous_program.md`).

### Verifier

Owns checks, residual reports, and discrepancy classification.

Start files:

- `scripts/validate_seed.py`
- `test/`
- `docs/residual_reports/`

First task:

- Add extraction/data validation tests before validating model numerics.

## Immediate Next Milestone

Complete Phase 1 enough to support the first model milestone:

- Render PDF page images.
- Extract Table II parameters.
- Extract `ccbar` and `bbbar` spectra from the original paper.
- Create raw extraction CSVs with page/table provenance.
- Log disagreements in `data/raw/extraction_audit.csv` (header row is in
  place; add one row per audited mismatch between extractors).

The first two bullets have a provisional start: Table II is digitized under
`data/raw/digitized_tables/table_ii_parameters/`, and spectrum screenshots are
available for visual extraction. They are not yet audited clean data.
Figures 3-9 are also digitized as first-pass label tables under
`data/raw/digitized_figures/`. Figure 3 has its original single-figure plot
script; Figures 4-9 are regenerated with
`python3 scripts/plot_spectrum_digitizations.py`.

## Do Not Trust By Default

- Later quoted GI values in `data/seed/`.
- Markdown extraction for numerical values.
- Intermediate Markdown files in `paper/text/` other than
  `godfrey_isgur_1985_prose.md`.
- Any sector marked complete without residuals and discrepancy classification.
