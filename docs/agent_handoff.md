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
  Inputs split **`GIParameters`** (switches and potential constants) from **`QuarkMassTable`**
  (`[masses]` in TOML); scripts call **`load_parameters_and_quark_masses`**, then
  **`attach_constituent_masses`** → **`Vector{ReferenceStateWithMasses}`** before
  **`compute_sector`** / **`compare`** (**`src/sector_comparison.jl`**; cached radial solves per
  **`RadialChannelKey`** = **`ConstituentMasses`** + `L`). See **`docs/code_architecture.md`**.
  Provisional parameters are in **`data/parameters.provisional.toml`**. Formula intent and flags are
  in **`docs/formula_map.md`**; sector residuals under **`docs/residual_reports/`**.
- `test/runtests.jl` encodes several convention checks (Coulomb derivative,
  fine-structure angular factors, reduced-radial expectations, smearing
  constant-preservation).
- Basis audit checkpoint: `GIParameters{FiniteDifferenceBasis}` remains the
  default solver path, and a first finite `GIParameters{HarmonicOscillatorBasis}`
  path now exists for paper-style basis comparison. The HO path assembles the
  spin-independent Hamiltonian in an oscillator subspace and reconstructs
  mesh wavefunctions for shared diagnostics. Run
  `julia --project=. scripts/audit_nonmixing_contact.jl` to regenerate
  `docs/residual_reports/nonmixing_scorecard.md`,
  `docs/residual_reports/contact_hyperfine_audit.md`, and
  `docs/residual_reports/basis_nonmixing_comparison.md`. Current read: FD and
  HO non-mixed residuals are highly consistent outside a few higher S-wave
  truncation-sensitive rows; the light-sector contact-hyperfine mismatch
  survives the basis comparison and remains the next audit target.
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
4. Read `docs/code_architecture.md` before changing `src/` or residual scripts (parameters vs masses vs sector rows).
5. Use `paper/text/godfrey_isgur_1985_prose.md` for prose search.
6. Use `paper/Godfrey-Isgur-1985.pdf` as the authority for every equation,
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

- `docs/code_architecture.md`
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

## 2026-05-05 Hyperfine Checkpoint

The contact-hyperfine mismatch audit found a concrete Appendix-A convention bug:
the spin-dependent momentum sandwich side exponent was implemented as
`1/4 + epsilon_i/2`, but the paper's energy-denominator limit requires
`1/2 + epsilon_i` on each side. With `epsilon_i = 0`, the two-sided sandwich now
turns the nonrelativistic `1/(m1*m2)` strength into `1/(E1*E2)`.

Implemented in `src/contact_hyperfine.jl` via
`gi_spin_dependent_side_exponent`, reused by `src/spin_fine_structure.jl`, and
guarded by a regression test in `test/runtests.jl`. Regenerated:

- `docs/residual_reports/nonmixing_scorecard.md`
- `docs/residual_reports/contact_hyperfine_audit.md`
- `docs/residual_reports/basis_nonmixing_comparison.md`

Result: open-flavor and light-sector non-mixed residuals improved sharply
without refitting parameters. The remaining stress is mostly the light
pseudoscalar ground states (`isovector 1^1S_0`, `strange 1^1S_0`), so the next
audit should focus on chiral/annihilation limitations and any remaining
contact-kernel ordering details rather than central-potential tuning.

## 2026-05-05 Contact Ordering Checkpoint

The next contact audit found that the remaining light-pseudoscalar edge was not
best addressed by retuning smearing or the published `epsilon_c`. The stronger
paper clue is ordering: GI's first diagonalization is already in fixed
`L,S,J` sectors, so S-wave contact hyperfine participates in the radial
eigenproblem before the later mixing/annihilation stages.

Implemented a finite-difference S-wave nonperturbative contact path in
`compare`: the active `contact_shift_GeV` is now the difference between the
central S-wave level and the level from diagonalizing `H_central + H_contact`
for the requested multiplicity. The perturbative expectation remains available
for non-S waves, non-FD basis diagnostics, and the legacy diagonal-contact
comparison.

Result after regenerating `scripts/audit_nonmixing_contact.jl`: isovector max
non-mixed residual fell from `133.4 MeV` to `55.0 MeV`, strange max from
`61.7 MeV` to `43.0 MeV`. Heavy sectors remain at the few-to-tens of MeV level,
with some expected tradeoff in charm-light singlet ground states. The remaining
edge now looks like the paper's known pseudoscalar/chiral-annihilation
sensitivity, not a broad contact blow-up.

## 2026-05-05 Light-Sector Contrast Checkpoint

The raw sector comparison now makes the debugging clue explicit. Heavy,
heavy-light, strange, and isovector sectors are already on the few-to-tens of
MeV scale in the refreshed `docs/residual_reports/scorecard.md`; raw isoscalar
is the outlier because the digitized Fig. 5 rows are mixed `n nbar / s sbar`
pairs, while the current plain residual report compares both partners to the
same unmixed `n nbar` prediction.

Added `docs/residual_reports/light_sector_audit.md` from
`scripts/audit_nonmixing_contact.jl`. Its two-branch sanity check compares the
lower isoscalar partner to an unmixed `n nbar` solve and the upper partner to an
unmixed `s sbar` solve. That drops the isoscalar diagnostic to `31.1 MeV` mean
absolute residual, with the large remaining failures concentrated in
`1^1S_0` and `2^1S_0` pseudoscalars. This supports the next physics target:
explicit isoscalar annihilation/flavor mixing, especially the GI P1/P2
pseudoscalar machinery, rather than a global light-sector retune.
