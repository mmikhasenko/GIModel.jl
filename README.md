# Godfrey-Isgur Reproduction

Local, reproducible reproduction of the Godfrey-Isgur relativized quark model
for meson masses.

The working plan lives in `docs/orchestrator_task.md`. The primary paper is
stored at `paper/Godfrey-Isgur-1985.pdf`. The seed database in `data/seed/` is
only a bootstrap from later sources and must be verified against the original
1985 paper before being treated as reference data.

## Handoff Map

Use these files to resume work quickly:

- `docs/orchestrator_task.md`: phase plan, validation gates, and acceptance
  criteria.
- `docs/agent_handoff.md`: current status, recommended subagent split, and
  immediate next milestone.
- `docs/source_inventory.md`: bibliographic metadata and source authority
  policy.
- `docs/extraction_notes.md`: running notes for PDF/text/image extraction.
- `docs/paper_navigation.md`: fast map to paper sections, tables, spectrum
  figures, and caution zones.
- `docs/formula_map.md`: required map from implementation terms to the paper.
- `docs/conventions.md`: spectroscopic, spin, sector, and basis conventions.
- `data/table_ii_parameters.csv`: central Table II parameter digitization.
- `data/parameters.provisional.toml`: provisional solver-facing parameter file
  copied from the Table II digitization.
- `data/reference_spectrum_charmonium.csv`: central Fig. 6 `ccbar` reference
  spectrum.
- `data/reference_spectrum_bottomonium.csv`: central Fig. 8 `bbbar` reference
  spectrum.
- `data/reference_spectrum_*.csv`: top-level copies of Fig. 3-9 model-label
  spectra for convenient use by scripts and tests.
- `data/seed/godfrey_isgur_seed_masses.csv`: bootstrap mass table, not final
  authority.
- `data/seed/godfrey_isgur_sources.csv`: manifest for the seed sources.
- `paper/Godfrey-Isgur-1985.pdf`: primary authority.
- `paper/text/godfrey_isgur_1985_prose.md`: preferred reading/search reference.
- `paper/text/pdftotext_layout.txt`: raw layout extraction for table work.
- `paper/text/pdftotext_bbox.html`: positional extraction for audit work.
- `paper/text/pdftotext_words.tsv`: word-position extraction for table
  reconstruction experiments.
- `paper/screenshots/spectrum_pages/`: rendered spectrum pages, PDF pages 6-10.
- `data/raw/digitized_tables/`: reference/provenance copies for table-specific
  raw snippets and structured transcriptions.
- `data/raw/digitized_figures/`: reference/provenance copies for figure-specific
  label transcriptions and replots.
- `scripts/validate_seed.py`: seed schema validation.
- `scripts/build_paper_prose.py`: rebuilds the prose-only paper reference.
- `scripts/plot_spectrum_digitizations.py`: regenerates clean Fig. 4-9
  comparison replots from digitized figure CSVs.

Intermediate generated paper references exist under `paper/text/` for
provenance and experiments, but they are not recommended handoff entry points.

## Central Data Targets

Use the top-level files in `data/` as the working inputs:

- `data/table_ii_parameters.csv`: Table II parameter digitization.
- `data/parameters.provisional.toml`: provisional parameter TOML derived from
  Table II.
- `data/reference_spectrum_charmonium.csv`: Fig. 6 `ccbar` target spectrum.
- `data/reference_spectrum_bottomonium.csv`: Fig. 8 `bbbar` target spectrum.
- `data/reference_spectrum_isovector.csv`: Fig. 3 isovector spectrum.
- `data/reference_spectrum_strange.csv`: Fig. 4 strange spectrum.
- `data/reference_spectrum_isoscalar.csv`: Fig. 5 isoscalar spectrum.
- `data/reference_spectrum_charmed.csv`: Fig. 7 charmed/charmed-strange
  spectra.
- `data/reference_spectrum_b_flavored.csv`: Fig. 9 bottom-light,
  bottom-strange, and bottom-charm spectra.

The deeper `data/raw/digitized_tables/` and `data/raw/digitized_figures/`
folders are kept as provenance/reference material. They preserve the original
per-table and per-figure working context, raw snippets, and replots, but routine
solver and validation scripts should start from the top-level `data/` files.

## Immediate Workflow

1. Run seed validation:

   ```bash
   python3 scripts/validate_seed.py
   ```

2. Run the current local spectrum baseline:

   ```bash
   julia --project=. scripts/run_baseline_solver.jl
   ```

   This writes:

   - `docs/residual_reports/ccbar_baseline.md`
   - `docs/residual_reports/bbbar_baseline.md`

   The current baseline is intentionally diagnostic: it uses the Table II
   quark masses, `b`, `c`, the Fig. 2 running Coulomb ansatz, semirelativistic
   kinetic energy, and the smeared S-wave contact hyperfine term. It does not
   yet include the full GI smearing/nonlocal potential, tensor interaction,
   spin-orbit terms, or mixing.

3. Run all top-level reference spectrum checks:

   ```bash
   julia --project=. scripts/run_all_spectrum_checks.jl
   ```

   This regenerates sector reports under `docs/residual_reports/` for every
   `data/reference_spectrum_*.csv` file. Heavy-light sectors use the
   `quark_content` column to choose unequal constituent masses.

4. Run Julia tests:

   ```bash
   julia --project=. test/runtests.jl
   ```

5. Rebuild paper text references if the PDF changes:

   ```bash
   pdftotext -layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_layout.txt
   pdftotext -bbox-layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_bbox.html
   pdftotext -raw paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_raw.txt
   pdftotext -tsv paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_words.tsv
   python3 scripts/build_paper_prose.py
   ```

6. Render page images for table verification when needed:

   ```bash
   pdftoppm -r 300 -png paper/Godfrey-Isgur-1985.pdf paper/pages/gi
   ```

7. Create raw extraction CSVs in `data/raw/`.
8. Promote verified rows into `data/clean/` with provenance preserved.

Current first-pass digitizations are copied to the top level of `data/`.
Reference/provenance copies remain under `data/raw/digitized_tables/` and
`data/raw/digitized_figures/`. Table I and Table III are retained only in the
raw reference folders for now; the central workflow needs Table II and the
Figure 3-9 model spectra.

## Authority Rules

- The original 1985 paper is the primary authority.
- The Markdown paper files are navigation aids, not authorities.
- Later quoted tables are seed data only.
- Every numerical value must carry source, table/page when available,
  extraction method, and confidence.
- Extraction data and cleaned physics data stay separate.
- Discrepancies are classified before code or data changes are made.
