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
- `data/raw/digitized_tables/`: table-specific raw snippets and provisional
  structured transcriptions.
- `data/raw/digitized_figures/`: figure-specific label transcriptions.
- `scripts/validate_seed.py`: seed schema validation.
- `scripts/build_paper_prose.py`: rebuilds the prose-only paper reference.
- `scripts/plot_spectrum_digitizations.py`: regenerates clean Fig. 4-9
  comparison replots from digitized figure CSVs.

Intermediate generated paper references exist under `paper/text/` for
provenance and experiments, but they are not recommended handoff entry points.

## Central Data Targets

The spectrum figures are the central numerical result to reproduce. Treat the
digitized Figure 3-9 model labels in `data/raw/digitized_figures/` as the
first-pass cross-check targets for any solver output, after auditing each row
against the original PDF image. The clean SVG/PNG replots are comparison aids;
the CSV files are the data source.

Table II is the central setup input. Its digitization at
`data/raw/digitized_tables/table_ii_parameters/table_ii_parameters.csv` records
the fitted model parameters and should be audited before any production solver
configuration is promoted to `data/clean/`.

## Immediate Workflow

1. Run seed validation:

   ```bash
   python3 scripts/validate_seed.py
   ```

2. Rebuild paper text references if the PDF changes:

   ```bash
   pdftotext -layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_layout.txt
   pdftotext -bbox-layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_bbox.html
   pdftotext -raw paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_raw.txt
   pdftotext -tsv paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_words.tsv
   python3 scripts/build_paper_prose.py
   ```

3. Render page images for table verification when needed:

   ```bash
   pdftoppm -r 300 -png paper/Godfrey-Isgur-1985.pdf paper/pages/gi
   ```

4. Create raw extraction CSVs in `data/raw/`.
5. Promote verified rows into `data/clean/` with provenance preserved.

Current first-pass digitizations live under `data/raw/digitized_tables/` and
`data/raw/digitized_figures/`. Table I and Table II have structured provisional
CSVs; Table III is a low-confidence visible-row transcription that needs
image-audited cleanup. Figures 3-9 have first-pass model-label CSVs and clean
comparison replots.

## Authority Rules

- The original 1985 paper is the primary authority.
- The Markdown paper files are navigation aids, not authorities.
- Later quoted tables are seed data only.
- Every numerical value must carry source, table/page when available,
  extraction method, and confidence.
- Extraction data and cleaned physics data stay separate.
- Discrepancies are classified before code or data changes are made.
