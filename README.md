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
- `docs/source_inventory.md`: bibliographic metadata and source authority
  policy.
- `docs/extraction_notes.md`: running notes for PDF/text/image extraction.
- `docs/formula_map.md`: required map from implementation terms to the paper.
- `docs/conventions.md`: spectroscopic, spin, sector, and basis conventions.
- `data/seed/godfrey_isgur_seed_masses.csv`: bootstrap mass table, not final
  authority.
- `data/seed/godfrey_isgur_sources.csv`: manifest for the seed sources.
- `paper/Godfrey-Isgur-1985.pdf`: primary authority.
- `paper/text/pdftotext_layout.txt`: raw layout text extraction.
- `paper/text/pdftotext_bbox.html`: positional text extraction.
- `paper/text/godfrey_isgur_1985.md`: searchable Markdown reference derived
  from layout text.
- `paper/text/godfrey_isgur_1985_polished.md`: polished Markdown reference for
  smoother reading and future agent interaction.
- `scripts/validate_seed.py`: seed schema validation.
- `scripts/build_paper_markdown.py`: rebuilds the Markdown paper reference.
- `scripts/polish_paper_markdown.py`: rebuilds the polished paper reference.

## Immediate Workflow

1. Run seed validation:

   ```bash
   python3 scripts/validate_seed.py
   ```

2. Rebuild paper text references if the PDF changes:

   ```bash
   pdftotext -layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_layout.txt
   pdftotext -bbox-layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_bbox.html
   python3 scripts/build_paper_markdown.py
   python3 scripts/polish_paper_markdown.py
   ```

3. Render page images for table verification when needed:

   ```bash
   pdftoppm -r 300 -png paper/Godfrey-Isgur-1985.pdf paper/pages/gi
   ```

4. Create raw extraction CSVs in `data/raw/`.
5. Promote verified rows into `data/clean/` with provenance preserved.

## Authority Rules

- The original 1985 paper is the primary authority.
- The Markdown paper files are navigation aids, not authorities.
- Later quoted tables are seed data only.
- Every numerical value must carry source, table/page when available,
  extraction method, and confidence.
- Extraction data and cleaned physics data stay separate.
- Discrepancies are classified before code or data changes are made.
