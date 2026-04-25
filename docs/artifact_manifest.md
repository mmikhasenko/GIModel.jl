# Godfrey--Isgur reproduction artifact manifest

Generated: 2026-04-25

## Files

- `data/seed/godfrey_isgur_seed_masses.csv`: seed mass table assembled from later papers quoting GI or GI-Original values. Contains 192 rows.
- `data/seed/godfrey_isgur_sources.csv`: source manifest with URLs and usage notes.
- `docs/orchestrator_task.md`: task brief for an orchestrator agent, including extraction, verification, implementation, and audit phases.
- `docs/paper_navigation.md`: fast map to high-value paper sections, figures, and extraction caveats.
- `scripts/validate_seed.py`: local schema check for the seed mass table.
- `paper/Godfrey-Isgur-1985.pdf`: local copy of the primary paper.
- `paper/text/pdftotext_layout.txt`: layout-preserving raw text extraction.
- `paper/text/pdftotext_bbox.html`: positional raw text extraction.
- `paper/text/godfrey_isgur_1985.md`: searchable Markdown reference derived from the layout extraction.
- `paper/text/godfrey_isgur_1985_polished.md`: conservative OCR-polished Markdown reference.
- `scripts/build_paper_markdown.py`: regenerates the Markdown reference.
- `scripts/polish_paper_markdown.py`: regenerates the polished Markdown reference.

## Important caveats

- The seed mass table is not a substitute for the original 1985 paper.
- Rows marked `quoted_GI_later_table` should be verified against the original Godfrey--Isgur PDF.
- Rows marked `GI_original_updated_*` come from open-heavy review tables and may not be exactly the 1985 published table values.
- Light-sector original GI values are intentionally absent from this seed table until direct PDF extraction is done.

## Suggested immediate next action

Execute Phase 0 and Phase 1 from `docs/orchestrator_task.md`, using the original PDF as primary authority.
