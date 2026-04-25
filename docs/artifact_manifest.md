# Godfrey--Isgur reproduction artifact manifest

Generated: 2026-04-25

## Files

- `data/seed/godfrey_isgur_seed_masses.csv`: seed mass table assembled from later papers quoting GI or GI-Original values. Contains 192 rows.
- `data/seed/godfrey_isgur_sources.csv`: source manifest with URLs and usage notes.
- `docs/orchestrator_task.md`: task brief for an orchestrator agent, including extraction, verification, implementation, and audit phases.
- `docs/paper_navigation.md`: fast map to high-value paper sections, figures, and extraction caveats.
- `scripts/validate_seed.py`: local schema check for the seed mass table.
- `paper/Godfrey-Isgur-1985.pdf`: local copy of the primary paper.
- `paper/text/godfrey_isgur_1985_prose.md`: preferred prose reference generated from raw extraction.
- `paper/text/pdftotext_layout.txt`: layout-preserving raw extraction for tables.
- `paper/text/pdftotext_bbox.html`: positional raw extraction for audit.
- `paper/text/pdftotext_words.tsv`: word-position raw extraction for table reconstruction experiments.
- `scripts/build_paper_prose.py`: regenerates the prose-only Markdown reference.

Intermediate generated references under `paper/text/` are retained for
provenance but are not recommended starting points.

## Important caveats

- The seed mass table is not a substitute for the original 1985 paper.
- Rows marked `quoted_GI_later_table` should be verified against the original Godfrey--Isgur PDF.
- Rows marked `GI_original_updated_*` come from open-heavy review tables and may not be exactly the 1985 published table values.
- Light-sector original GI values are intentionally absent from this seed table until direct PDF extraction is done.

## Suggested immediate next action

Execute Phase 0 and Phase 1 from `docs/orchestrator_task.md`, using the original PDF as primary authority.
