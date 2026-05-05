# Godfrey--Isgur reproduction artifact manifest

Generated: 2026-04-25

## Files

- `data/seed/godfrey_isgur_seed_masses.csv`: seed mass table assembled from later papers quoting GI or GI-Original values. Contains 192 rows.
- `data/seed/godfrey_isgur_sources.csv`: source manifest with URLs and usage notes.
- `docs/orchestrator_task.md`: task brief for an orchestrator agent, including extraction, verification, implementation, and audit phases.
- `docs/agent_handoff.md`: current project status, subagent split, and immediate next milestone.
- `docs/paper_navigation.md`: fast map to high-value paper sections, figures, and extraction caveats.
- `scripts/validate_seed.py`: local schema check for the seed mass table.
- `paper/Godfrey-Isgur-1985.pdf`: local copy of the primary paper.
- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`: preferred
  reading/search reference generated from rendered page and column images.
- `paper/vision_ocr/pages/`: one Markdown transcription per PDF page.
- `paper/vision_ocr/page_images/` and `paper/vision_ocr/column_crops/`:
  rendered provenance images for audit.
- `paper/vision_ocr/usage.jsonl`: API usage/provenance log for the vision OCR
  run.
- `paper/screenshots/spectrum_pages/`: rendered PNG pages containing spectrum figures.
- `data/raw/digitized_tables/`: one-folder-per-table raw snippets and provisional structured transcriptions.
- `data/raw/digitized_figures/`: one-folder-per-figure model-label transcriptions.
- `scripts/vision_ocr_paper.py`: regenerates the vision-OCR Markdown reference.

The older `pdftotext` artifacts were removed after the vision-OCR replacement
because they damaged equation typography.

## Important caveats

- The seed mass table is not a substitute for the original 1985 paper.
- Rows marked `quoted_GI_later_table` should be verified against the original Godfrey--Isgur PDF.
- Rows marked `GI_original_updated_*` come from open-heavy review tables and may not be exactly the 1985 published table values.
- Light-sector original GI values are intentionally absent from this seed table until direct PDF extraction is done.

## Suggested immediate next action

Execute Phase 0 and Phase 1 from `docs/orchestrator_task.md`, using the original PDF as primary authority.
