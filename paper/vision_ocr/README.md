# Vision OCR Source Policy

This folder contains OCR and image-provenance artifacts for the Godfrey-Isgur
1985 paper. The current OCR pass was useful, but the page-level strategy was
not fully successful: the paper is mostly two-column text, while many equations,
figures, captions, and tables span the full page. Full-page OCR and fixed
left/right column crops therefore each miss important structure.

## Current Authority Order

1. `../Godfrey-Isgur-1985.pdf` is the primary source for final values.
2. `godfrey_isgur_1985_vision_ocr.md` is the best current Markdown reading and
   search source. It contains manual improvements, especially for tables, that
   are not present in the original per-page OCR files.
3. `page_images/` and `column_crops/` are visual provenance for audit.
4. `pages/` is archived historical output from the initial OCR project. These
   files are useful for provenance and comparison only; do not treat them as an
   independent source of truth.
5. `usage.jsonl` is OCR-run metadata, not source text.

## Future Shape

The preferred next structure is chapter/block based rather than page based:

- split the paper into chapter-level Markdown files,
- keep long objects such as tables, full-width equations, figures, and captions
  as separate named blocks,
- assemble the full paper with Quarto includes so the full Markdown/HTML/PDF
  view is generated from the same blocks,
- keep structured manual digitizations under `data/raw/digitized_tables/`, not
  inside OCR reading text.

This keeps language-model context manageable without creating several competing
versions of the paper.
