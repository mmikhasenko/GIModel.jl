# Digitized Tables

One folder per table, keeping raw extraction snippets separate from structured
transcriptions.

## Current Status

- `table_i_confinement/`: digitized from PDF page 1.
- `table_ii_parameters/`: digitized from PDF page 5; parameter symbols still
  need visual verification against the PDF/page image.
- `table_iii_isoscalar_mixings/`: raw extraction and provisional visible-row
  transcription from PDF pages 11-12; this table needs image-audited cleanup
  before promotion to `data/clean/`.

## Rules

- Keep raw page text in each table folder.
- Keep structured CSV/TOML separate from raw text.
- Mark uncertain symbols and OCR-normalized entries explicitly.
- Promote only audited rows to `data/clean/`.

