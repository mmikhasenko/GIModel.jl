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
- **`table_v_strong_decays.csv`** (220 rows), **`table_vi_photon_decays.csv`**
  (79 rows), **`table_vii_annihilation_em.csv`** (61 rows): canonical,
  from-page-image transcriptions of Tables V/VI/VII. Schema and footnotes are
  documented in **[`README_canonical_tables.md`](README_canonical_tables.md)**.
- `table_v_strong_decays/`: **superseded** piecemeal `*.provisional.csv` files
  (folded into `table_v_strong_decays.csv`). Kept only because
  `audit_table_v_decays.jl` and `audit_table_v_2s_1d_decays.jl` still read them;
  they will be removed once those audits are rewired to the canonical file. Do
  not use for new work — see that folder's `DIGITIZATION_NOTES.md` header.

## Rules

- Keep raw page text in each table folder.
- Keep structured CSV/TOML separate from raw text.
- Mark uncertain symbols and OCR-normalized entries explicitly.
- Canonical table CSVs are transcribed row-by-row from `page_images/` and cite
  the source page in each row's `page` column; prefer them over any earlier
  provisional extraction.
- Promote only audited rows to `data/clean/`.

