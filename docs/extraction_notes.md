# Extraction Notes

Record PDF text extraction, page-image extraction, table identification, and
all row-level disagreements here.

Expected generated files:

- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`
- `paper/vision_ocr/pages/*.md`
- `paper/vision_ocr/page_images/*.png`
- `paper/vision_ocr/column_crops/*.png`
- `paper/screenshots/spectrum_pages/*.png`
- `data/raw/digitized_tables/*`
- `data/raw/tables_from_page_images.csv`
- `data/raw/extraction_audit.csv` (header-only template until text-vs-image rows are logged)

## 2026-04-25 Setup

- Downloaded the primary paper to `paper/Godfrey-Isgur-1985.pdf`.
- Verified PDF metadata with `pdfinfo`: 43 pages, unencrypted, PDF 1.4.
- Generated initial text-extraction artifacts and a placeholder page-image root.

These artifacts were useful for early navigation but damaged equations. They
were removed after the 2026-05-05 vision-OCR replacement; use
`paper/vision_ocr/` instead.

## 2026-04-25 Spectrum Screenshots And First Tables

- Rendered PDF pages 6-10 to `paper/screenshots/spectrum_pages/`.
- Created `data/raw/digitized_tables/` with one folder per table.
- Digitized Table I into
  `data/raw/digitized_tables/table_i_confinement/table_i_confinement.csv`.
- Digitized Table II into
  `data/raw/digitized_tables/table_ii_parameters/table_ii_parameters.csv` and
  a provisional TOML copy.
- Added a low-confidence visible-row transcription for Table III at
  `data/raw/digitized_tables/table_iii_isoscalar_mixings/table_iii_visible_rows.provisional.csv`.

Table II and especially Table III still need visual verification before clean
promotion.

## 2026-04-25 Figure 3 Label Digitization

- Added `data/raw/digitized_figures/fig03_isovector_mesons/`.
- Transcribed Fig. 3 model-state labels into
  `data/raw/digitized_figures/fig03_isovector_mesons/figure_03_labels.csv`.
- Scope is model mass-bar labels only; experimental shaded bands and axis ticks
  are excluded.
- A crop check corrected the `0-+` label `3^1S_0(1.88)` and confirmed
  `1--` label `3^3S_1(2.00)`.
- Added clean comparison plot outputs:
  `data/raw/digitized_figures/fig03_isovector_mesons/figure_03_replot.svg`
  and `.png`, generated from the CSV without experimental hatching or bands.

## 2026-05-05 Vision OCR Pass

- Added `scripts/vision_ocr_paper.py`, a reproducible OpenAI vision-OCR pipeline:
  render each PDF page at 300 dpi, create overlapping left/right column crops,
  send a low-detail full-page image plus high-detail column crops to `gpt-4.1`,
  and cache one Markdown file per page.
- Generated `paper/vision_ocr/pages/page-001.md` through `page-043.md` and the
  aggregate `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`.
- Saved provenance images under `paper/vision_ocr/page_images/` and
  `paper/vision_ocr/column_crops/`; API usage is logged in
  `paper/vision_ocr/usage.jsonl`.
- Actual billed pass, including one repeated Appendix-A pilot page: 143,537
  input tokens and 72,014 output tokens with `gpt-4.1`, about $0.86 at the
  current $2/M input and $8/M output rates. The latest one-output-per-page set is
  140,309 input tokens and 70,189 output tokens, about $0.84.
- The result is substantially better than `pdftotext`/Tesseract for equations,
  but dense formulas still require visual audit before being treated as
  authoritative. The known Appendix-A (A13) drift from the pilot was manually
  corrected against the saved crop in `page-037.md` and the aggregate.

## 2026-04-25 Figures 4-9 Label Digitization

- Added first-pass model-state label CSVs for Figures 4-9 under
  `data/raw/digitized_figures/`.
- Added `scripts/plot_spectrum_digitizations.py` to regenerate clean
  comparison replots for Figures 4-9 from those CSV files.
- Scope remains model mass-bar labels only; experimental shaded bands, quoted
  errors, and omitted near-degenerate states from captions are not drawn.
- Figure 5 is particularly crowded and has several `medium` or `low`
  confidence rows that should receive crop-by-crop audit before promotion to
  clean data.
