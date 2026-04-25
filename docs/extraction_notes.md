# Extraction Notes

Record PDF text extraction, page-image extraction, table identification, and
all row-level disagreements here.

Expected generated files:

- `paper/text/pdftotext_layout.txt`
- `paper/text/pdftotext_raw.txt`
- `paper/text/pdftotext_bbox.html`
- `paper/text/pdftotext_words.tsv`
- `paper/text/godfrey_isgur_1985_prose.md`
- `paper/pages/*.png`
- `paper/screenshots/spectrum_pages/*.png`
- `data/raw/digitized_tables/*`
- `data/raw/tables_from_pdf_text.csv`
- `data/raw/tables_from_page_images.csv`
- `data/raw/extraction_audit.csv`

## 2026-04-25 Setup

- Downloaded the primary paper to `paper/Godfrey-Isgur-1985.pdf`.
- Verified PDF metadata with `pdfinfo`: 43 pages, unencrypted, PDF 1.4.
- Generated `paper/text/pdftotext_layout.txt` with `pdftotext -layout`.
- Generated `paper/text/pdftotext_bbox.html` with `pdftotext -bbox-layout`.
- Generated `paper/text/pdftotext_raw.txt` with `pdftotext -raw`.
- Generated `paper/text/pdftotext_words.tsv` with `pdftotext -tsv`.
- Generated `paper/text/godfrey_isgur_1985_prose.md` with
  `scripts/build_paper_prose.py`.

The raw extraction has much better reading order than the layout extraction.
The prose Markdown is currently the preferred reading/navigation file. Table
extraction must still be checked against the PDF and rendered page images.
Intermediate generated Markdown experiments are retained in `paper/text/` for
provenance, but should not be used as default references.

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
