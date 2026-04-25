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
