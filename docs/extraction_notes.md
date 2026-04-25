# Extraction Notes

Record PDF text extraction, page-image extraction, table identification, and
all row-level disagreements here.

Expected generated files:

- `paper/text/pdftotext_layout.txt`
- `paper/text/pdftotext_bbox.html`
- `paper/text/godfrey_isgur_1985.md`
- `paper/text/godfrey_isgur_1985_polished.md`
- `paper/pages/*.png`
- `data/raw/tables_from_pdf_text.csv`
- `data/raw/tables_from_page_images.csv`
- `data/raw/extraction_audit.csv`

## 2026-04-25 Setup

- Downloaded the primary paper to `paper/Godfrey-Isgur-1985.pdf`.
- Verified PDF metadata with `pdfinfo`: 43 pages, unencrypted, PDF 1.4.
- Generated `paper/text/pdftotext_layout.txt` with `pdftotext -layout`.
- Generated `paper/text/pdftotext_bbox.html` with `pdftotext -bbox-layout`.
- Generated `paper/text/godfrey_isgur_1985.md` with
  `scripts/build_paper_markdown.py`.
- Generated `paper/text/godfrey_isgur_1985_polished.md` with
  `scripts/polish_paper_markdown.py`.

The Markdown files are search references only. Table extraction must still be
checked against the PDF and rendered page images.
