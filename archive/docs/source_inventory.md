# Source Inventory

## Primary Source

- S. Godfrey and N. Isgur, *Mesons in a relativized quark model with
  chromodynamics*, Phys. Rev. D 32, 189-231 (1985).
- DOI: `10.1103/PhysRevD.32.189`
- Status: primary authority for final reference values.
- Local PDF: `paper/Godfrey-Isgur-1985.pdf`
- Source URL used for local copy:
  `https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevD.32.189/fulltext`
- PDF metadata: 43 pages, letter page size, PDF 1.4.

## Derived Paper References

- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`: preferred local
  Markdown reference for reading, search, and future agent interaction. It is
  currently the most complete OCR-derived text and includes manual improvements,
  especially around tables.
- `paper/vision_ocr/pages/`: archived historical one-file-per-page OCR output.
  These files remain in place for provenance links, but they are not an
  independent source of truth.
- `paper/vision_ocr/page_images/` and `paper/vision_ocr/column_crops/`: rendered
  provenance images for equation/table audit.
- `paper/vision_ocr/usage.jsonl`: API usage/provenance log for the vision OCR
  run.

The next preferred cleanup is chapter/block based rather than page based: keep
chapters as manageable Markdown includes and store long objects (tables,
full-width equations, captions, and figures) as separate named blocks that
Quarto can assemble into the full paper view.

The previous `pdftotext` references were removed because they damaged equation
typography. The PDF and rendered vision-OCR crops remain the authority for
equations, signs, table alignment, and state labels.

## Seed Sources

The seed files are stored in `data/seed/`.

- `data/seed/godfrey_isgur_seed_masses.csv`
- `data/seed/godfrey_isgur_sources.csv`

These values are bootstrap data assembled from later papers quoting GI or
GI-original values. They are not final reference values until checked against
the primary PDF.

## Extraction Policy

Each extracted numerical value must record:

- source paper,
- page number,
- table number or table label,
- extraction method,
- confidence,
- raw row provenance.
