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

- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`: preferred reference for
  reading, search, and future agent interaction.
- `paper/vision_ocr/pages/`: one Markdown transcription per PDF page.
- `paper/vision_ocr/page_images/` and `paper/vision_ocr/column_crops/`: rendered
  provenance images for equation/table audit.
- `paper/vision_ocr/usage.jsonl`: API usage/provenance log for the vision OCR
  run.

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
