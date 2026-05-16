# Paper Navigation

Use this as a fast map into the Godfrey-Isgur paper while coding or extracting
tables. For reading and search, start with
`paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`; it is the most complete
current OCR-derived Markdown source and has better table cleanup than the
initial page files. Always verify equations and table values against
`paper/Godfrey-Isgur-1985.pdf` or the saved `paper/vision_ocr/page_images/` and
`paper/vision_ocr/column_crops/`.

## Core Model

- PDF pages 1-2: introduction and motivation.
- PDF pages 2-4: soft-QCD Hamiltonian, relativized kinetic term, potential
  structure, smearing, running coupling, and numerical procedure.
- PDF page 5: Table II, published spectroscopy parameters.
- PDF pages 35-36: Appendix A, effective relativistic potentials and smearing
  prescriptions.
- PDF pages 37-38: Appendix B, wave-function and flavor conventions.

## Spectrum Figures

- PDF pages 5-6: light nonstrange and strange mesons.
- PDF page 7: isoscalar mesons and Table III mixing information.
- PDF page 8: charmonia and charmed mesons, including `c qbar` and `c sbar`
  mixing-angle notes.
- PDF page 9: bottomonia and b-flavored mesons, including `b qbar`, `b sbar`,
  and `b cbar` mixing-angle notes.

## Mixing And Unequal Masses

- PDF page 4: tensor and antisymmetric spin-orbit off-diagonal effects are
  treated by diagonalizing a mass matrix after the main Hamiltonian solve.
- PDF page 4: antisymmetric spin-orbit mixing appears only when quark masses
  are unequal.
- PDF page 8: charmed-sector mixing angles for `c dbar`, `c ubar`, and
  `c sbar`.
- PDF page 9: b-flavored mixing angles for `b qbar`, `b sbar`, and `b cbar`.
- PDF pages 37-38: flavor basis and state-vector conventions.

## Tables To Extract First

- Table II, PDF page 5: model parameters.
- Table III, PDF page 11 in the PDF text flow: mixed-isoscalar composition.
- Figures 6-9, PDF pages 8-9: first target spectra for `ccbar`, `cqbar`,
  `csbar`, `bbbar`, `bqbar`, `bsbar`, and `bcbar`.

Rendered spectrum pages are available under `paper/screenshots/spectrum_pages/`.
Digitized spectrum labels begin under `data/raw/digitized_figures/`.

## Caution Zones

- Markdown references are convenient for search but not reliable enough for
  final numbers.
- `paper/vision_ocr/pages/` is archived historical page OCR from a layout pass
  that was not fully successful. Prefer the full-paper Markdown unless auditing
  OCR history.
- The vision-OCR pass is much better than the removed `pdftotext` artifacts, but
  dense equations and crowded tables still need crop-level audit.
- The paper uses isospin symmetry for heavy-light `Q qbar` doublets; precision
  `u` vs `d` splittings are outside the baseline 1985 reproduction.
