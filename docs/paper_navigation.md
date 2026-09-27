# Paper Navigation

> **Detail layer.** The manifest/dashboard (`docs/paper_manifest/*.toml`;
> `cd docs && make dashboard`) now anchors every equation, table, and figure to
> its PDF page of the article (PDF page = journal page − 188) — use it to jump
> from a unit to its page. This file stays as the human **reading map**.

Use this as a fast map into the Godfrey-Isgur paper while coding or extracting
tables: S. Godfrey and N. Isgur, *Phys. Rev. D* **32**, 189 (1985),
[doi:10.1103/PhysRevD.32.189](https://doi.org/10.1103/PhysRevD.32.189). The
article is not distributed with this repository; always verify equations and
table values against the published PDF.

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

Digitized spectrum labels begin under `GIPaper/data/raw/digitized_figures/`;
the recomputed tables and figures, in the paper's layout, are indexed in
[GIPaper/docs/paper_tables](../GIPaper/docs/paper_tables/README.md).

## Caution Zones

- OCR transcriptions are convenient for search but not reliable enough for
  final numbers; dense equations and crowded tables need an image-level audit.
- The paper uses isospin symmetry for heavy-light `Q qbar` doublets; precision
  `u` vs `d` splittings are outside the baseline 1985 reproduction.
