# Godfrey-Isgur Reproduction

Local reproduction of the Godfrey-Isgur relativized quark model for meson
masses. The original paper at `paper/Godfrey-Isgur-1985.pdf` is the authority;
OCR Markdown and extracted CSVs are navigation/provenance aids.

## Current Map

- `src/GIModel.jl` is the Julia package entry point.
- `test/runtests.jl` is the main correctness gate.
- `scripts/verify_project.sh` runs the current full gate.
- `scripts/data_checks.py` is the only Python data/check entry point:
  - `python3 scripts/data_checks.py promote-clean`
  - `python3 scripts/data_checks.py validate`
  - `python3 scripts/data_checks.py score-annihilation`
- `scripts/run_all_spectrum_checks.jl` regenerates sector residual reports and
  the compact scorecard.
- `scripts/analyze_heavy_quarkonium.jl` regenerates heavy-quarkonium diagnostics.
- `scripts/audit_nonmixing_contact.jl` regenerates the non-mixing/contact
  scorecards.
- `docs/formula_map.md` maps active code paths to paper equations.
- `docs/paper_gap_ledger.md` is the current “what remains vs the paper” list.
- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` is the preferred
  searchable Markdown rendering of the paper.
- `docs/code_architecture.md`, `docs/conventions.md`, and
  `docs/appendix_a_from_paper.md` describe the implementation conventions.

Extraction utilities are intentionally separate from the core gate:
`scripts/vision_ocr_paper.py`, `scripts/plot_figure3_digitization.py`, and
`scripts/plot_spectrum_digitizations.py`.

## Data

Working solver inputs:

- `data/parameters.provisional.toml`
- `data/reference_spectrum_*.csv`

Promoted audited data:

- `data/clean/masses.csv`
- `data/clean/mixings.csv`
- `data/clean/parameters.toml`

Raw provenance stays under `data/raw/` and `paper/vision_ocr/`.

## Common Commands

```bash
python3 scripts/data_checks.py validate
julia --project=. test/runtests.jl
julia scripts/run_all_spectrum_checks.jl
julia scripts/analyze_heavy_quarkonium.jl
```

Full local gate:

```bash
bash scripts/verify_project.sh
```

## Authority Rules

- The 1985 paper is authoritative.
- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` is the preferred OCR
  Markdown context; archived page OCR under `paper/vision_ocr/pages/` is
  historical provenance only.
- Raw extraction and clean physics data stay separate.
- Every promoted numerical value should retain provenance and confidence.
