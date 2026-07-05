# Godfrey-Isgur Reproduction

Local reproduction of the Godfrey-Isgur relativized quark model for meson
masses. The original paper at `paper/Godfrey-Isgur-1985.pdf` is the authority;
OCR Markdown and extracted CSVs are navigation/provenance aids.

The repository holds two Julia packages:

- **GIModel** (repository root) — pure computation. Mesons are specified by
  quark flavors (`Meson(mq, :c, :b)`), levels to compute by explicit
  `n^{2S+1}L_J` multiplets (`spectrum_levels`), and `compute_spectrum` returns
  an organized `Spectrum` with contribution breakdowns and intra-meson mixing.
  It knows nothing about the paper comparison.
- **GIPaper** (`GIPaper/`) — the comparison layer. It owns the digitized paper
  data, maps reference CSV rows to mesons (`reference_meson`, no fallback
  masses), runs `compare_reference`, applies the Table III annihilation
  prescriptions, and writes residual reports.

## Current Map

- `src/GIModel.jl` is the computation package entry point.
- `GIPaper/src/GIPaper.jl` is the comparison package entry point.
- `test/runtests.jl` gates the pure numerics; `GIPaper/test/runtests.jl` gates
  the reference comparison.
- `scripts/verify_project.sh` runs the current full gate (both packages).
- `GIPaper/scripts/data_checks.py` is the only Python data/check entry point:
  - `python3 GIPaper/scripts/data_checks.py promote-clean`
  - `python3 GIPaper/scripts/data_checks.py validate`
  - `python3 GIPaper/scripts/data_checks.py score-annihilation`
- `GIPaper/scripts/run_all_spectrum_checks.jl` regenerates sector residual
  reports and the compact scorecard.
- `GIPaper/scripts/analyze_heavy_quarkonium.jl` regenerates heavy-quarkonium
  diagnostics.
- `GIPaper/scripts/audit_nonmixing_contact.jl` regenerates the
  non-mixing/contact scorecards.
- `docs/formula_map.md` maps active code paths to paper equations.
- `docs/paper_gap_ledger.md` is the current “what remains vs the paper” list.
- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` is the preferred
  searchable Markdown rendering of the paper.
- `docs/code_architecture.md`, `docs/conventions.md`, and
  `docs/appendix_a_from_paper.md` describe the implementation conventions.

Extraction utilities are intentionally separate from the core gate:
`GIPaper/scripts/vision_ocr_paper.py`, `GIPaper/scripts/plot_figure3_digitization.py`,
and `GIPaper/scripts/plot_spectrum_digitizations.py`.

Concluded material from the reproduction phase (poster, Quarto report, Table III
forensics, early research notes) lives under `archive/` — see
[archive/README.md](archive/README.md).

## Data

Model configuration (GIModel):

- `data/parameters.provisional.toml`

Paper reference data (GIPaper):

- `GIPaper/data/reference_spectrum_*.csv`

Promoted audited data:

- `GIPaper/data/clean/masses.csv`
- `GIPaper/data/clean/mixings.csv`
- `GIPaper/data/clean/parameters.toml`

Raw provenance stays under `GIPaper/data/raw/` and `paper/vision_ocr/`.

## Common Commands

```bash
python3 GIPaper/scripts/data_checks.py validate
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
julia GIPaper/scripts/run_all_spectrum_checks.jl
julia GIPaper/scripts/analyze_heavy_quarkonium.jl
```

Pure-model usage without any reference data:

```julia
using GIModel
params, mq = load_parameters_and_quark_masses("data/parameters.provisional.toml")
spec = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(2))
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
