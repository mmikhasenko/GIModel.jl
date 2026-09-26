# GIModel.jl

Local reproduction of the Godfrey-Isgur relativized quark model for meson
masses. The original paper at `paper/Godfrey-Isgur-1985.pdf` is the authority;
OCR Markdown and extracted CSVs are navigation/provenance aids.

GIModel computes meson spectra and wavefunctions with independent harmonic-oscillator
and finite-difference solvers. It has three external runtime dependencies
(KrylovKit, QuadGK and SpecialFunctions), plus Julia standard libraries.
Julia 1.11 or later is required; release checks use Julia 1.11.

From a checkout, install into your own Julia environment:

```julia
using Pkg
Pkg.develop(path="/absolute/path/to/GIModel.jl")
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
spec = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(1))
```

The preset path is independent of your working directory. The ordinary spectrum
uses [18 numerical physics inputs](docs/model_inputs.md); optional isoscalar mixing
and transition calculations have additional inputs. For wave operations and
worked examples, see [examples](examples/README.md).

The repository has five first-class deliverables:

- **GIModel** (repository root) — pure computation. Mesons are specified by
  quark flavors (`Meson(mq, :c, :b)`), levels to compute by explicit
  `n^{2S+1}L_J` multiplets (`spectrum_levels`), and `compute_spectrum` returns
  an organized `Spectrum` with contribution breakdowns and intra-meson mixing.
  It owns the Schrödinger solver, physical-state components, wave
  representations, and generic radial/momentum overlap operations. It knows
  nothing about decay operators or the paper comparison. Isoscalar annihilation
  mixing remains here because it changes the mass spectrum.
- **QuarkModelTransitions** (`QuarkModelTransitions/`) — bare quark-model
  transition operators and observables consuming GIModel states and waves. It
  owns channels, helicities, partial waves, coherent state composition,
  matrix elements, radiative and annihilation observables, charge radii, widths,
  and the frozen Table IV/V reference backend.
- **GIPaper** (`GIPaper/`) — the comparison layer. It owns the digitized paper
  data, maps reference CSV rows to mesons (`reference_meson`, no fallback
  masses), runs `compare_reference`, applies the Table III annihilation
  prescriptions, and writes residual reports.
- **Report** (`report/`, separate Overleaf repository) — the code-free, high-level account,
  *Recomputing the Godfrey--Isgur Relativized Quark Model*. It is a RevTeX
  (PRD) paper whose local sources live in `report/`:
  motivation, computation and architecture, implementation (the FD-then-HO
  route and the corrections that taught the physics), usage, four decades of
  post-1985 developments, and outlook. Clone it into `report/` with
  `git clone https://git@git.overleaf.com/6aa00b740749f422ca1291fe report`
  (git-ignored here). Its `archive/` holds the previous Quarto edition, the one
  with the sector-spectrum and wavefunction plots.
- **Learning track** (`LearningTrack/`) — **AGI: Agentic Godfrey--Isgur**, a
  nine-sheet pen-and-paper course from elementary radial quantum mechanics to
  the complete paper-order spectrum algorithm, plus a question-led bridge to a
  future hands-on numerical course. Each theory problem includes a worked
  solution and a post-solution concept check; the LaTeX sources build to ten
  standalone PDFs with `make`.

## Current Map

- `src/GIModel.jl` is the computation package entry point.
- `QuarkModelTransitions/src/QuarkModelTransitions.jl` is the transition
  package entry point and depends one-way on GIModel.
- `GIPaper/src/GIPaper.jl` is the comparison package entry point; its
  [README](GIPaper/README.md) defines the package boundary and separate script
  environment.
- `test/runtests.jl` gates the pure numerics;
  `QuarkModelTransitions/test/runtests.jl` gates transition physics; and
  `GIPaper/test/runtests.jl` gates the reference comparison.
- `scripts/verify_project.sh` runs the current full gate (all three packages).
- `GIPaper/test/data_validation.jl` owns package-level CSV, TOML, and
  provenance invariants and runs as part of `Pkg.test()`.
- `GIPaper/extraction/data_checks.py promote-clean` remains the write-oriented
  data-promotion utility; scorecard generation is Julia-native via
  `GIPaper/scripts/score_annihilation.jl`.
- `GIPaper/scripts/run_all_spectrum_checks.jl` regenerates sector residual
  reports and the compact scorecard.
- `GIPaper/scripts/analyze_heavy_quarkonium.jl` regenerates heavy-quarkonium
  diagnostics.
- `GIPaper/scripts/audit_nonmixing_contact.jl` regenerates the
  non-mixing/contact scorecards.
- `examples/` collects worked examples of GIModel driven as a physics tool —
  public API only, no paper data — in their own environment
  (`examples/Project.toml`: GIModel + QuarkModelTransitions + CairoMakie + PlutoUI). See
  [examples/README.md](examples/README.md), which doubles as the "how is this
  package used" walkthrough. Three curated examples, plus
  `examples/played_with_model.jl`, the uncurated scratch notebook (the one file
  there that uses GIPaper, and runs on its environment):
  - `examples/adaptive_ho_refinement.jl`, a Pluto notebook that shows how the
    mesh-free HO solver optimizes `β`, enlarges the basis until two successive
    energy checks pass, records its convergence certificate, and fails loudly
    at unresolved limits; it ends with an interactive comprehension quiz.
  - `examples/heavy_quark_transition.jl`, a Pluto notebook that dials `m_Q` from
    charm to bottom — `Q Q̄` / `Q q̄` / `Q s̄` switch, the level scheme in
    `n^{2S+1}L_J` and in `J^P` (axis following the spectrum, so only the shape
    moves), and an optional recorded GIF of either sweep.
  - `examples/chi_c_annihilation_widths.jl`, a script computing `χ_c0`/`χ_c2`
    two-gluon and two-photon widths and splitting their ratio into the `15/4`
    spin algebra times the J-dependent distortion of the wave at the origin.
- [Rate input ledger](GIPaper/docs/rate_input_ledger.md) traces the numerical
  inputs to Tables V, VI and VII, with generated per-row records.
- [Discoverability](docs/discoverability.md) describes the help conventions and links
  to the [public API documentation graph](docs/discoverability_graph.md).
- `docs/formula_map.md` maps active code paths to paper equations.
- `docs/original_1985_algorithm_audit.md` compares the paper's three-stage,
  mesh-free HO spectrum algorithm with the completed implementation.
- `docs/incident_followup_audit.md` records the forensic review of the abandoned
  migration, the additional stale infrastructure found, and what was repaired or
  deliberately left as history.
- `docs/work_plan.md` records closure of the Table VI photon-decay work: all
  79 canonical rows are computed, with quantitative residuals retained.
- `diff_support/` contains the pre-implementation audit, equations, numerical
  probes, backend research, and risk register for possible parameter
  differentiation.
- `docs/paper_gap_ledger.md` is the current “what remains vs the paper” list.
- `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` is the preferred
  searchable Markdown rendering of the paper.
- `docs/code_architecture.md`, `docs/conventions.md`, and
  `docs/appendix_a_from_paper.md` describe the implementation conventions.

Extraction utilities are intentionally separate from the core gate:
`GIPaper/extraction/vision_ocr_paper.py`,
`GIPaper/extraction/plot_figure3_digitization.py`, and
`GIPaper/extraction/plot_spectrum_digitizations.py`.

Concluded material from the reproduction phase (poster, Table III forensics,
early research notes) lives under `archive/` — see
[archive/README.md](archive/README.md).

## Possible Improvements Beyond the 1985 Paper

The strong-decay tables in the paper use an analytic SU(6), single-oscillator-
scale approximation with `beta = 0.40 GeV`. `QuarkModelTransitions` now also
applies the underlying Eq. (19) operator `g σ·q + h σ·p'` directly between
resolved physical meson wavefunctions. The native path retains radial nodes,
state-dependent length scales, coherent spectroscopic/flavor mixing, and the
solver-native HO or FD representation; it returns helicities, all allowed
partial waves, and an Eq. (C2) width. This is deliberately separate from the
frozen Table IV/V reproduction because it goes beyond the approximation used
for the paper's numerical strong-decay results.

Another possible extension is differentiable mass prediction with respect to
continuous model parameters. The preliminary route-specific numerical evidence
and the risks around eigenvalue crossings, adaptive beta selection, and state
identity are preserved under `diff_support/`; they are research notes, not a
scheduled implementation.

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
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=QuarkModelTransitions -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
julia GIPaper/scripts/run_all_spectrum_checks.jl
julia GIPaper/scripts/analyze_heavy_quarkonium.jl
```

Pure-model usage without any reference data:

```julia
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
spec = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(2))
```

GIModel/GIPaper package checks: `bash scripts/verify_packages.sh`.

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
