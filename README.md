# GIModel.jl

A Julia implementation of the Godfrey–Isgur relativized quark model for meson
spectra, wavefunctions and decays: S. Godfrey and N. Isgur, *Phys. Rev. D*
**32**, 189 (1985), [doi:10.1103/PhysRevD.32.189](https://doi.org/10.1103/PhysRevD.32.189).

GIModel computes meson spectra and wavefunctions with independent
harmonic-oscillator and finite-difference solvers. It has three external runtime
dependencies (KrylovKit, QuadGK and SpecialFunctions), plus Julia standard
libraries. Julia 1.11 or later is required; release checks use Julia 1.11.

## Quick start

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

## Packages

The repository holds three Julia packages with a one-way dependency chain,
GIPaper → QuarkModelTransitions → GIModel (GIPaper also uses GIModel directly).

- **GIModel** (repository root) — pure computation. Mesons are specified by
  quark flavors (`Meson(mq, :c, :b)`), levels by explicit `n^{2S+1}L_J`
  multiplets (`spectrum_levels`), and `compute_spectrum` returns an organized
  `Spectrum` with contribution breakdowns and intra-meson mixing. It owns the
  Schrödinger solvers, physical-state components, wave representations, and
  generic radial/momentum overlap operations. It knows nothing about decay
  operators or the paper comparison. Isoscalar annihilation mixing remains here
  because it changes the mass spectrum.
- **QuarkModelTransitions** ([`QuarkModelTransitions/`](QuarkModelTransitions/README.md))
  — quark-model transition operators and observables consuming GIModel states
  and waves: channels, helicities, partial waves, coherent state composition,
  matrix elements, radiative and annihilation observables, charge radii,
  widths, and the frozen Table IV/V reference backend.
- **GIPaper** ([`GIPaper/`](GIPaper/README.md)) — the comparison layer. It owns
  the digitized paper data, maps reference rows to mesons (`reference_meson`,
  no fallback masses), runs `compare_reference`, applies the Table III
  annihilation prescriptions, and writes residual reports. Its
  [paper tables](GIPaper/docs/paper_tables/README.md) recompute every numbered
  table and spectrum figure of the 1985 paper in the paper's own layout.

To install all three into one environment:

```julia
using Pkg
root = "/absolute/path/to/GIModel.jl"
Pkg.develop([PackageSpec(path=root),
    PackageSpec(path=joinpath(root, "QuarkModelTransitions")),
    PackageSpec(path=joinpath(root, "GIPaper"))])
```

Also in the repository:

- [`examples/`](examples/README.md) — worked examples of GIModel as a physics
  tool, in their own environment: adaptive HO refinement, a charm-to-bottom
  heavy-quark sweep, and χ_c two-gluon/two-photon widths.
- [`LearningTrack/`](LearningTrack/README.md) — a pen-and-paper course from
  radial quantum mechanics to the complete spectrum algorithm; the LaTeX
  sources build with `make`.

## Documentation

- [Model inputs](docs/model_inputs.md) — the 18 spectrum inputs and where each lives.
- [Code architecture](docs/code_architecture.md) and [conventions](docs/conventions.md).
- [Formula map](docs/formula_map.md) — active code paths against paper equations;
  [paper navigation](docs/paper_navigation.md) — a page map of the article;
  [Appendix A](docs/appendix_a_from_paper.md) and its
  [equation audit](docs/appendix_a_equation_audit.md).
- [Paper manifest](docs/paper_manifest/) — every equation, table and figure
  with its implementation status, checked by `GIPaper/scripts/check_manifest.jl`.
- [Original 1985 algorithm audit](docs/original_1985_algorithm_audit.md) — the
  paper's mesh-free HO spectrum algorithm against the implementation.
- [Paper gap ledger](docs/paper_gap_ledger.md) — what remains relative to the paper.
- [Discoverability](docs/discoverability.md) — help conventions and the
  [public API documentation graph](docs/discoverability_graph.md).
- [Rate input ledger](GIPaper/docs/rate_input_ledger.md) — numerical inputs to
  Tables V, VI and VII, with generated per-row records.
- [Strong decays tutorial](docs/strong_decays_tutorial.qmd) (Quarto).

## Beyond the 1985 paper

The strong-decay tables in the paper use an analytic SU(6), single-oscillator-
scale approximation with `beta = 0.40 GeV`. `QuarkModelTransitions` also
applies the underlying Eq. (19) operator `g σ·q + h σ·p'` directly between
resolved physical meson wavefunctions. The native path retains radial nodes,
state-dependent length scales, coherent spectroscopic/flavor mixing, and the
solver-native HO or FD representation; it returns helicities, all allowed
partial waves, and an Eq. (C2) width. It is deliberately separate from the
frozen Table IV/V reproduction because it goes beyond the approximation used
for the paper's numerical strong-decay results.

Transition calculations use `matrix_element(final, operator, initial)` and
`decay_width`, with masses and momenta derived from the GIModel states.
`mass_correction_factor` provides separate comparisons at another external mass
or momentum while holding the waves fixed. The
[transition guide](QuarkModelTransitions/README.md) covers photon and strong
emission, leptonic currents, annihilation, units, and correction factors.

## Data

- Model configuration (GIModel): `data/parameters.provisional.toml`.
- Paper reference spectra (GIPaper): `GIPaper/data/reference_spectrum_*.csv`.
- Promoted audited data: `GIPaper/data/clean/masses.csv`,
  `GIPaper/data/clean/mixings.csv`, `GIPaper/data/clean/parameters.toml`
  (a Table II transcription, not a solver configuration).
- Raw extraction evidence: `GIPaper/data/raw/`.

## Common commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=QuarkModelTransitions -e 'using Pkg; Pkg.test()'
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
julia GIPaper/scripts/run_all_spectrum_checks.jl
julia GIPaper/scripts/paper_tables/generate.jl
```

`Pkg.test()` runs each package's quick suite, including
[Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) quality checks; this is what
CI runs on every push and pull request. GIModel's solver-convergence sweeps
(`test/heavy/`, a few minutes) run only on request:

```bash
GI_HEAVY_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

Fast GIModel/GIPaper check: `bash scripts/verify_packages.sh`. Full gate, all
three packages with the heavy tests plus every generated report:

```bash
bash scripts/verify_project.sh
```

## Authority rules

- The 1985 paper is authoritative. It is not distributed here; cite and verify
  against the published article.
- Raw extraction and clean physics data stay separate.
- Every promoted numerical value retains provenance and confidence.

## License

MIT — see [LICENSE](LICENSE). Each package directory carries the same license
file. The 1985 article itself is not covered and is not distributed here.
