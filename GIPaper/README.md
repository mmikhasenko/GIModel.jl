# GIPaper

GIPaper is the paper-reproduction and comparison layer for GIModel.

- `GIModel` owns model parameters, physics types, numerical solvers, spectra,
  amplitudes, widths, and other reusable computations.
- `GIPaper` owns digitized 1985-paper data, paper-row adapters, reference-state
  matching, residual reports, and reproduction audits.
- `GIPaper/src/` depends only on GIModel's exported API. Plotting and audit-only
  dependencies live in the separate `GIPaper/scripts/` environment.
- `GIPaper/extraction/` contains the isolated Python OCR/data-provenance tools;
  neither Julia package depends on them.

Run the package tests with:

```sh
julia --project=GIPaper -e 'using Pkg; Pkg.test()'
```

Run a reproduction script directly; each script activates
`GIPaper/scripts/Project.toml` itself:

```sh
julia GIPaper/scripts/reproduce_table_v.jl
```
