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

Run a reproduction script directly; plotting audits activate
`GIPaper/scripts/Project.toml` themselves:

```sh
julia GIPaper/scripts/reproduce_table_v.jl
```


The complete radiative-decay audit needs only the GIPaper environment:

```sh
julia GIPaper/scripts/audit_table_vi_photon_decays.jl
GI_TABLE_VI_SOLVER=fd julia GIPaper/scripts/audit_table_vi_photon_decays.jl
julia GIPaper/scripts/compare_table_vi_solvers.jl
```

`load_table_vi()` provides canonical transition identities and preserves the
printed prediction qualifiers. `load_table_vi_states()` supplies explicit
spectroscopic assignments and kinematic mass provenance. The audit writes one
computed result per canonical row; no reference number is duplicated in its
transition implementation.

## Ownership and single-source rules

- `GIModel/src/` owns reusable numerical physics. It must not import paper rows,
  historical assignments, residual thresholds, or publication prose.
- `GIPaper/data/` owns reference values and external-input provenance.
- `GIPaper/src/` owns reusable adapters from those references to `GIModel` APIs.
- `GIPaper/scripts/` owns comparison policy, residual classification, and the
  detailed diagnosis attached to deviations.
- `GIPaper/docs/residual_reports/` is generated output and is the authoritative
  numerical audit record. Hand-written ledgers link to it instead of copying
  row counts or per-row diagnoses.
- Consumer publications may snapshot these reports and summarize their meaning,
  but must not implement an independent comparison path.
