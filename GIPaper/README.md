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

## Ten-sector spectrum plots

Fetch the PDG reference table separately, then render either the complete or
the column-compressed spectrum:

```sh
python3 GIPaper/scripts/fetch_pdg_mesons.py
julia GIPaper/scripts/plot_ten_meson_spectra.jl
julia GIPaper/scripts/plot_ten_meson_spectra_simplified.jl
```

Both renderers write PNG and PDF versions under
`GIPaper/scripts/spectrum_plots/`. The complete plot contains the original
three $S$, two $P/D$, and leading $F/G$ families. The simplified companion is
restricted to $1S$--$3S$, $1P$--$2P$, $1D$, and $1F$ with $J\leq4$; it combines
the allowed $C$ partners for $1^+$ and $2^-$ in every self-conjugate sector.
The unique $0^{-+}$, $0^{++}$, $1^{--}$, and $2^{++}$ columns retain their
physical $C$ signs. The $J=3,4$ columns combine parity and charge conjugation
completely and are labeled by $J$ alone. In every displayed column,
experimental points above the highest retained calculated level are omitted.


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
- `GIPaper/docs/investigations/` contains reproducible mechanism studies. These
  link to the authoritative residual report, separate model operators or test
  hypotheses, and do not create an independent paper-comparison path. See its
  [study index](docs/investigations/README.md).
- Consumer publications may snapshot these reports and summarize their meaning,
  but must not implement an independent comparison path.
