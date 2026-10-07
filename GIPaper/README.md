# GIPaper

Reference data and comparisons with Godfrey and Isgur (1985).
GIModel computes the spectra; `GIModel.QuarkModelTransitions` computes the
transition amplitudes. This package connects those calculations to paper rows.

- `data/`: paper spectra, tables, state assignments and comparison prescriptions.
  `later_godfrey_e1_widths.csv` holds E1 widths from Godfrey's later GI-model
  papers, the benchmark for unequal-mass E1 that the 1985 paper lacks.
  `mass_inputs/` records pinned experimental masses for decay kinematics;
  `provenance/` preserves the original transcription notes and source records.
- `src/`: Julia readers and paper-specific state matching.
- `checks/`: direct comparisons for spectra, mixing and Tables III, V, VI and VII.
- `reports/`: recorded comparison results, including known discrepancies.
- `test/runtests.jl`: a small reference-data and numerical smoke test.

Use Julia 1.11 or later. From the repository root:

```sh
julia --project=GIPaper -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia GIPaper/checks/run_all_spectrum_checks.jl
julia GIPaper/checks/audit_table_iii_mixings.jl
julia GIPaper/checks/audit_mixing_angles.jl
julia GIPaper/checks/reproduce_table_v.jl
julia GIPaper/checks/audit_table_vi_photon_decays.jl
julia GIPaper/checks/audit_heavy_light_e1.jl
julia GIPaper/checks/audit_table_vii.jl
```

Paper comparisons use `OscillatorSolver()`, following the original
investigation. This is independent of GIModel's faster finite-difference
default for general calculations. `compare_reference` defaults to HO; an
explicit FD solver or legacy mesh keywords select an independent cross-check.
Report headers must name the solver actually used. Existing FD reports remain
historical cross-checks until regenerated with HO.

Each comparison writes its results under `reports/`. These are comparisons,
not a claim that every paper value is reproduced: residuals, unavailable inputs
and historical conventions remain explicit. The full checks can take longer
than the small package test suite. Reports retained from earlier audits are
historical snapshots until regenerated.

Paper numbers and their qualifiers live in the data files. Table-specific
calibration and mixing choices live in `data/table_policy.toml`. Model
parameters come from `GIModel.default_parameters_path()`; the transcribed
Table II file records the paper, not a separate solver configuration.

Extraction, plotting, publication generation and exploratory studies are no
longer maintained here. Their previous versions remain in Git history.
