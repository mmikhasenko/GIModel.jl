# Report refresh after the all-L contact correction

Date: 2026-09-27. Core baseline: `4c84a9e` (all-L smeared contact), plus the
existing uncommitted working-tree changes. No parameter fit or new core
physics change is part of this refresh.

## Dependency map

| Report item | Authoritative calculation / action |
|---|---|
| P-wave residual figure and causal narrative | Canonical seven-sector residual reports; `report/sources/pwave_residuals.jl` now reads values directly |
| Full and simplified ten-sector spectra | `plot_ten_meson_spectra.jl`, then simplified rendering from its fresh full CSV; fixed PDG snapshot |
| Vector separation densities | `examples/density_candidates.jl precompute`, then report renderer; includes changed D-wave partners |
| Density masses and radii | Caption refreshed directly from the new TOML caches |
| Tables V--VII and numerical input traces | `trace_rate_inputs.jl`, followed by `check_input_traces.jl` |
| Decay-channel census | Recalculate fixed-sector masses on 450/900 point grids, rerun thresholds and independent Python reconstruction |
| Threshold-ratio figure | Verify and rerender archived central-P/D, contact-resummed-S diagnostic; unaffected by the all-L correction by construction |
| Pipeline and conceptual illustration | No numerical model output; preserved |

The census and archived threshold study were available only in the local
`codex/a-class-universality` worktree (363fc54). Their report dependencies were
copied here without changing that worktree. Full threshold recomputation still
requires that branch's topology/coupling-coefficient API; it is not claimed as
a fresh full-spin-wave calculation. Rendering preserves the original numerical
source hashes and records the renderer hash separately.

## Reproduction

Run `bash scripts/refresh_report.sh` from the GIModel repository. It uses the
current working tree, including existing transition-code changes, and leaves
the PDG snapshot fixed. Julia 1.11.5 is used for the report gate and figures; the census was run
with Julia 1.11.6. Both use single-threaded BLAS and the existing package
environments. Python 3.11+, `rsvg-convert`, and local LaTeX are also required.
The standalone native editor cannot resolve this multi-file report's figures
and bibliography; use the documented local `latexmk -pdf` build.

## Interpretation

The large light P-wave discrepancy is now a demonstrated missing-contact bug,
not an unresolved cancellation diagnosis. At this refresh checkpoint the report retained
the remaining spin-orbit and tensor-composition discrepancies (the same-J question was
subsequently closed by the literature audit) and distinguished convergence,
implementation checks, and reproduction of GI's published amplitudes. Table VI
coverage is 79/79; this is not a claim that every amplitude agrees.

## Execution results

- Decay census: 205 -> 201 isospin-plus-charge-conjugation groups, split as
  27 isovector, 40 isoscalar, 45 strange, 49 charmed, 40 bottom-flavored.
  Both grids give the same open/closed decisions. Independent Python
  reconstruction matches all 1680 final-channel keys, cuts and groupings.
- Archived threshold diagnostic: all 51,840 raw pieces, 9,720 comparisons,
  3,240 certificates and 7,290 plot contributions pass consistency checks.

- Package tests: GIModel 1041/1041; QuarkModelTransitions 598/598;
  GIPaper 2652/2652.
- Recomputed full ten-sector spectrum: 330 rows; both full and simplified
  PNG/PDF figures refreshed from the new table.

- All three density caches recomputed; each model-source fingerprint matches
  this checkout. Maximum angular-normalization error is 2.79e-8, below 1e-6.
  Masses and radii are unchanged at the publication's three-decimal precision.

- Density publication image rerendered from the fresh caches.

- Native HO convergence and high-resolution FD comparator audits regenerated
  successfully; FD certifies 8 calibration sectors and 8 mixed states.

- Tables V--VII and common/row input traces regenerated from current sources.
  Table VI still covers 79/79 rows. Its E1 median magnitude ratio is 1.009
  (previous report 1.028); the two magnetic moments quoted in the report
  remain +0.688 and -1.200 at displayed precision.

- Input traces pass source-fingerprint, ownership and row-coverage checks:
  Table V 220 rows, Table VI 79, Table VII 61.
- Independent FD Table VI run and the existing excited-eta/tensor-photon
  trace-replay diagnostics regenerated as well.
- Appendix-A HO comparison regenerated: 72 channels, mean absolute
  HO/FD mass difference 0.2 MeV and maximum 1.2 MeV.

- Refreshed Table VI HO/FD median relative differences: M1 0.099%,
  E1 0.047%, M2 0.094%.

- Complete `scripts/verify_project.sh` finished successfully, including contact,
  mixing, realistic-factor, heavy-quarkonium, spectrum and annihilation reports.
  Manifest check: 103 units, 8 fragments, zero broken links.
- Local `latexmk -pdf` build succeeds (11 pages), with no LaTeX warnings,
  undefined references or overfull boxes. All pages were visually inspected,
  with changed spectrum, density and status pages rechecked after regeneration.
- Both repository diffs pass `git diff --check`. Computed-text refresh is
  idempotent. Merging is deferred until the follow-up investigation.

The known HO quadrature warning is still emitted by the numerical pipeline;
its basis/observable gates pass, but that does not certify the requested
quadrature tolerance for every internal matrix. The pre-existing package
manifest compatibility warnings were left visible; no dependency upgrade or
parameter adjustment was made to suppress them.

Subsequent closure: the large same-J caption mismatch is a historical reference
discrepancy, corroborated by later GI-model calculations; see the
[target ledger](../../data/mixing_angle_targets.md). The ψ(3.82) tensor
admixtures remain a separate open issue. Its
threshold figure remains the verified restricted-wavefunction diagnostic,
not a newly computed full-spin-wave result.
