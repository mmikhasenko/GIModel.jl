# Release editing plan

Updated 2026-09-26. This is the active checklist derived from the [release claims audit](release_claims_audit.md). The audit remains the historical evidence; this document tracks edits and closure. Numerical provenance is maintained in the [rate input ledger](../GIPaper/docs/rate_input_ledger.md).

## Agreed direction

GIModel is close to release. Its pass should remove obsolete material, clarify the shipped configuration and finish package presentation, without reopening the validated solver architecture. QuarkModelTransitions needs the main numerical-input cleanup. GIPaper needs a corresponding pass over reference data, table prescriptions, scripts and generated outputs.

The dependency direction is intentional: GIPaper → QuarkModelTransitions → GIModel, with GIPaper also directly using GIModel. Spectrum annihilation mixing remains in GIModel; decay observables belong to QuarkModelTransitions. The ordinary q/s/c/b spectrum Hamiltonian has **18 numerical physics inputs: 12 configurable values and six fixed Gaussian coupling values**. Optional flavor-annihilation prescriptions add inputs; transition observables add more. Numerical tolerances and quantum-number choices are recorded separately.

A folder pass finishes with a small, reviewable change and an explicit list of retained, edited, moved and deleted files. An old date or an unfashionable filename is not evidence that a file is unused. Conversely, redundant obsolete implementations do not need a new archive merely to avoid deleting them: Git preserves their history.

## Already completed — do not repeat

- [x] Verify GIModel's independence from both companion directories, including an isolated-copy spectrum smoke test.
- [x] Count and trace the 18 spectrum inputs.
- [x] Move radiative, leptonic, gluonic, two-photon and charge-radius observables and their tests to QuarkModelTransitions; retain generic wave operations in GIModel.
- [x] Update affected callers, example dependencies and equation-manifest paths.
- [x] Remove inactive `Lambda_MeV` from the active configuration; retain its Table II provenance in GIPaper.
- [x] Expose current electromagnetic and strong-decay defaults in `QuarkModelTransitions/src/observable_inputs.jl`.
- [x] Generate V/VI/VII input traces and preserve numerical report results through the move.
- [x] Run the three test suites, reference-map check and migrated annihilation example at the migration checkpoint.

These are completed changes in the working tree, not a claim that the whole release has passed its final gate. Concurrent native-transition edits require a fresh final checkpoint. In particular, the historical audit's “native matrix_element is not connected” finding is superseded: the current native route exists. Reassess its documented coverage against the current implementation rather than reopening the old missing-method task.

## Folder order

| Pass | Folder / scope | Main output | Completion evidence |
|---|---|---|---|
| G1 | Root metadata and `data/` | Clear installation and audited runtime preset | Preset/loader checks; minimal example outside repository cwd |
| G2 | `src/` | Disposition of every core file and exported diagnostic | Caller search plus core tests for actual removals/edits |
| G3 | `test/`, `examples/`, root `scripts/` | Current tests, examples and verification commands | Relevant suites/examples and link checks |
| G4 | Root `docs/` | One current navigation path, historical records clearly separated | No stale source/API links or competing active plans |
| T1 | `QuarkModelTransitions/src/` | Complete, explicit observable-input ownership | Input/consumer ledger; parameter propagation checks |
| T2 | `QuarkModelTransitions/test/`, `docs/`, README | Tested public contract matching current capabilities | Transition suite and runnable public examples |
| P1 | `GIPaper/data/` and `src/` | Canonical reference/policy records and adapters | Provenance, schema and assignment checks |
| P2 | `GIPaper/scripts/` | One forward-computation route per table; obsolete scripts removed | Rate/report regeneration and numeric comparisons |
| P3 | `GIPaper/docs/`, `extraction/` | Current generated outputs and reproducible provenance tools | Source hashes, coverage and regeneration instructions |
| R1 | Non-release material, `LearningTrack/` | Explicit release inclusion and historical status | Broken-link scan; no runtime dependency on historical work |
| R2 | Whole repository | Release candidate | Final clean-environment tests and full relevant gate |

Start with **G1**. Complete GIModel's short cleanup passes before the larger T/P passes. Do not combine all folders into one broad refactor.

## G1 — root metadata and active inputs

- [x] Review `Project.toml`, `Manifest.toml`, `.gitignore` and README. Declare the tested Julia compatibility range; document how companions are installed from this monorepo. Test an ordinary downstream environment, not just the developer checkout.
- [ ] Resolve release licensing with the repository owner if it has not already been chosen; add the selected license rather than inventing one.
- [x] Replace the active configuration's `provisional_raw_digitization` description and obsolete promotion-path note with reviewed provenance. Decide the final preset filename here and update `default_parameters_path()` and callers together; avoid keeping two identical active files.
- [x] Explicitly distinguish the runnable GI preset from GIPaper's clean transcription of Table II. Prevent loading the latter from silently choosing a different central Hamiltonian or disabling sandwiches. Define and test the intended loader error/preset behavior.
- [x] Define accepted configuration keys and missing-field policy. Reject misspelled/inactive physics keys with useful errors, while preserving documented programmatic diagnostic constructors. Do not silently turn an incomplete production file into a different model.
- [x] Make the 18-input inventory discoverable from the core documentation, including the six values in `src/constants.jl`. Do not add a second numerical source or turn fixed Gaussian values into fitted knobs solely for bookkeeping.
- [x] Make README start with GIModel installation and one runnable spectrum/wave example. Link companions, the report and learning material after the primary user path.

Done when the shipped preset has one clear identity, the documented example works from another directory, and inactive/incorrect configuration cannot masquerade as an effective input. Changing the preset name or metadata must not change masses.

## G2 — core source, file by file

- [x] Review every include/export in `src/GIModel.jl` against actual callers and the public documentation. Record supported public API, internal helpers and intentionally retained diagnostics.
- [x] Inspect the comparator families in `central_potential_dispatch.jl`, `radial_1d_coulomb_smear.jl`, `appendix_a_derivative_potential.jl`, `smearing_appendix_a.jl`, `appendix_a_status.jl` and their parameter-method exports. For each, identify a current test/audit/user purpose. Delete superseded methods and their plumbing only when that purpose has ended; retain independent numerical checks with a concise description. These are review candidates, not an approved blanket deletion list.
- [ ] Remove obsolete migration comments and stale algorithm descriptions from solver, spectrum and wave files; preserve useful equation references and numerical failure explanations. Review the compatibility name `mock_momentum_wave` as part of public API cleanup, without creating duplicate transform implementations.
- [x] Confirm no moved observable implementation/constants remain in the core. Keep `pseudoscalar_annihilation.jl` and `flavor_mixing.jl` while they implement mass-spectrum mixing; document their extra inputs outside the ordinary 18-number contract.
- [x] Resolve the native transition package's qualified derivative-overlap call through an intentional public wave-operation API. Recheck the current implementation before editing, because transition development is concurrent.
- [x] Clarify the scope of HO energy certificates and quadrature warnings. Preserve existing numerical behavior during editorial cleanup; any proposed change from warning to failure is a separate behavioral change with a demonstrated test case.

Done when each surviving file has a current role, deletions have no unresolved callers, and the independent FD/HO validation routes survive. No wholesale solver rewrite is required for this pass.

## G3/G4 — tests, examples, scripts and core documentation

- [ ] In `test/`, remove tests solely for deleted APIs alongside those APIs. Retain independent physics, normalization, convergence and boundary checks; keep observable tests under QuarkModelTransitions.
- [ ] In `examples/`, classify curated examples and scratch work. Review `played_with_model.jl` and example data individually; move ongoing exploration to research or delete superseded material after checking its unique content. Update examples to the new package ownership.
- [x] In root `scripts/`, keep a documented fast package check and the full reproduction gate. Remove dead command paths; ensure the full gate checks all three packages and the input-trace workflow at the appropriate stage.
- [ ] In `docs/`, reconcile README, architecture, formula map and discoverability pages. Correct claims that all solvers/defaults are interchangeable or that every observable uses only the spectrum inputs.
- [x] Reclassify `work_plan.md` as its existing Table VI completion record, not a competing release plan. Review reproduction/incident/equation audit documents for unique evidence; link historical records from one index, and remove redundant current-status summaries.
- [ ] Review generated HTML/dashboard artifacts separately from editable sources and establish what is shipped versus regenerated. Do not delete ignored local renders or separate repositories as part of a tracked-source sweep.

Done when the user can find one current API path, one release plan and clear verification commands. Check links and run only the examples/suites affected by that batch.

## T1 — QuarkModelTransitions numerical-input cleanup (main work)

- [ ] Inventory every numerical physics input in `observable_inputs.jl`, `annihilation_widths.jl`, `mock_meson_overlaps.jl`, `radiative_decays.jl`, `strong_decays.jl`, `pseudoscalar_emission.jl` and flavor-state conventions. Classify each as universal constant/unit convention, phenomenological input, calibration datum, exact algebra, derived quantity or numerical control.
- [ ] For each empirical/default input record symbol, value, units, owner, source/equation, precision, calibration channel (if any), consumers and override path. Use the current rate ledger as evidence, not a second set of constants.
- [ ] Consolidate alphaEM, GF, hbar-c and the nuclear-magneton mass convention in one authoritative transition location with provenance. Decide the precision policy explicitly; do not silently modernize values and call the resulting shifts cleanup.
- [ ] Make the 0.7/0.5 electromagnetic exponents, 0.2 charge-radius exponent and 0.40 GeV recoil scale propagate explicitly through the public observable APIs. Check higher-level wrappers: an overridable low-level kernel is insufficient if its public wrapper silently restores defaults. Keep variational HO beta distinct from phenomenological recoil beta.
- [ ] Separate the generic strong-decay operator/model parameters from the paper's calibration anchors (+12.4, −11.0) and `LeadingS0` prescription. Put reference calibration values/provenance in GIPaper, with an explicit compatibility decision for the frozen backend; avoid duplicated defaults in both packages.
- [ ] Keep native g/h and constituent routing explicit. Classify fixed eta/eta-prime flavor assumptions and phase conventions; do not replace documented algebra with per-channel numerical corrections.
- [ ] Inventory numerical cutoffs/tolerances separately. Remove stale controls only after proving they are unused; do not merge distinct integrations just because their default sample counts happen to match.

Done when every observable's numerical leaves have one owner and can be recovered from the inputs used for that call. Verify selected overrides through the full wrapper route, compare baseline rates, and retain calibration/validation separation. Cleanup should not change baseline physics without an explicit, explained correction.

## T2 — transition public API and documentation

- [ ] Reconcile README, phase documents and current native implementation with the agent developing it. The public native method and general orbital work have progressed beyond the original audit; identify actual remaining coverage gates from current tests.
- [ ] Remove only superseded development checkpoints that add no unique validation evidence. Keep frozen Table V reproduction visibly distinct from native solved-wave predictions.
- [ ] Provide runnable examples for the supported operator routes, with units, normalization, emitted-field role, closed channels and mass-source choices explicit.
- [ ] Extend boundary checks to the intended public GIModel wave API; verify that no paper-row imports or resonance-label corrections enter generic operator code.

Done when advertised capabilities match tested ones. Additional physics development is owned by its existing workstream, not automatically a release-cleanup prerequisite.

## P1 — GIPaper reference data and adapters

- [ ] Inventory `data/` by role: immutable raw evidence, canonical transcriptions, runtime assignments/calibration policy, pinned experimental inputs and historical records. For apparent duplicates, identify their consumers and provenance chain before deleting or promoting either file.
- [x] Keep historical Lambda and alpha-critical in Table II reference data, clearly marked as historical/derived relative to the active Gaussian runtime profile. Do not restore them as extra spectrum inputs.
- [x] Use the mass registry's explicit provenance approach for non-mass input policy: the Table V angles 34/33/−41 degrees, calibration amplitudes, Table VI +0.01 muN addition, isovector assignments, ideal-flavor prescriptions and legacy validation anchors.
- [ ] Replace label-switch policy in scripts with explicit validated records where appropriate. Preserve the distinction between paper-prescribed inputs, computed mixing and comparison-only targets; a named exception is not a universal operator rule.
- [ ] Review `src/` adapters for unused paths, private core access, duplicate parsing and fallback masses. Keep unknown/unavailable assignments explicit and maintain the allowed dependency on both computation packages.

Done when every external value in a rate trace resolves to a source record, and deleting an obsolete record cannot silently trigger a fallback.

## P2/P3 — GIPaper scripts, reports and extraction

- [ ] Classify every script as canonical report generator, independent validation, active investigation or obsolete/superseded work. Review overlapping spectrum plotters, mixing demos and old diagnostic audits by their actual outputs, not their filenames.
- [ ] Keep one forward prediction path for each canonical table. Preserve independent convergence checks; remove obsolete alternate calculations whose sole effect is reproducing a retired implementation.
- [x] Remove unused `NGRID`/`RMAX` in the HO Table VII harness after confirming no current consumer. Correct stale headers and “zero parameters” statements to name the active inputs and wave prescriptions.
- [ ] Make comparison thresholds, inverse implied-mass diagnostics and historical momenta visibly downstream of predictions. They must not feed calibration or modify the forward rate silently.
- [x] Add input-trace integrity checking: canonical row coverage, source fingerprints, declared input ownership and regenerated values must agree. Review literal numerical-control entries in `trace_rate_inputs.jl` so future changes cannot leave a plausible but stale inventory.
- [x] After concurrent source work settles, regenerate traces and reports together. Check numerical diffs, unavailable-row status, units and calibration labels. The current traces are snapshots, not guaranteed current after another agent changes their source files.
- [ ] In `docs/`, retain a single authoritative report per configuration and clear links from hand-written summaries. Retire stale duplicate reports; preserve explicit HO/FD comparisons.
- [ ] In `extraction/`, retain source evidence, checksums and necessary promotion/import tools. Remove superseded utilities only after the clean-data provenance no longer depends on them.

Done when every shipped result has a reproducible command and an input record, every unavailable row has a reason, and hand-written prose does not claim more than the current generator computes.

## R1/R2 — repository closure

- [x] Move non-release material out of the repository (2026-09-27): the archive, differentiation-support notes, decay research, later-paper study, paper source archive, executable paper remix and local scratch now live in a separate, git-ignored research workspace with its own repository. No released file links into it. The paper-layout table drivers moved into GIPaper (`GIPaper/scripts/paper_tables/`, outputs in `GIPaper/docs/paper_tables/`); Table III vector/tensor provenance now resolves inside `GIPaper/data/raw/`. The report stays a separate, git-ignored checkout.
- [ ] Treat `LearningTrack/` as a separate educational deliverable with its own build.
- [ ] Review package contents and documented installation from a fresh environment; resolve manifest policy and dependency compatibility for each package.
- [ ] Run all three package suites, the affected example smoke checks, reference/link checks, rate traces and the full relevant `scripts/verify_project.sh` gate on the final candidate. Attribute existing warnings; do not treat a test pass as proof of every numerical certificate.
- [ ] Record the release revision, validation commands/results and any intentionally unsupported capability in one final checklist. Apply version/tag/publish steps only as part of the actual release request.

## Paper handoff

The existing **Discrepancies** task owns `report/main.tex`. Requested edit: show the 18 numerical spectrum inputs in Table 1, with the three running-coupling weights in its caption beginning “In addition ...”; put the three Gaussian scales in the body, and connect the Fig. 2 discussion to Table 1 including its caption. Lambda and the derived alpha-critical value do not increase the count. The wording must distinguish the ordinary Hamiltonian inputs from any optional isoscalar or transition inputs in the plotted calculation. That task compiles/verifies the paper; this cleanup task does not concurrently edit it.

## Record after each folder pass

| Pass | Files kept / edited / moved / deleted, and reason | Validation | Remaining issue | Status |
|---|---|---|---|---|
| G1 | Metadata, strict preset loader, README, input inventory; filename retained for compatibility | 1,022 core tests; fresh `/tmp` downstream spectrum | Owner license choice | Implemented; licensing open |
| G2 | Generic derivative API, comparator status/docs; diagnostic and independent solver routes retained | Core suite, caller inventory below | No physics algorithm changes | Reviewed |
| G3/G4 | Fast/full gate commands and current documentation navigation | Shell syntax and diff whitespace checks | Full gate after parallel transition work | Implemented |
| P1 | Policy/provenance record, stricter mass registry, explicit script inputs | 2,652 GIPaper tests; fresh downstream load | Legacy validation edition unknown, explicitly recorded | Implemented |
| P2/P3 | Dead helper/control removal, trace integrity workflow | See final verification record below | Concurrent sources must settle for release-wide gate | Implemented |

Use this table for completion evidence rather than adding another status document. Review proposed deletions with their caller/provenance evidence in the folder's change; do not run blanket directory deletion commands.

## GIModel/GIPaper folder disposition — 2026-09-26

The transition package and active report/research investigations are owned by
separate pipelines. No transition implementation edits belong to this pass.

| Files / group | Disposition and current purpose |
|---|---|
| Root `Project.toml`, README, runtime preset | Edited Julia minimum, downstream installation, input inventory and reviewed preset metadata. Kept the single historical preset filename to preserve companion callers. |
| `src/parameters.jl`, `quark_mass_table.jl` | Strict file schema and finite-value checks; constructors retain diagnostic use. Optional complete isoscalar section remains separate. |
| `src/constants.jl`, `running_coupling.jl` | Retained fixed Gaussian input source; no Lambda evaluation or additional fit. |
| `src/model_objects.jl`, `momentum_waves.jl`, `quark.jl`, `meson.jl` | Retained generic states/waves/flavors; derivative overlap intentionally exported. `mock_momentum_wave` retained as a public compatibility transform name. |
| `src/central_potential_dispatch.jl`, `smearing_appendix_a.jl`, `radial_1d_coulomb_smear.jl`, `appendix_a_derivative_potential.jl`, `appendix_a_status.jl` | Retained diagnostic constructions: `test/central_potentials.jl` exercises typed dispatch and independent convolution identities. Corrected stale preset/precedence descriptions. |
| `src/radial_grid.jl`, `hamiltonian.jl`, `harmonic_oscillator_basis.jl`, `channel_solver.jl`, `fixed_channel_solver.jl`, `sector_solver.jl`, `solver_options.jl` | Retained FD/HO implementations, typed controls and native caching. Numerical behavior unchanged; energy certification is not an observable certificate. |
| `src/contact_hyperfine.jl`, `spin_fine_structure.jl`, `state_mixing.jl`, `spectrum.jl`, `pseudoscalar_annihilation.jl`, `flavor_mixing.jl` | Retained spectrum pipeline and mass mixing, including optional inputs beyond the ordinary 18. |
| `src/GIModel.jl`, `test/` | Core public boundary retained after operator migration; new tests cover strict parameter-file failures. |
| `examples/` | Curated solver/transition examples retained. `played_with_model.jl` retains unique interactive paper-comparison exploration and is explicitly uncurated; density demonstrations/data are separate examples, not runtime inputs. No blanket deletion of figures or local renders. |
| Root `scripts/` | Added fast GIModel/GIPaper test gate. Full gate regenerates V/VI/VII with traces together, then checks their integrity. Doc-graph generator retained. |
| Root `docs/` | Added model-input entry point; corrected formula-map typed dispatch and HO support. `work_plan.md` already identifies itself as the Table VI completion record and remains historical evidence. HTML is ignored/regenerated; source documents retained. |
| `GIPaper/data/clean/parameters.toml` | Relabeled as reference transcription. Lambda/derived alpha-critical preserved as paper evidence; full runtime loader rejects it. |
| `GIPaper/data/raw/`, `seed/`, reference CSVs | Retained extraction evidence, promotion inputs and canonical comparison datasets. They are not interchangeable runtime presets. |
| `GIPaper/data/mass_inputs/`, `src/mass_inputs.jl` | Pinned experimental inputs and archived historical comparisons retained. Loader rejects duplicate keys, invalid measurements, unresolved average components and unknown assignments. No fallback introduced. |
| `GIPaper/data/table_policy.toml`, `src/table_policy.jl` | Added source-bearing paper calibration/mixing/footnote records and explicitly comparison-only legacy anchors. Original edition of legacy validation widths remains unknown and is recorded as such. |
| Other `GIPaper/src/` adapters | Retained reference identity, meson mapping, Table III/VI/V parsing, comparison and report responsibilities. GIPaper's dependency on both computation packages is intentional. |
| `reproduce_table_v.jl` | Deleted unused parent-flavor/level mapping functions and unused mixing metadata; loads paper policy and explicitly passes calibration amplitudes. |
| `audit_table_vi_photon_decays.jl`, `audit_table_vii.jl` | Loads paper policy; removed unused VII FD mesh constants and inaccurate zero-parameter/unsupported-subtable headers. |
| `trace_rate_inputs.jl`, `check_input_traces.jl` | Same forward report computations; source-change guard, source/result fingerprints and canonical-row coverage gate. |
| Remaining canonical report scripts | `run_all_spectrum_checks`, `analyze_heavy_quarkonium`, `audit_table_iii_mixings`, `audit_mixing_angles`, `audit_realistic_factors`, `score_annihilation` retained for distinct reports. |
| Independent numerical audits | `audit_ho_convergence`, `audit_fd_convergence`, `audit_w6_ho_order`, `audit_appendix_a_ho_comparison`, `audit_nonmixing_contact`, `compare_table_vi_solvers` retained for independent comparisons. |
| Plot/demo scripts | Sector/ten-sector/full/simplified plotters, basis/offdiagonal plots and mixing demo retained: different outputs and presentation scopes, not duplicate prediction backends. |
| `investigate_*.jl` | Active research retained, including concurrent eta/tensor work. Not promoted into canonical rate prescriptions. |
| PDG fetch/import and `extraction/` tools | Retained explicit refresh/promotion/OCR/visual validation roles; offline package tests do not fetch or re-extract data. |

The remaining release-wide gate must wait for the parallel transition changes.
Licensing remains an owner decision. This pass does not claim that legacy
validation anchors have acquired missing historical provenance, or that optional
research deliverables have received their own release review.

### Verification for this pass

- Julia 1.11.6: GIModel 1,022/1,022 tests; GIPaper 2,652/2,652 tests.
- Existing 30-GeV coarse-grid and diffuse-HO quadrature warnings retained; no solver behavior changed.
- Fresh temporary downstream environment: standalone GIModel spectrum from `/tmp`, then local companion installation and reference/policy loading passed.
- Tables V/VI/VII regenerated through `trace_rate_inputs.jl`; tracked numerical reports and Table VI CSV unchanged (Table VI Markdown run date only).
- Input-trace integrity passed: 220 Table V, 79 Table VI and 61 Table VII canonical rows; source and result fingerprints verified.
- Reference manifest: 103 units, 8 fragments, zero broken links. Documentation graph: 118 core links, unchanged.
- Script shell syntax and `git diff --check` passed.
- The final all-package/full-reproduction release gate remains deferred until the separately owned transition changes settle. These results are a cleanup checkpoint, not a release tag.
