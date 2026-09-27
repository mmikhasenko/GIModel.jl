# Follow-up audit of the unfinished radial abstraction

Date: 2026-08-07

This audit follows the attached incident report describing an abandoned A2/A3
radial-wave migration. It distinguishes defects in the repository from correct
numerical work that happened to land in the same history.

## Verdict

The report's central claim was correct: the repository had implemented much of
the HO mathematics but discarded native HO coefficients at the solver boundary,
then made downstream code depend on a reconstructed mesh. The architecture was
therefore not complete even though planning notes said that it was.

That boundary is now repaired. Both solvers return `ChannelRadialSolution`, each
level is a native `RadialWave`, and consumers dispatch on radial/momentum
operations. No second solution, representation, or physical-state hierarchy was
introduced.

## What was retained

The numerical work in the later differentiation-support commit was not garbage.
The following changes are coherent and independently verified:

- parametric model components and exact constituent-mass storage;
- cache-key rounding separated from physics values;
- concrete FD solver/eigensolver dispatch;
- numeric-generic central Hamiltonian assembly; and
- the FD and fixed-beta HO Hellmann--Feynman probes.

The differentiation-support notes (since moved out of the repository into the
unreleased research workspace) are research and verification material, not
runtime code. Their historical audit documents intentionally retain the types
seen at their original snapshot; their implementation-status note states the
current native solution types.

## Damage and residue found

| Finding | Impact | Disposition |
| --- | --- | --- |
| `ChannelRadialSolution` emulated a three-tuple and virtual `.eigenvectors`/`.r` fields | Preserved the mesh-shaped API and allowed new code to keep depending on it | Removed completely |
| HO central solve reconstructed a reporting mesh before returning | Threw away the representation needed by the original algorithm | Now returns `OscillatorWave` coefficients natively |
| Spin, mixing, annihilation, decay, and report code unpacked `(values, vectors, r)` | Duplicated representation-specific code throughout the repository | Migrated to `radial_wave` and shared operations |
| Contact-state inactive cases returned empty tuple components | Could be interpreted as an algorithm fallback | Now return `nothing`; active cases return `ChannelRadialSolution` |
| Tests recreated solution objects and tested compatibility shims | Protected the wrong abstraction | Replaced with representation-independent invariants |
| `docs/code_architecture.md` described removed GIModel files and APIs that had moved to GIPaper | Misled the next worker about ownership and call flow | Rewritten from the live module |
| `docs/conventions.md` and `docs/formula_map.md` retained the same pre-split ownership model | Reintroduced dead names during follow-up work | Corrected to the GIModel/GIPaper split |
| `GIPaper/extraction/data_checks.py` read deleted `src/sector_comparison.jl` and searched for the old `compare` spelling | The implementation score path could fail or score zero for the wrong reason | Pointed at `GIPaper/src/comparison.jl` and current tests |
| The differentiation-support baseline probe still called `first(channel_solution(...))` and claimed HO reconstructed a mesh | The audit probe itself broke once the compatibility shim was removed | Updated and rerun successfully |
| Several paper scripts and one example used removed basis markers, loose solver keywords, tuple unpacking, or cached eigenvector matrices | Scripts outside unit-test coverage would rot silently | Migrated; all edited Julia scripts parse and both package test suites pass |
| `OscillatorSolver` still carried `ngrid`/`rmax` fields described as a reporting mesh | Kept a mesh-shaped concept in the HO calculation even after native waves existed | Fields, constructors, display text, and `with_mesh(::OscillatorSolver, ...)` removed; plots own explicit sampling grids |
| A mesh-bearing contact overload silently ignored its `r` argument for HO | Preserved the appearance that an HO solve could consume an FD operator grid | Mesh overload is FD-only; the representation-free overload dispatches directly to native `fixed_channel_solution` |
| `FDOriginP2Smearing` and label-inferred flavor coherence remained as deprecated compatibility paths | Kept a mesh-only approximation and stringly typed physics alive for old callers | Removed from the public API; `isoscalar_coherent` is mandatory |
| Pseudoscalar mixing used a second eigenvector-phase helper based on the largest component | Made phase conventions differ across otherwise identical `MixingBlock` calculations | Removed; the shared first-basis-component convention is used consistently |
| Annihilation selected one radial wave from a tensor-mixed physical state | Discarded the rest of the signed physical composition or failed once mixed `radial_wave` access became strict | It now coherently projects all matching `(L,S,J)` radial components; no representative wave or mesh fallback exists |
| `StateMixing` stored a scalar angle and maximum off-diagonal on every member | Duplicated shared-block data and produced `NaN` once mixing blocks correctly contained several radial levels | Removed both stored fields; maximum coupling is derived from the shared block and reports derive an explicitly documented pair projection |
| GIPaper comparison code reconstructed a supposed `2x2` block from two states | Indexed the wrong masses/components when the shared production block contained more radial levels | It now projects the requested two basis rows from the one shared `MixingResult` while retaining the full physical eigensystem |
| Adaptive HO exposed `Inf * 0 = NaN` in infinite-interval high-order oscillator integrals | Converged 64-72-state waves could fail only when later observables were evaluated | Native coordinate and momentum quadratures now stop beyond the analytic Gaussian tail; low-order results agree with the old integral to machine precision and high-order regressions cover both spaces |
| Current-status ledgers and generated report prose still described fixed 24-state, mesh-reconstructed, first-order-FD execution | They contradicted the implemented architecture and even contradicted a computed 2/4 sign result with “all four signs agree” | Report generators were corrected; the gap ledger now names only PA-17, PA-18, and FD-COMP, and generated provenance describes adaptive native fixed-sector solves |
| The adaptive controller silently raised the advertised initial `nbasis` to `8*nlevels` | Made solver configuration and provenance misleading | Removed; `nbasis` is the actual initial variational size and a request for more eigenlevels fails with instructions to raise it |

## Duplicate-entity audit

The repaired design reuses the existing nouns:

- method choice: `RadialSolver` (`FiniteDifferenceSolver`, `OscillatorSolver`);
- one radial level: `RadialWave` (`MeshWave`, `OscillatorWave`);
- one momentum representation: `MomentumWave` (`MeshMomentumWave`,
  `OscillatorMomentumWave`);
- a channel eigensolution: `ChannelRadialSolution`;
- cache and provenance: `SectorComputation`;
- mixing algebra: `MixingBlock` and `MixingResult`; and
- staged physical output: `CentralSpectrum`, `CorrectedSpectrum`, and
  `MixedSpectrum`.

The obsolete `RadialWaveOnUniformMesh`, `MockMomentumWave`,
`FDOriginP2Smearing`, label-inference, and HO mesh-setting surfaces were removed;
this pre-release repository does not carry a backward-compatibility layer.
`BasisState` and `FineStructureMultiplet` are intentionally not consolidated:
the former identifies a particular radial level (`n,L,S,J,label`), while the
latter identifies a fixed sector (`L,S,J`) before radial eigenlevels exist.
The remaining direct `.vectors` access in `audit_table_iii_mixings.jl` compares
the low-level `MixingResult` itself with the paper's printed eigenvectors; it is
not an observable consumer and does not construct a second physical state.

## Repository scheduling and follow-up at the time of this audit

The surviving responsibilities are intentionally separated:

- `docs/work_plan.md` was the sole executable plan and contained only
  completion of the Table VI photon-decay audit. That work is now complete;
  the file records its closure and verification.
- `docs/code_architecture.md` records the settled public contract and extension
  rules; regression tests pin the numerical invariants that must not recur.
- `docs/original_1985_algorithm_audit.md` records the completed paper-vs-code
  algorithm comparison, not task status.
- `docs/paper_manifest/*.toml` records per-paper-unit reproduction status.
- `scripts/verify_project.sh` is the repository integration gate, including
  report generation order and coverage.

Do not create another planning file for each implementation slice. Update the
one plan and the relevant manifest entry in the same change as code and tests.

## Deliberately untouched history

- Branch `table-v-reproduction` points at `05a86cf`, two commits behind `main`.
  It is historical, not part of the active implementation.
- Two stashes contain edits against the removed path
  `src/GIModel/spin_fine_structure.jl`. They were not applied or deleted because
  stash deletion is destructive and their intent cannot be inferred safely.
- The Table III forensics scripts are pre-migration. They are preserved outside
  the repository as historical evidence and are excluded from the live
  consumer contract.

## Verification performed

- GIModel `Pkg.test()`: passed.
- GIPaper `Pkg.test()`: passed.
- Paper manifest check: 103 units, 8 fragments, 0 broken links.
- Type-stability probe: concrete `ChannelRadialSolution{MeshWave}` and
  `ChannelRadialSolution{OscillatorWave}` returns.
- FD central-Hamiltonian derivative probe: passed.
- Fixed-beta HO derivative probe: passed.
- Audit baseline probe: repaired and passed.
- Every edited Julia paper/example script: syntax parsed successfully.
- FD-COMP: independent fixed-domain/fixed-spacing sweeps passed for eight
  six-level calibration sectors, eight mixed states, and three cross-state
  transition functionals.
- `scripts/verify_project.sh`: passed, including all report writers and the
  103-unit manifest anti-drift check.

The post-incident PA-17--PA-18 audit is also complete. It found four additional
physics defects that the attached worker report did not identify: the contact
term multiplied `alpha_s(r)` by a second Gaussian instead of taking the exact
Laplacian of the smeared Coulomb kernel; A16 duplicated a Coulomb derivative;
the tensor normalization was four times too large; and fitted spin bridge
parameters obscured those errors. The literal A15-A16 implementation now has no
bridge parameters, and final flavor-mixed observables consume the one shared
`MixedSpectrum` composition instead of reconstructing vectors locally.

The original spectrum-algorithm queue and its optional FD-COMP infrastructure
follow-up are therefore closed. FD is a numerically certified independent
comparator, explicitly outside the native-HO paper path. The separate whole-
paper manifest tracked partial Eq. (19) and Table VI decay units at the time of
this audit. The current manifest supersedes that historical status and records
complete 79-row Table VI coverage.
