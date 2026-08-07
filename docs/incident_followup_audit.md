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

`diff_support/` is research and verification material, not runtime code. Its
historical audit documents intentionally retain the types seen at their original
snapshot; `diff_support/implementation_status.md` now states the current native
solution types.

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
| `GIPaper/scripts/data_checks.py` read deleted `src/sector_comparison.jl` and searched for the old `compare` spelling | The implementation score path could fail or score zero for the wrong reason | Pointed at `GIPaper/src/comparison.jl` and current tests |
| `diff_support/probes/audit_baseline.jl` still called `first(channel_solution(...))` and claimed HO reconstructed a mesh | The audit probe itself broke once the compatibility shim was removed | Updated and rerun successfully |
| Several paper scripts and one example used removed basis markers, loose solver keywords, tuple unpacking, or cached eigenvector matrices | Scripts outside unit-test coverage would rot silently | Migrated; all edited Julia scripts parse and both package test suites pass |

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

The obsolete `RadialWaveOnUniformMesh` and `MockMomentumWave` aliases were
removed; this pre-release repository does not carry a backward-compatibility
surface. The remaining semantic overlap between `FineStructureMultiplet` and the `L,S,J` subset of
`BasisState` is recorded for later consolidation; changing it is not required
to finish the radial representation boundary.

## Repository scheduling and follow-up

The work is intentionally split among existing repository mechanisms:

- `docs/paper_algorithm_work_plan.md` is the dependency-ordered executable
  queue. PA-01--PA-07 are complete; PA-08 is the next available physics unit.
- `docs/engineering_work_plan.md` records cross-cutting code invariants and
  mistakes that must not recur.
- `docs/original_1985_algorithm_audit.md` records the paper-vs-code algorithm
  gap, not task status.
- `docs/paper_manifest/*.toml` records per-paper-unit reproduction status.
- `scripts/verify_project.sh` is the repository integration gate, including
  report generation order and coverage.

Do not create another planning file for each implementation slice. Update the
dependency board and the relevant manifest entry in the same change as code and
tests.

## Deliberately untouched history

- Branch `table-v-reproduction` points at `05a86cf`, two commits behind `main`.
  It is historical, not part of the active implementation.
- Two stashes contain edits against the removed path
  `src/GIModel/spin_fine_structure.jl`. They were not applied or deleted because
  stash deletion is destructive and their intent cannot be inferred safely.
- `archive/table_iii_forensics/` contains pre-migration scripts. It is preserved
  as historical evidence and is excluded from the live consumer contract.

## Verification performed

- GIModel `Pkg.test()`: passed.
- GIPaper `Pkg.test()`: passed.
- Type-stability probe: concrete `ChannelRadialSolution{MeshWave}` and
  `ChannelRadialSolution{OscillatorWave}` returns.
- FD central-Hamiltonian derivative probe: passed.
- Fixed-beta HO derivative probe: passed.
- Audit baseline probe: repaired and passed.
- Every edited Julia paper/example script: syntax parsed successfully.
- `scripts/verify_project.sh`: passed, including all report writers and the
  103-unit manifest anti-drift check. The final phase-dispatch consolidation was
  then rechecked by both package suites, the Table III audit, and manifest check.

The remaining PA-08--PA-18 work is not cleanup from this incident. It is the
documented implementation of the full 1985 paper-order spin Hamiltonians,
mixed-state components, and final headline-path certification.
