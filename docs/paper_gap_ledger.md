# Paper Gap Ledger

This is the concise physics-facing list of what is still missing relative to
the 1985 Godfrey–Isgur calculation. Task status and dependency order live in
[`paper_algorithm_work_plan.md`](paper_algorithm_work_plan.md); equation/table
status lives in `paper_manifest/*.toml`. This file does not duplicate those
queues.

For the full three-stage algorithm audit and object/method map, see
[`original_1985_algorithm_audit.md`](original_1985_algorithm_audit.md).

## Implemented algorithm

- Both numerical representations return the same `ChannelRadialSolution`
  contract with native waves: `MeshWave` for FD and `OscillatorWave` for HO.
- `fixed_channel_solution` assembles central, contact, symmetric spin-orbit,
  and diagonal tensor matrices before diagonalizing each fixed `(L,S,J)`
  sector. The HO implementation never constructs an operator mesh.
- The HO controller continuously refines beta and enlarges the basis until all
  requested energies pass 0.1 MeV on two successive refinements. Its achieved
  `OscillatorConvergence` certificate is stored on the channel solution;
  endpoint railing and basis-cap exhaustion are errors.
- Stage 2 builds complete tensor and unequal-mass antisymmetric spin-orbit
  blocks from all compatible requested radial states.
- Stage 3 builds general Eq. (16) and P1/P2 isoscalar annihilation blocks.
  `physical_components` recursively exposes the final signed, flavor-tagged,
  solver-native wave composition.
- FD remains an independent implementation/comparator. It is not part of the
  original paper algorithm and is not a hidden fallback for HO.

## Remaining gaps

### PA-17 — final-state observable consumers

Some flavor-sensitive electromagnetic and decay reports still assemble their
own flavor/annihilation coefficient vectors. They must consume the final
`physical_components` composition directly so masses, coefficients, phases,
and waves all come from the same physical state.

Acceptance: no observable maintains a parallel flavor eigenvector or selects a
representative wave from a mixed state; ground and radial isoscalar rows are
covered by end-to-end tests.

### PA-18 — literal paper-strength certification

The active parameters still contain `k_spin_orbit = 0.48` and
`k_tensor = 0.42`. These are repository bridge factors, not parameters in the
paper. Native HO matrices and staging are implemented, but the headline paper
mode cannot be called literal until those factors are removed or isolated in an
explicitly named comparator mode and the resulting formula normalization is
validated.

Acceptance: the native-HO spectrum and observable gate passes with paper inputs,
adaptive convergence certificates, final physical-state consumers, and no
unexplained spin bridge scale.

### FD-COMP — independent FD convergence report

FD should be validated on its own terms by varying `ngrid` and `rmax`, then
compared with the converged HO result as diagnostic evidence. FD agreement is
useful, but it is not an acceptance condition for PA-18.

Acceptance: internally converged FD masses and wave-sensitive observables are
reported separately from the paper certification.

## Explicitly not missing

The following were gaps in older reports but are no longer open work: native HO
waves, Appendix-A position/momentum matrices, contact/spin-orbit/tensor fixed-
sector assembly, multi-radial tensor and antisymmetric spin-orbit mixing,
isoscalar flavor/radial annihilation, final state composition, and automatic HO
basis/beta convergence. Historical descriptions of mesh-reconstructed HO waves
or first-order-only production are not current architecture.
