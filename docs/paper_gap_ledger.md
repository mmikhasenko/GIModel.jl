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

## Completed final integration

### PA-17 — final-state observable consumers

Complete. Flavor-sensitive electromagnetic reports now call
`physical_state_amplitude` or `physical_transition_amplitude` on the final
`Spectrum`. Their kernels see signed `physical_components`; no report owns a
second annihilation eigenvector, state-dependent phase, or representative wave.

### PA-18 — literal paper-strength certification

Complete. The non-paper `k_spin_orbit` and `k_tensor` parameters and call
arguments were deleted. During the formula audit, three compensating errors
were exposed and corrected directly from A15-A16:

- the contact density is the Laplacian of the smeared Coulomb potential,
  `sum(alpha_k * delta_tau_k)`, not `alpha_s(r) * delta_sigma`;
- A16 contains the self-pair scalar-confinement derivatives, not a second
  Coulomb derivative; and
- the conventional Pauli `S12` angular matrix element multiplies the A15
  spin bracket with `1/12`, not `1/3`.

The A15 vector self/pair masses and A16 scalar self masses are assembled
separately before their fixed angular contractions. The same literal kernels
serve fixed-sector and later antisymmetric/tensor mixing elements.

The native-HO convergence gate passes six radial levels in q/s/c/b `1S0` and
`3P2` sectors without quadrature warnings. The end-to-end Table-VII gate keeps
all 16 gluonic signs with median magnitude ratio 1.05, all 26 leptonic signs,
and all 8 clean two-photon signs. The two excited isoscalar-pseudoscalar signs
remain a visible model discrepancy, not an implementation fallback.

### FD-COMP — independent FD convergence report

FD should be validated on its own terms by varying `ngrid` and `rmax`, then
compared with the converged HO result as diagnostic evidence. FD agreement is
useful, but it is not an acceptance condition for PA-18.

Acceptance: internally converged FD masses and wave-sensitive observables are
reported separately from the paper certification.

This is the next spectrum-infrastructure item. It is comparator validation, not
unfinished work in the original 1985 HO algorithm.

## Outside the spectrum-algorithm queue

The paper manifest still marks two decay-side units partial: Eq. (19)'s deeper
wavefunction treatment beyond the SU(6)/single-beta SHO limit, and completion
of the remaining Table VI rows (including the cancellation-sensitive
`Upsilon'' -> eta_b gamma` sign). They block a literal “whole paper reproduced”
claim, but they are not missing stages in the HO meson-spectrum algorithm.

## Explicitly not missing

The following were gaps in older reports but are no longer open work: native HO
waves, Appendix-A position/momentum matrices, contact/spin-orbit/tensor fixed-
sector assembly, multi-radial tensor and antisymmetric spin-orbit mixing,
isoscalar flavor/radial annihilation, final state composition, and automatic HO
basis/beta convergence. Historical descriptions of mesh-reconstructed HO waves
or first-order-only production are not current architecture.
