# Paper Coverage Ledger

This is the concise physics-facing record of coverage relative to the 1985
Godfrey–Isgur calculation. Table VI now has complete 79-row coverage, recorded in
[`work_plan.md`](work_plan.md); equation/table status lives in
`paper_manifest/*.toml`. This file does not duplicate that plan.

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
and all 8 clean two-photon signs. Excited isoscalar-pseudoscalar magnitude residuals
remain visible model discrepancies, not implementation fallbacks.

### FD-COMP — independent FD convergence report (complete)

FD is validated on its own terms by varying `ngrid` at fixed `rmax` and varying
`rmax` at fixed spacing. The six-level q/s/c/b `1S0` and `3P2` calibration gates
masses, RMS radii, and smeared-origin factors; separate blocks cover unequal-
mass P mixing, direct P/F tensor mixing, and allowed/cancellation-sensitive
transition functionals. The certified precision profile is `(ngrid, rmax) =
(2400, 32.0)`; `(450, 24.0)` remains the faster historical report setting.

The audit also found and fixed an FD observable-adapter defect: `pmax = π/h`
grew under coordinate refinement while the momentum sample count stayed fixed,
silently coarsening the momentum quadrature. The transform now caps the
independently checked physical range at 60 GeV.

Results and the complete sweeps are in
[`fd_comparator_convergence.md`](../GIPaper/docs/residual_reports/fd_comparator_convergence.md).
HO agreement remains diagnostic evidence, not an acceptance condition for
PA-18 or a dependency of the original 1985 algorithm.

## Closed: same-J spin-orbit mixing investigation

The missing L>0 smeared contact interaction is fixed, and solving/reporting now
share its operator definition. The remaining large mismatch with the GI85
caption angles is a historical reference discrepancy: later GI-model
calculations corroborate the corrected result. The historical cause is unknown;
matching those captions is no longer an open implementation task. The original
1985 reference policy and all 13 caption targets remain unchanged.

The [target ledger](../GIPaper/data/mixing_angle_targets.md) adds six GK91
angles alongside the originals. The [generated comparison](../GIPaper/docs/residual_reports/mixing_angles.md)
retains numerical residuals and projected singlet probabilities. Small residuals
against later tables are optional precision follow-up. The ψ(3.82) S–D
admixtures, annihilation and photon-decay residuals remain separate open issues.

## Outside the spectrum-algorithm queue

Table VI is complete as an implementation/accounting audit: all 79 rows use
native fixed-channel waves and canonical targets. Excited eta and deeply
cancelled transitions retain explicit numerical residuals. The old
`Upsilon'' -> eta_b gamma` sign discrepancy was an encoding error;
the printed -0.004 agrees with the calculation.
Applying Eq. (19) directly to the calculated physical wavefunctions goes
beyond the paper's numerical SU(6)/single-beta SHO treatment. It is now a
separate native `QuarkModelTransitions` workflow, not a reproduction gap or a
fallback used by the frozen Table-IV/V reports.

## Explicitly not missing

The following were gaps in older reports but are no longer open work: native HO
waves, Appendix-A position/momentum matrices, contact/spin-orbit/tensor fixed-
sector assembly, multi-radial tensor and antisymmetric spin-orbit mixing,
isoscalar flavor/radial annihilation, final state composition, and automatic HO
basis/beta convergence. Historical descriptions of mesh-reconstructed HO waves
or first-order-only production are not current architecture.
