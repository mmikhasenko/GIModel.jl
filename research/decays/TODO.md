# Strong-decay and resonance TODO

This is the actionable companion to
[quark_model_resonances.md](quark_model_resonances.md). Items remain open until
their implementation, tests, and documentation are complete.

The proposed public API and dispatch boundaries are specified in
[transition_matrix_element_api.md](transition_matrix_element_api.md).
The review and implementation sequence is specified in
[matrix_element_api_implementation_plan.md](matrix_element_api_implementation_plan.md).

## Immediate audit: what `matrix_element` means

- [x] Correct the units documented for
  `matrix_element(::StrongDecayAmplitude)`. In the current factorization,
  `coefficient * reduced` is dimensionless; the analytic `spatial_overlap`
  carries `MeV^(1/2)`, and only `total` has amplitude units.
- [x] Audit the `StrongDecayAmplitude` docstring, which currently describes all
  decomposition fields as `MeV^(1/2)` even though the factors have different
  dimensions.
- [x] Rename the one-argument method to
  `reduced_matrix_element(::StrongDecayAmplitude)` in the `0.3.0` breaking
  commit, with no deprecated alias; reserve three-argument `matrix_element` for
  a complete transition.
- [x] Add dimensional and factorization tests asserting
  `total == reduced_matrix_element(a) * spatial_overlap` and
  `decay_width(a) == abs2(total)`.
- [x] Replace `convention::Symbol` with `TableIVPolynomial()` and
  `LeadingS0()` dispatch types.
- [x] Make the tutorial distinguish the reduced interaction/spin-flavor factor
  from the complete Table V amplitude and from solver-native helicity
  amplitudes and their derived `⟨BC;k,LS|T|A⟩` projections.

## Existing infrastructure to reuse

- [ ] Inventory the exact public building blocks and record their conventions:
  `radial_overlap`, `momentum_overlap`, `momentum_wave`,
  `physical_state_amplitude`, and `physical_transition_amplitude`.
- [ ] Treat the existing Table IV/V implementation as the analytic SHO
  reference, not as code to replace.
- [ ] Treat the Appendix-D mock-meson overlaps as examples of applying
  solver-native HO and FD waves to transition operators.

## WP1: solver-native strong-decay amplitudes

- [ ] Implement GI Eq. (19) pseudoscalar emission as the first native operator;
  keep `^3P0` behind its later three-route spatial-integral prototype.
- [ ] Define the complete off-shell helicity vector and derive every allowed
  `M[A → BC; k, L, S] = ⟨BC;k,LS|T|A⟩` by a tested Jacob--Wick projection.
- [ ] Make `TwoMesonChannel` contain daughters only; derive rather than request
  `(L,S)` and identical-particle exchange symmetry.
- [ ] Factor every pure-basis result into exact algebraic coefficients times a
  shared, wave-keyed basis of spatial integrals.
- [x] Introduce eagerly resolved `PhysicalState` values; no matrix-element
  layer may retain or inspect a `Spectrum`.
- [x] Support legacy quasi-two-body labels through a mass-only `ReferenceState`
  that native operators reject explicitly.
- [x] Specify threshold behavior: closed on-shell widths vanish, while a
  closed-channel vertex requires explicit complex `CMKinematics` and is never
  silently zeroed.
- [ ] Evaluate the spatial integrals with the actual initial and daughter GI
  wave functions and combine mixed physical states coherently.
- [ ] Add a three-state coherent composer for one parent and two daughters. It
  must conjugate final-state coefficients and retain every component triple in
  the amplitude decomposition.
- [ ] Use one physical mass per external state for kinematics; never recompute
  phase space from individual component masses.
- [ ] Move identical-daughter normalization and exchange phases from row data
  into the typed `TwoMesonChannel` construction.
- [ ] Add phase-invariance and interference tests for mixed initial and final
  states.
- [ ] Include a `mixing_only` regression row whose pure-basis coefficient is
  zero but physical amplitude is nonzero.
- [ ] Recover the existing equal-`β` analytic Table IV/V amplitudes in the
  appropriate limit.
- [ ] Compare exact waves with equal-`β` and rms-matched SHO surrogates to
  isolate node and scale effects.
- [ ] Separate calibration channels from validation channels.

## Later work

- [ ] Keep continuum dressing explicitly labelled as an optional post-GI
  extension; it must not change the 1985 reproduction path or its parameters.
- [ ] Use the momentum-dependent vertices in self-energy integrals.
- [ ] Keep loop-induced, energy-dependent mixing in the propagator matrix; do
  not encode it as an ordinary `MixedState` or apply it twice.
- [ ] Introduce a separate complex `ResonancePole`/`PoleState` result with
  explicit left/right-vector and residue normalization conventions.
- [ ] Dress and mix the positive-parity `D_s` states through `DK` and `D*K`.
- [ ] Continue amplitudes to complex poles and construct unitary real-axis
  scattering amplitudes.
- [ ] Compare the coupled-channel Hamiltonian with finite-volume lattice
  spectra.
- [ ] Stress-test the framework on overlapping `2S`/`1D` light-vector
  resonances.
