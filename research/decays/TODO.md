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

- [x] Implement GI Eq. (19) pseudoscalar emission as the first native operator;
  keep `^3P0` behind its later three-route spatial-integral prototype.
- [ ] Define the complete off-shell helicity vector and derive every allowed
  `M[A → BC; k, L, S] = ⟨BC;k,LS|T|A⟩` by a tested Jacob--Wick projection.
- [x] Generate the Appendix-C vector-pseudoscalar projection from Clebsch--Gordan
  coefficients and reproduce every sign and magnitude in Table XI through
  parent `J=5`.
- [x] Implement normalized sparse Appendix-B flavor states, exact spectator
  selection, and separate quark/antiquark flavor contractions.
- [x] Implement the coupled two-spin-1/2 states and spherical Pauli matrix
  elements needed by the elementary-emission operator.
- [x] Make `TwoMesonChannel` contain daughters only; derive rather than request
  `(L,S)` and identical-particle exchange symmetry.
- [x] Factor the Phase-III elementary-emission target set into algebraic
  coefficients times a
  shared, wave-keyed basis of spatial integrals. The Eq. (19) term record now
  separates spin, flavor, topology, spherical contraction, and orbital label;
  full spectroscopic recoupling produces helicity and partial-wave coefficient
  matrices over the same integral columns. Numerical integral values remain
  Phase IV.
- [x] Introduce eagerly resolved `PhysicalState` values; no matrix-element
  layer may retain or inspect a `Spectrum`.
- [x] Support legacy quasi-two-body labels through a mass-only `ReferenceState`
  that native operators reject explicitly.
- [x] Specify threshold behavior: closed on-shell widths vanish, while a
  closed-channel vertex requires explicit complex `CMKinematics` and is never
  silently zeroed.
- [x] Evaluate the spatial integrals with the actual initial and daughter GI
  wave functions and combine mixed physical states coherently.
- [x] Add a three-state coherent composer for one parent and two daughters. It
  must conjugate final-state coefficients and retain every component triple in
  the amplitude decomposition.
- [x] Use one physical mass per external state for kinematics; never recompute
  phase space from individual component masses.
- [x] Move identical-daughter normalization and exchange phases from row data
  into rules derived from the typed `TwoMesonChannel` and partial wave.
- [x] Add phase-invariance and interference tests for mixed initial and final
  states.
- [x] Include a `mixing_only` regression fixture whose pure-basis coefficient
  is zero but physical amplitude is nonzero.
- [ ] Recover the existing equal-`β` analytic Table IV/V amplitudes in the
  appropriate limit.
- [ ] Compare exact waves with equal-`β` and rms-matched SHO surrogates to
  isolate node and scale effects.
- [ ] Separate calibration channels from validation channels.

### Investigation: what the SU(6)/SHO decay approximation costs

Paper-reading note (2026-09-26): Eq. (19) is not absent from the original
analysis. It motivates the elementary pseudoscalar-emission operator and is
reduced with SU(6), pure spectroscopic states, and common harmonic-oscillator
waves to obtain the Table IV relations. The numerical Table V analysis then
explicitly *foregoes* enforcing all reduced amplitudes in terms of the same
`g,h`: it treats the Table IV classes phenomenologically and adopts
`A' ≃ A'' ≃ A_c ≃ A0 ≃ A` and `S ≃ D ≃ P ≃ S_c` for a two-strength fit.
The native `PseudoscalarEmission` route is therefore a beyond-paper test: it
applies Eq. (19) directly to the resolved GI wavefunctions and mixed states.

The calibration is global, but it should not be described as one fitted pair
of `g,h` used numerically in every sector. Table V uses exactly two input
decays: `rho -> pi pi` fixes the leading `A` strength and
`B -> [omega pi]_S` fixes the leading `S` strength. All other rows are then
predictions of the two-strength approximation, together with the assumed
class ratios, masses, mixing conventions, and phase-space prescription.
Light and strange channels share this calibration. Charm is not refitted:
the table retains explicit `[A_c/A]` or `[S_c/S]` factors and uses the
charm-specific recoil/form factor, while the numerical estimate adopts the
same approximate reduced strength. Table V contains no bottom-meson decay
block, however, so universality through the bottom sector is an extrapolation,
not something tested by this analysis. Its footnote `f` means that finite
widths of final-state particles were included; it is not a fit marker.

The decay-model `β=0.40 GeV` is not the variational `β` used to represent and
solve the GI Hamiltonian in an oscillator basis. The latter is a numerical
basis coordinate whose residual influence must disappear with basis
convergence. The former is a new physical approximation introduced only after
the calculated waves are discarded in favor of one common SU(6)/SHO family.
The paper chooses `0.40 GeV` as a compromise among SHO scales fitted separately
to the exact pion, rho, A2, and g coordinate- and momentum-space sizes (quoted
fits span `0.33--0.76 GeV`). Consequently the frozen Table V reproduction may
use `β=0.40 GeV`, but the native-wave amplitude must contain no independently
refitted decay `β`.

For the normalized three-dimensional SHO convention of Appendix B,
`ψ ∝ exp(-β²r²/2)`, define `C_nL = 2n + L + 3/2`. Then
`<r²>_SHO = C_nL/β²` and `<p²>_SHO = C_nL β²`, with `r` the relative
quark-antiquark coordinate. Fitting one realistic wave separately gives
`β_r = sqrt(C_nL/<r²>)` and `β_p = sqrt(<p²>/C_nL)`. Here `β,p` are in GeV,
`r` in `GeV^-1`, and conversion to `fm²` multiplies `<r²>` by
`(ℏc)² ≃ (0.19733 GeV fm)²`. For `1S`, `1P`, and `1D`, respectively,
`C_nL = 3/2, 5/2, 7/2`. The two fitted scales agree only when the realistic
wave has the corresponding SHO shape; their mismatch is itself a measure of
the approximation.

- [ ] In the pure, equal-mass, equal-`β` SHO limit, derive with the native
  contraction machinery every representative Table IV class and verify
  `A=(g+h/4)β`, `A'=(g-h/4)β`, `A''=(g+h/8)β`, `A0=gβ`, and the printed
  `S,D,P` momentum polynomials, including all signs and Table V spin-flavor
  coefficients.
- [ ] Optional reference validation: recover the exact signed `A qtilde^4`
  coefficients of the `1^3F_4` Table V block in the common-SHO limit, including
  `-1/sqrt(241920)`. This is not a prerequisite for native predictions: the
  production path evaluates Eq. (19) directly with the resolved wavefunctions
  and performs the angular projection numerically, without factoring the result
  into the paper's compact SHO coefficient times `A qtilde^4`. Keep the
  digitized coefficients only in the frozen Table V reproduction backend.
- [ ] Keep separate (a) fractions inside the reduced classes, such as `h/4`
  and `h/8`, from (b) the channel-dependent square-root coefficients multiplying
  `A`, `A'`, etc. in Table V.
- [ ] Infer `g,h` from the same two calibration anchors used for `A,S0`, where
  `S0=3hβ`, and compare the resulting predicted class ratios with the paper's
  unity approximation. In this strict reduction,
  `A'=A-S0/6`, `A''=A-S0/24`, and `A0=A-S0/12`. As an immediate scale check,
  the current reproduction values `A=1.644`, `S0=3.291` give
  `A'/A=0.666`, `A''/A=0.917`, and `A0/A=0.833`; if isolated from all other
  effects, their squared width factors relative to the unity approximation are
  about `0.444`, `0.840`, and `0.694`.
- [ ] Quantify successively: the paper's common-`β`/unity-ratio baseline; exact
  Eq. (19) on the same SHO states; rms-matched state-dependent SHO surrogates;
  native unmixed GI waves; and native physical mixed states. This separates
  class-relation, size, node, relativistic-wave, and mixing effects.
- [ ] Audit every native result to ensure that an oscillator-solver basis `β`
  enters only through a converged representation of the wavefunction and is
  never reused as the physical `β` of the frozen decay form factor.
- [ ] State the remaining model limitation prominently: applying Eq. (19) to
  native parent and surviving-meson waves improves the two-wave spectator
  matrix element, but the emitted pseudoscalar is still treated as pointlike.
  A complete pair-creation amplitude is linear in all three meson wavefunctions
  and is a separate future operator, not something recovered by choosing `β`.
- [ ] Report amplitude ratios and partial-width ratios channel by channel,
  reserving the two calibration channels for calibration and testing the
  approximation on independent Table V channels.
- [ ] Build a provenance view of Table V grouped by reduced matrix-element
  class and flavor sector. Trace every predicted width to one of the two fit
  anchors, its explicit class ratio (`A'/A`, `D/S`, `A_c/A`, etc.), and any
  mixing, recoil, form-factor, or finite-final-width assumption.
- [ ] Compare three distinct hypotheses without conflating them: Eq. (19) with
  one universal `g,h`; the paper's two global reduced strengths with unity
  class ratios; and genuinely sector-dependent refits. Use independent rows
  for validation, and treat bottom as a new test rather than part of the
  original Table V evidence.

### Native Eq. (19) versus Table V amplitude atlas

This is the end-to-end validation and first physics application of the new
matrix-element API. The observable is not the special threshold moment of
Eqs. (20)--(21), but the complete on-shell partial-wave amplitude for each
supported Table V `M* -> M + P` row at the same external masses and momentum.
Keep three layers distinct:

1. the published Table V baseline (`LeadingS0`, common `beta=0.40 GeV`, fitted
   `A,S0`, unity reduced-class ratios);
2. the unreduced Eq. (19) operator on pure common-`beta` SHO states, with
   `g,h` inferred once from the same two calibration channels;
3. the same Eq. (19) operator and transported `g,h` on resolved native GI
   physical states and wavefunctions;
4. as a separately labelled comparison, native GI waves with `g,h` refitted
   to the same two anchors and no other rows.

Layer 1 to layer 2 isolates the paper's reduced-class and dropped-momentum
approximations. Layer 2 to layer 3 isolates realistic radial shapes, nodes,
spin-dependent wave distortion, unequal constituent masses, and coherent GI
mixing with the operator parameters held fixed. Layer 4 compares the models at
equal calibration cost and tests whether realistic waves improve the pattern
of independent channels rather than merely changing the overall scale. An
optional rms-matched SHO layer may later separate size changes from
non-Gaussian shape and node effects.

The raw native partial waves are dimensionless Appendix-C amplitudes, whereas
the Table V column is a width-normalized amplitude in `MeV^(1/2)`. Compare them
through
`a_L = sqrt(1000*q/(2*pi*(2J_i+1))) * H_L`, so that
`Gamma_L[MeV] = abs2(a_L)` and the sum over partial waves equals the public
native `decay_width`. Absolute amplitude signs require the documented common
external-state phase convention; widths and within-channel partial-wave
ratios are phase-safe. Use the same physical masses at all three layers so a
wavefunction comparison is not contaminated by changed phase space.

What is already wired: the complete 220-row canonical Table V transcription
and frozen reproduction backend; the generic Eq. (19) flavor-spin-angular
contraction; arbitrary orbital spatial integrals on analytic SHO, native HO,
and FD waves; coherent physical-state mixing; on-shell kinematics; partial-wave
projection for the Table V spin-zero and vector daughters; and native widths.
The analytic `rho -> pi pi` limit and the `A1 -> rho pi` S/D ratio already
recover Table IV in focused tests.

What is not wired: a canonical Table V state-assignment registry, automatic
construction of explicit charged/isospin flavor states, one report that runs
both backends on the same rows, a public/directly inspectable conversion to the
common width-amplitude convention, and line-shape integration for footnote-`f`
unstable daughters. The current realistic-factor audit is only the limited
Eqs. (20)--(21) diagnostic and is not this calculation.

- [ ] Add a page-verified Table V state registry mapping every parent and
  surviving daughter to `(n,L,2S+1,J)`, ordered flavor content, physical mass,
  and any Table III/Fig. 4 mixing prescription. Do not infer identities from
  display labels.
- [ ] Add a width-normalized partial-wave amplitude view shared by
  `PseudoscalarEmission` and `TableVReference`, with units and phase behavior
  explicit. Verify that its squared magnitudes sum to `decay_width` for both
  normalizations.
- [ ] Expose the exact linear coupling decomposition of the native operator,
  `H_L(q;g,h) = g*C_g,L(q) + h*C_h,L(q)`. The present API can already obtain
  it by evaluating `PseudoscalarEmission(1,0,masses)` and `(0,1,masses)`; add
  one discoverable coupling-basis helper rather than requiring fitted `g,h` or
  symbolic algebra. Verify recombination against arbitrary numerical `g,h`.
- [ ] For the old `A`-type classes, compare the phase-aligned ratio
  `C_h/C_g` directly with the SHO values `+1/4`, `-1/4`, `+1/8`, and `0`.
  This ratio cancels every common barrier, normalization, and SHO Gaussian.
  Also report `C_g^native/C_g^SHO` and `C_h^native/C_h^SHO` at the same `q`;
  unlike their ratio, these diagnose the overall spatial-overlap change.
- [ ] Separate physical-point and shape comparisons: evaluate every channel at
  its on-shell `q`, but also scan native and SHO coefficients on a common `q`
  grid. Divide only the exact threshold barrier `q^L`; do not present the SHO
  Gaussian as a model-independent factor of the realistic matrix element.

#### Flavor-scoped `A`-class universality survey

Table IV assigns the same `A=(g+h/4)beta` reduced amplitude to seven distinct
spectroscopic partial-wave transitions:

1. `1^3S1 -> 1^1S0 + P`, relative `L=1`;
2. `1^3P2 -> 1^1S0 + P`, relative `L=2`;
3. `1^3P2 -> 1^3S1 + P`, relative `L=2`;
4. `1^3P1 -> 1^3S1 + P`, relative `L=2`;
5. `1^1P1 -> 1^3S1 + P`, relative `L=2`;
6. `1^3D3 -> 1^1S0 + P`, relative `L=3`;
7. `1^3D3 -> 1^3S1 + P`, relative `L=3`.

Do not cross these seven rows mechanically with all ten spectrum sectors. The
Table-IV universality claim belongs to the light SU(6), common-oscillator
calculation. Use the five paper/repository sector colors `isovector`,
`isoscalar`, `strange`, `charmed`, and `bottom-flavored`, while retaining the
exact flavor pair inside each trace. Their mass-routing content is:

1. light-light (`n nbar`, `n sbar`, and `s sbar`), retaining the paper's SU(3)
   equal-mass approximation as one reference and the actual `n-s` mass routing
   as a separate test;
2. charm-light (`c nbar`, `c sbar`), using the explicit heavy-light recoil
   coefficient represented in Table IV by
   `A_c = [g + (1/2) m_l/(m_l+m_c) h] beta`;
3. bottom-light (`b nbar`, `b sbar`), as the controlled extension
   `A_b = [g + (1/2) m_l/(m_l+m_b) h] beta`, clearly labelled as beyond the
   Table-V numerical study.

Exclude charmonium, bottomonium, and `B_c`: they are not instances of light
pseudoscalar emission from the light constituent in the heavy-light formula.
Light isoscalar and isovector labels may select different physical flavor
compositions, but they do not constitute independent constituent-mass sectors
for the radial-overlap comparison.

Evaluate threshold-reduced off-shell kernels and common-`q` scans for each
applicable cell; attach physical on-shell channels only where a real Table V
decay exists. For every unequal-mass state, retain quark-emission and
antiquark-emission momentum routing separately. Compare against the exact
common-beta SHO evaluation with the same masses and routing; demand `1/4` only
in the equal-mass limit.

The unequal-mass factor is not an extra phenomenological correction. Starting
from Eq. (19), resolving the emitting constituent coordinate and momentum into
center-of-mass and relative variables introduces the constituent mass
fractions. For an `A`-type SHO sandwich with emitter `e` and spectator `s`, the
coefficient becomes
`g + (1/2) m_e/(m_e+m_s) h`, reducing to `g+h/4` at equal masses. The native
operator must obtain this through momentum routing and the wave integral, not
through an `A_c` or `A_s` hard-coded branch. The paper's selective choice to
show/use this correction for charm while keeping strange rows in the SU(3)
equal-mass approximation is the additional approximation to test. For strange
channels there is no single universal `A_s`: the fraction depends on whether
the strange or nonstrange constituent emits, so topology remains part of the
trace.

- [ ] Define a stable trace identity for every result:
  `(sector, flavor_pair, transition_id, parent_basis, daughter_basis,
  relative_L, emitter_topology, q, wave_backend, operator_piece)`. Store the
  resolved physical-state labels/masses/mixing coefficients, solver and wave
  provenance, `C_g`, `C_h`, SHO comparators, and all derived ratios alongside
  it. Never collapse transitions or sectors before writing this long-form
  record.
- [ ] Compute the threshold coefficients
  `lim(q->0) C_g(q)/q^L` and `lim(q->0) C_h(q)/q^L`, plus a certified common-`q`
  grid. Validate the limits by decreasing nonzero `q`, rather than dividing an
  exactly vanishing amplitude at `q=0`.
- [ ] Present views generated from the same trace: seven transition panels
  comparing the five sectors, and a primary sector-colored histogram of
  `U_hg = (C_h/C_g)_native / (C_h/C_g)_SHO`. Because the denominator uses the
  same masses and emitter routing, all sectors have the common SU(6)/SHO
  expectation `U_hg=1`; draw one vertical reference line there. Complement it
  with native/SHO magnitude ratios for `C_g` and `C_h`, also referenced to one.
  Preserve node/sign information and links from every plotted entry to its
  trace row. If either denominator is numerically small, use the phase-aligned
  normalized coupling-vector comparison and flag the entry instead of showing
  an unstable ratio.
- [ ] Infer one `g,h` pair from the same `rho -> pi pi` and
  `B -> [omega pi]_S` anchors, using the common-`beta` Eq. (19) bridge. Record
  the inferred values. Obtain them from the two complete finite-momentum SHO
  matrix elements as a linear two-parameter solve, not by silently identifying
  the paper's dropped-polynomial `S0` with the full second amplitude.
- [ ] Report both native calibration policies: transported SHO `g,h`, which
  isolates the raw wavefunction change, and a native two-anchor refit, which
  compares validation predictions at equal parameter count. Never mix the two
  in one ratio or refit to validation rows.
- [ ] Start with a phase-audited validation panel rather than all rows:
  `rho -> pi pi` (calibration and P wave), `A1 -> rho pi` (simultaneous S/D),
  one scalar `P P` S wave, one radial-node transition, one strange channel,
  one charmed unequal-mass channel, and one coherently mixed state.
- [ ] For every panel row, evaluate all labelled layers at identical external
  masses and report signed width-amplitudes, partial widths, layer ratios, and
  a term-level decomposition. Cross-check native HO against FD for the same
  solved states before interpreting a difference physically.
- [ ] Expand to every Table V amplitude row supported by resolved states.
  Classify rather than hide closed channels, missing experimental masses,
  mixing-only rows, quasi-two-body daughters, and footnote-`f` channels whose
  quoted width requires a daughter lineshape convolution.
- [ ] Produce the paper-facing figures only after the validation panel passes:
  reference versus native signed amplitudes by reduced-amplitude class;
  `Gamma_native/Gamma_TableV` on a logarithmic scale; and layer-2/layer-1
  versus layer-3/layer-2 panels separating reduction effects from realistic
  wavefunction effects. Mark the two calibration rows visibly and do not use
  near-zero amplitude ratios without also plotting absolute amplitudes.

### Paper Eqs. (20)--(21): realistic-wave factors

This exercise is already implemented by
`GIPaper/scripts/audit_realistic_factors.jl` and reported in
`GIPaper/docs/residual_reports/realistic_factors.md`. The current calculation
finds `R_A = 1.50, 1.81, 1.97` for decay `L=2,3,4` (paper `1.5`,
`1.7--1.8`, `2.1`) and `R_S=1.16` (paper `1.2--1.3`).

The operator interpretation should be made explicit. For a type-`A` amplitude,
the explicit emitted momentum in `σ·q` supplies one power of `q` and one rank-1
angular tensor. The threshold-leading plane-wave multipole is therefore
`j_(L-1)(qr) ~ (qr)^(L-1)`, giving the expected `q^L` barrier times the
internal moment `<f|r^(L-1)|i>`. In the common-SHO reduction, the recoil piece
that joins the same `A` class reduces to this same leading spatial structure.
For the type-`S` amplitude (`L=0`, internal `P -> S`), the plane wave starts at
one and `g σ·q` vanishes at threshold, while `h σ·p'` survives. Hence the
leading comparison contains the first-order vector momentum operator `p`, not
the positive quantity `p²`; in momentum space its radial integral has measure
`p² dp` and one operator power, i.e. `∫p³ Φ_S Φ_P dp`. The `g` term returns at
order `q²` in the structure-dependent `S(q)` polynomial.

Table V makes this harder to read than the prose suggests: its single
`"Realistic" factor` column interleaves both ratios without identifying their
kind. The factor must be inferred from the reduced-amplitude class and the
spin of the surviving daughter. The `B -> [omega pi]_S` row is deliberately
blank because it is the fitted type-`S` reference channel, with a `^3S1`
surviving daughter; it defines the denominator of Eq. (21), so its relative
factor is one. Its neighboring `B -> [omega pi]_D` row carries `(1.5)`, but
that is an `A`-class `D`-wave amplitude and therefore an Eq. (20) factor. The
actual Eq. (21) entries occur for type-`S` decays with a surviving `^1S0`
daughter, notably the `1^3P0` block (`delta2 -> eta pi`, `epsilon -> pi pi`,
etc.), where Table V prints `(1.3)`, and the excited/scalar strange rows where
it prints `(1.2)`.

- [ ] Add this threshold/multipole derivation to the generated realistic-factor
  report so Eqs. (20)--(21) are explained rather than merely transcribed.
- [ ] Make the digitized/report presentation label each nonempty realistic
  factor as `R_A` or `R_S`; do not require readers to reconstruct the kind from
  the row's amplitude class and daughter spin.
- [ ] Remove stale hard-coded prose and manifest metrics for the earlier
  `R_A=1.65,2.03,2.23` run; the report table currently contains the newer
  `1.50,1.81,1.97` values while a paragraph and `decays.toml` still quote the
  old results.
- [ ] Keep Eqs. (20)--(21) labelled as rough leading-operator diagnostics, not
  multiplicative corrections to the full finite-`q` native Eq. (19) amplitude,
  which retains the complete spherical-Bessel and derivative structure.

### Photon emission: Eq. (22) versus the mock-meson implementation

The paper presents Eq. (22) as the leading nonrelativistic one-body quark
current, but explicitly says that the Table VI numbers use a hybrid
mock-meson prescription so that the relativistic `m/E` ambiguity can be
included. The current `QuarkModelTransitions` implementation follows that
Appendix-D/Table-VI prescription rather than evaluating Eq. (22) literally.
Its M1 kernel is `mock_meson_overlap`, with pointwise momentum weight
`(1/m_i)(m_i/E_i)^0.7` and the mock-mass normalization. Its E1/M2 kernel is
`mock_meson_radial_moment`, with the global factor
`|m_i/sqrt(<E_i>_x<E_i>_y)|^0.5` multiplying `<x|r^n|y>`. The exponents 0.7
and 0.5 are phenomenological fits quoted by the paper, not consequences of
Eq. (22). They are distinct from the optional footnote-g recoil form factor
`exp(-q^2/(16 beta^2))` and from Appendix-A Hamiltonian smearing.

The primitive exported API remains available for kernel-level work:
`m1_transition_moment`, `e1_transition_amplitude`,
`spin_flip_photon_amplitude`, and the special
`m1_recoil_moment`; `m1_radiative_width` converts an M1 moment to a width.
The high-level path is now `PhotonEmission <: TransitionOperator`, evaluated as
`matrix_element(final_state, operator, initial_state)`. State spectroscopy
selects one of five typed kernels (`DirectM1`, `HinderedM1`, `AllowedE1`,
`SpinFlipE1`, or `SpinFlipM2`), while explicit flavor components determine the
standard current. `RadiativeAmplitude` retains the selected class,
`recoil_order`, and every coherently summed component term, and converts to a
width through the same `decay_width` name used by strong transitions. The
transition package contains no Table VI lookup: the complete 79-row audit owns
its footnote-to-order policy and a resolved-current adapter for legacy coarse
`(:q,:q)` states.

- [x] Design a typed photon-emission operator that dispatches through
  `matrix_element`, while retaining the existing primitive overlap functions
  as inspectable kernels. The final object is the surviving `PhysicalState`
  directly; the emitted photon is part of the operator/current, so no fake
  meson or one-element channel is introduced.
- [x] Make the choice between literal nonrelativistic Eq. (22), the published
  mock-meson hybrid, and any future current model explicit in the operator
  configuration and result provenance; never call all three simply Eq. (22).
- [x] Expose the `m/E` prescription coherently at the high level (including
  its fitted exponents), and report separately the mock-meson factor, recoil
  correction, and optional footnote-g form factor.

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
