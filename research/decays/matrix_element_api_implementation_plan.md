# Transition matrix-element API redesign: implementation plan

Status: implementation in progress; Phases 0--3 and package extraction complete,
2026-09-23; Phase 4 active

Companion documents:

- [transition_matrix_element_api.md](transition_matrix_element_api.md)
- [matrix_element_api_review.md](matrix_element_api_review.md)

## 1. Approved direction

Implement the redesign in stages with these decisions fixed up front:

1. `matrix_element(final, operator, initial; kinematics)` follows
   bra/operator/ket order and returns the complete transition.
2. Helicity amplitudes are the computational primitive. All allowed partial
   waves are derived together by an explicit Jacob--Wick projection.
3. `TwoMesonChannel` contains the daughters, not a user-selected `(L,S)`.
4. The calculation factorizes into pure angular/flavor/color coefficients
   times a small, shared basis of spatial integrals.
5. `physical_state(spectrum, label)` eagerly resolves mass, discrete GI mixing,
   waves, identity, and provenance. Matrix-element code never receives a
   `Spectrum`.
6. Relativistic two-body normalization is canonical for new amplitudes.
   `GITableVNormalization` is an explicit reference convention with its
   Appendix-C conversion.
7. The current Table IV/V calculation remains a frozen reference backend.
8. The first solver-native operator is GI pseudoscalar emission, Eq. (19).
   `^3P0` follows only after that architecture is validated.
9. Symbolic computation may generate and check finite algebra offline, but is
   not a runtime dependency.
10. Continuum dressing remains a separate post-GI layer consuming these
    off-shell vertices; it is not represented as a static `MixedState`.
11. `GIModel` owns the Schrödinger solver, spectra, waves, and generic overlap
    primitives. `QuarkModelTransitions` owns operators, matrix elements,
    channels, partial waves, and widths, with a one-way dependency on GIModel.

This revision adopts the review's structural findings rather than preserving
the earlier partial-wave-per-call design.

## 2. Scope

### In scope

- A coherent public API for complete strong-decay transition amplitudes.
- Eagerly resolved physical states and reference-only daughter states.
- Helicity amplitudes and their complete allowed partial-wave projections.
- Decomposition by physical-state component, topology, algebraic coefficient,
  and spatial integral.
- Coherent mixing of parent and both daughter states.
- Explicit on-shell, off-shell, and below-threshold behavior.
- A byte-identical legacy Table V adapter.
- Native HO and FD evaluation of GI Eq. (19).
- A later, separately gated `^3P0` spatial prototype.

### Out of scope

- Refitting the GI Hamiltonian or changing the 1985 reproduction parameters.
- Self-energy resummation, unitary scattering, or pole extraction.
- Three-body dynamics and unstable-daughter spectral functions. Stable
  reference objects for quasi-two-body Table V labels are nevertheless in
  scope because the dataset contains `epsilon`, `delta_2`, and `kappa`.
- Claiming a universal strong-decay operator.
- Runtime computer algebra.

## 3. Current-state facts

`QuarkModelTransitions/src/strong_decays.jl` provides `StrongDecayModel`, `DecayChannel`,
`StrongDecayAmplitude`, and `decay_amplitude`. Its one-argument
`matrix_element(a::StrongDecayAmplitude)` returns only `coefficient * reduced`.
That quantity is dimensionless; `spatial_overlap` and `total` carry the legacy
amplitude units. Only one test and one tutorial line call the misleading
one-argument method.

Reusable infrastructure already includes `RadialWave`, `OscillatorWave`,
`MeshWave`, `radial_overlap`, `momentum_wave`, `momentum_overlap`, state masses,
and `physical_components` with signed discrete-mixing coefficients.

Important constraints are:

- strong decay is one-to-two and may use three independently computed spectra;
- GI Eq. (19) is a two-wave elementary-emission operator;
- `^3P0` is a three-composite-meson shifted convolution, not three calls to
  `radial_overlap`;
- existing Table V coefficients and identical-particle factors remain an oracle
  until generic algebra derives them;
- the canonical CSV has 220 rows, including quasi-two-body daughters and rows
  whose physical amplitude arises only through mixing.

## 4. Non-negotiable invariants

### Physics

- Mixed-state amplitudes are coherent sums, never averages of widths.
- Daughter coefficients are complex-conjugated; the parent coefficient is not.
- Every component term uses the same physical external masses and momentum.
- Overall external-state phases cannot change a width; relative phases remain.
- Identical-particle normalization and exchange phases are applied exactly once.
- Allowed helicities and `(L,S)` follow from quantum numbers, parity, and
  identical-particle symmetry.
- All partial waves from one helicity calculation are retained and correlated.
- Forbidden channels are documented exact zeros; unsupported channels produce
  targeted errors.
- A closed on-shell channel has zero width, but its off-shell vertex is never
  silently replaced by zero.

### Numerical

- Analytic SHO formulas and direct numerical SHO evaluation agree where both
  exist.
- Native HO and FD routes agree within a convergence-derived tolerance.
- No silent interpolation, resampling, flavor fallback, or normalization change.
- Numerical controls and convergence results are stored in provenance.

### API

- `matrix_element` returns helicity and every allowed partial-wave amplitude.
- `decay_width` delegates to the normalization object and sums the full
  partial-wave set exactly once.
- Public values have concrete fields; decomposition and provenance are not
  stored as `Any`.
- Long-form `show(io, MIME"text/plain", amp)` exposes states, conventions,
  factors, units, and provenance.
- Existing Table V reports remain byte-identical under
  `bash scripts/verify_project.sh` until a separately approved data change.

## 5. Public API target

```julia
parent = physical_state(parent_spectrum, "1^3P_0")
d1 = physical_state(d_spectrum, "1^1S_0")
d2 = physical_state(k_spectrum, "1^1S_0")

final = TwoMesonChannel(d1, d2)
operator = PseudoscalarEmission(g, h, quark_masses)

amp = matrix_element(final, operator, parent; kinematics=OnShell())

amp.helicity
partial_waves(amp)
amp[PartialWave(0, 0)]
decay_width(amp)
```

For a width query that may be closed, a convenience overload resolves
`OnShell()` before attempting the vertex:

```julia
decay_width(final, operator, parent) # returns zero below threshold
```

Off shell, including an explicitly selected analytic continuation:

```julia
amp_k = matrix_element(final, operator, parent;
                       kinematics=CMKinematics(k_GeV)) # real or complex
```

An optional `partial_wave=PartialWave(L,S)` selector may be provided as a
convenience over the complete result; it is never the primitive calculation.

### Core domain types

```julia
abstract type TransitionOperator end
abstract type StrongDecayOperator <: TransitionOperator end

struct StateComponent{B,C,W}
    basis::B
    coefficient::C
    wave::W
end

struct PhysicalState{C,P}
    label::String
    mass_GeV::Float64
    components::C
    provenance::P
end

struct ReferenceState{P}
    label::String
    mass_GeV::Float64
    provenance::P
end

struct TwoMesonChannel{A,B}
    first::A
    second::B
end

struct PartialWave
    relative_L::Int
    channel_spin::Int
end

abstract type TransitionKinematics end
struct OnShell <: TransitionKinematics end
struct CMKinematics{T<:Number} <: TransitionKinematics
    momentum_GeV::T
end

abstract type AmplitudeNormalization end
struct RelativisticTwoBodyNormalization <: AmplitudeNormalization end
struct GITableVNormalization <: AmplitudeNormalization end
```

`ReferenceState` lets the legacy backend represent a label and mass without
inventing a radial wave. Native operators deliberately have no method for it.

The result contains concrete collections rather than an untyped scalar bag:

```julia
struct TransitionAmplitude{T,O,I,F,K,N,P}
    operator::O
    initial::I
    final::F
    kinematics::K
    normalization::N
    helicity::Vector{Pair{HelicityLabel,T}}
    partial_wave_amplitudes::Vector{Pair{PartialWave,T}}
    terms::Vector{TransitionTerm{T}}
    provenance::P
end
```

`TransitionTerm`, `StateComponent`, and the exact provenance type remain
internal initially. `HelicityLabel` supports general daughter helicities; it
does not bake the paper's `h_0,h_1` special case into the core type.

Width dispatch is readable and normalization-owned:

```julia
decay_width(amp::TransitionAmplitude) =
    partial_width(amp.normalization, amp)

partial_width(::RelativisticTwoBodyNormalization, amp) = ...
partial_width(::GITableVNormalization, amp) = ...
```

There is no universal `abs2(total)` fallback.

### Discoverability and export budget

Export only the user path: matrix-element and width verbs, physical-state and
channel constructors, partial-wave inspection/projection, kinematics,
normalizations, result type, and public operators. Keep topology enumeration,
coefficient generation, spatial integrals, component types, and term types
internal until a second operator proves those extension points useful. Use the
existing long-form `show`; add a narrowly named `term_table(amp)` only if the
display becomes too large. Do not add a generic `explain` export.

## 6. Internal calculation graph

```text
Spectrum + label
      |
      v
physical_state -> PhysicalState (mass + components + waves + provenance)
      |
      v
matrix_element(final, operator, initial; kinematics)
      |
      +--> allowed helicities / allowed (L,S)
      +--> component triples x quark-line topologies
      |       +--> pure coefficient decomposition xi_mu
      |       `--> shared spatial integrals I_mu(k; waves)
      +--> coherent helicity vector H
      +--> Jacob--Wick projection A^(LS) = P H
      `--> TransitionAmplitude -> partial_width(normalization, amplitude)
```

### State resolution and coherent composition

`physical_state(spec, label)` calls `physical_components` once and captures the
physical mass once. The internal one-ket/two-bras composer evaluates

```math
\mathcal M_{A\to BC}
=\sum_{a,b,c}c^A_a(c^B_b)^*(c^C_c)^*K(b,c;a)
```

and retains every component triple.

### Pure algebra times shared spatial integrals

For each component triple and topology, use

```math
\mathcal M=\sum_{t,\mu}\xi_{t\mu}(\text{quantum numbers})
I_{t\mu}(k;w_A,w_B,w_C).
```

The internal boundary is:

```julia
orbital_decomposition(op, topology, final_basis, initial_basis)
spatial_integral(op, label::OrbitalLabel, waves..., k)
```

The spatial cache key contains operator parameters, relevant topology, orbital
label, wave identities, momentum, and numerical controls—not the full external
spectroscopic identities. One integral may serve several `J` values, mixed
components, helicities, and partial-wave projections.

Table IV's `A, A', A'', A0, S, D, P` classes become elementary-emission
`OrbitalLabel`s; they are not assumed universal for `^3P0`.

### Helicity projection

The operator evaluator produces helicity amplitudes in a declared convention.
The separate projection layer provides:

```julia
allowed_partial_waves(final, initial)
partial_wave_projection(final, initial)
```

The second value is a testable linear map from the complete helicity vector to
the allowed `(L,S)` vector. Appendix C/Table XI are fixtures of the general
Jacob--Wick implementation, not operator code.

### Threshold policy

`OnShell()` requests real two-body kinematics:

- above threshold, compute the unique real center-of-mass momentum;
- at or below threshold, `decay_width` is exactly zero;
- the state/operator convenience `decay_width(final, op, initial)` performs
  that threshold check before constructing an amplitude;
- `matrix_element(...; kinematics=OnShell())` below threshold raises
  `ClosedChannelError` rather than pretending the vertex is zero;
- a closed-channel vertex requires explicit `CMKinematics(k::Complex)`.

Real off-shell `CMKinematics(k)` remains valid. Each operator documents its
analytic-continuation domain.

### Identical daughters

`TwoMesonChannel` canonicalizes daughter ordering. Once `(L,S)` is derived,
exchange symmetry uses `(-1)^(L+S)` together with internal quantum numbers to
remove forbidden waves and apply normalization once. Mixed daughters are
compared after component expansion, not merely by display label. The exact
component-identity key is a Phase 2 deliverable.

### Conventions and generated algebra

Replace `convention::Symbol` with dispatch types:

```julia
abstract type ReducedAmplitudeConvention end
struct TableIVPolynomial <: ReducedAmplitudeConvention end
struct LeadingS0 <: ReducedAmplitudeConvention end
```

Use `PartialWaveFunctions.jl` for Clebsch--Gordan and Wigner-d values after a
convention audit. Generate finite 6j/9j tables offline where clearer. Every
generator has declared inputs, a phase ledger, deterministic output,
source-equation references, and a stale-output check. Runtime symbolic
evaluation and `@generated` functions are excluded absent profiling evidence.
The ledger begins from the existing `fix_outer_phase`,
`fix_annihilation_phase!`, and oscillator-basis phase rules. It must explicitly
verify that the angular `(-1)^L` implied by the B31--B36 spherical-harmonic
convention enters exactly once through the Condon--Shortley algebra.

## 7. Implementation phases and gates

### Phase 0 — repair and freeze the legacy contract

Status: complete, 2026-09-20.

- Rename `matrix_element(::StrongDecayAmplitude)` to
  `reduced_matrix_element` in one breaking commit; add no deprecated alias.
- Bump package version `0.2.x -> 0.3.0` in that commit.
- Correct amplitude units and factorization documentation.
- Replace `:table_iv`/`:leading` with `TableIVPolynomial()`/`LeadingS0()`.
- Add factorization, dimensions, and convention-dispatch tests.

Gate: `bash scripts/verify_project.sh` passes; Table V residual reports are
byte-identical; search finds no old one-argument call or convention symbol; the
export count grows only by deliberate replacement names.

### Phase 1 — domain objects, projection skeleton, reference adapter

Status: complete, 2026-09-20. General Jacob--Wick projection remains the
explicit Phase 3 deliverable; Phase 1 supplies the interface and the unique
spin-zero-daughter projection.

- Implement the domain objects, concrete amplitude, and `physical_state`.
- Implement channel canonicalization, allowed-wave enumeration, and the
  projection interface.
- Encode the below-threshold policy now.
- Wrap the Table V path unchanged, including quasi-two-body reference daughters.
- Implement compact and long-form display.

Gate: typed and old reference calls agree on every supported Table V row;
quasi-two-body rows remain reproducible; native operators reject
`ReferenceState`; threshold behavior is explicit; inference/allocation checks
show concrete result storage; no unapproved exports appear.

### Phase 2 — coherent physical states and identical particles

- [x] Implement the internal three-state composer and complete term decomposition.
- [x] Finish identical mixed-daughter canonicalization and exchange symmetry.
- [x] Enforce one physical mass per external state.

Gate: synthetic complex-coefficient tests cover every combination of mixed
external states; overall-phase invariance and relative-phase interference pass;
a `mixing_only` fixture is zero in the pure basis and nonzero physically;
identical factors appear once; the identity rule enters the convention ledger.
The phase adds no public names.

### Phase 2.5 — package boundary

- [x] Extract transition-domain objects and the frozen Table IV/V backend into
  `QuarkModelTransitions`.
- [x] Keep the dependency one-way: `QuarkModelTransitions -> GIModel`.
- [x] Remove matrix-element, channel, and width exports from GIModel.
- [x] Update GIPaper to consume both packages and retain its complete test gate.

Gate: package-boundary tests pass; GIModel contains no transition dependency;
the transition and GIPaper suites pass; the Table V report remains byte-identical.

### Phase 3 — algebra and helicity-to-partial-wave vertical slice

Status: complete, 2026-09-23. The completion record is
`QuarkModelTransitions/docs/phase3_completion.md`.

- [x] Complete the state/helicity/flavor/topology/Jacob--Wick phase ledger.
  Radial, orbital, spin, flavor, physical-state, identical-daughter, Eq. (19)
  topology, spherical-contraction, and vector-pseudoscalar projection
  conventions are recorded. Constituent-mass `q'` and derivative-integral
  conventions are explicitly Phase-IV spatial inputs.
- [x] Audit `PartialWaveFunctions.jl` against Appendix C, every Table XI row
  through parent `J=5`, and projection orthonormality.
- [x] Implement normalized sparse Appendix-B flavor states and separate quark
  and antiquark emission contractions, including exact spectator/OZI zeros.
- [x] Implement the two-spin-1/2 coupled basis and spherical Pauli components,
  with sign-sensitive quark/antiquark and Hermiticity tests.
- [x] Implement the pure coefficient/orbital-integral decomposition, with
  provenance for every topology term and no numerical wave evaluation. The
  Eq. (19) spin-flavor/topology seam and its spherical-component orbital labels
  are implemented; coupling full spectroscopic states and shared radial labels
  is implemented as one coefficient matrix over a shared integral basis.
- [x] Derive the complete symbolic helicity vectors and all partial waves for
  the target set from that decomposition.

Targets: `rho -> pi pi`; one channel with two partial waves from one helicity
vector; one S-wave structure-sensitive decay; one unequal-mass charmed decay;
one isoscalar mixed decay; one forbidden channel; one identical-daughter
channel; and `f' -> pi pi` as a `mixing_only` row (or `E -> delta_2 pi` if its
source conventions make the sharper fixture).

Gate: exact signs/magnitudes match the coefficient oracle; Appendix C/Table XI
matrices agree; projection and recoupling identities pass; instrumentation
shows S/D amplitudes share the same integral basis; zeros are derived.
Only the previously budgeted projection inspection functions are exported.

### Phase 4 — solver-native GI pseudoscalar emission

Status: active, 2026-09-23. The first spatial checkpoint is recorded in
`QuarkModelTransitions/docs/phase4_spatial_checkpoint.md`.

- Implement Eq. (19) with `g,h`, analytic SHO, and native numerical waves.
- Compute the full off-shell helicity vector and all partial waves.
- Derive the Appendix C conversion to relativistic normalization.
- Keep `A,S0,beta`, `LeadingS0`, and Table IV polynomials explicitly
  reproduction/calibration-only. Relate them to `g,h` only within the paper's
  single-beta SU(6) assumptions.
- Update `docs/paper_manifest/appendix_cd.toml`: native evaluation un-folds C1,
  C3, and Table XI; the legacy backend remains directly partial-wave based.

Gate: equal-beta SHO formulas reproduce; native HO/FD agree under a stored
convergence certificate; node-sensitive examples refine stably; calibration
and validation sets are separate; the measured tolerance becomes the
regression threshold only after this study. The operator type is the only
planned new public name.

### Phase 5 — `^3P0` spatial prototype

- Fully specify operator, routing, recoil, topology phases, color, and norms.
- Determine its minimal orbital-integral basis.
- Compare direct multidimensional quadrature, an analytically reduced integral,
  and closed-form SHO results, including supported complex momenta.

Gate: all routes agree on declared benchmarks; cache keys share work without
conflating waves or numerical controls; unsupported domains fail explicitly;
no prototype-only spatial interface is exported.

### Phase 6 — public `^3P0` sector validation

- Stabilize internal dispatch from Phase 5 evidence.
- Cover a chosen heavy-light calibration and independent validation sector.
- Produce widths and branching fractions with sensitivity scans.
- Compare exact waves with equal-beta and rms-matched SHO surrogates.

Gate: no per-row formulas; results reproduce from parameters and state
identities; model limitations accompany conclusions; exports are justified by
user workflows.

### Phase 7 — documentation and release completion

- Rewrite the tutorial around resolved states, helicities, partial waves, and
  normalization-aware widths.
- Add an operator decision tree and document the post-GI continuum boundary.
- Run documentation-graph and public-link audits.

Gate: a user reaches a runnable example through docstrings; every export has a
purpose and test; verification and documentation builds pass.

## 8. Test matrix

| Concern | Required evidence |
|---|---|
| Legacy stability | verification script; byte-identical Table V reports |
| API contract | constructor validation, no ambiguities, concrete fields |
| Projection | Appendix C/Table XI fixtures; several waves from one helicity vector |
| Algebra | triangle/parity rules, exact identities, signs, flavor normalization |
| Mixing | conjugation, interference, `mixing_only`, phase invariance |
| Identical states | canonical order, exchange phase, forbidden-wave removal, one factor |
| Threshold | open on shell, zero closed width, explicit continued vertex |
| Spatial numerics | SHO analytic/numeric, HO/FD, nodes, refinement |
| Normalization | Appendix C conversion and width equivalence |
| Performance | cache reuse across waves/components and momentum-grid cost |

Performance never substitutes for numerical convergence.

## 9. Compatibility policy

The package is pre-1.0 and the old one-argument name has only two repository
call sites. The rename is therefore a clean `0.3.0` break, without an ambiguous
deprecation period. `StrongDecayModel`, `DecayChannel`,
`StrongDecayAmplitude`, and `decay_amplitude` remain. The canonical CSV adapter
and reports stay unchanged until typed replacements have full coverage.

## 10. Principal risks

| Risk | Mitigation | Gate |
|---|---|---|
| Partial wave treated as input | helicity primitive plus one projection | 1, 3 |
| Spatial cache keyed by full states | pure coefficients times wave/label-keyed integrals | 3, 5 |
| Ambiguous normalization | typed normalization and Appendix C conversion | 1, 4 |
| Closed channel silently zeroed | zero width only; explicit complex-k vertex | 1 |
| Phase mismatch | phase ledger and sign-sensitive fixtures | 2–4 |
| Identical mixed daughters double counted | canonical channel and exchange tests | 2 |
| Quasi-two-body rows lost | `ReferenceState` accepted by legacy only | 1 |
| `^3P0` reduced to radial overlaps | mandatory three-route prototype | 5 |
| Shifted-momentum interpolation bias | refinement and SHO comparisons | 5 |
| Combinatorial mixed-state cost | exact-zero pruning and shared integral cache | 2–5 |
| Type/export proliferation | concrete values and export budget | all |
| Runtime CAS dependency | deterministic offline generation | 3 |
| `A,S0` identified with `g,h` | limited-assumption calibration statement | 4 |
| Static and continuum mixing conflated | separate propagator/pole layer | docs |

## 11. Settled and evidence-gated decisions

Settled: argument order; eager `PhysicalState`; internal composer; relativistic
canonical normalization; GI Eq. (19) first; `PartialWaveFunctions.jl` for
CG/Wigner-d after tests; offline finite 6j/9j generation as needed; singleton
legacy conventions; immediate `0.3.0` rename; and reproduction-only status for
`A,S0,beta` outside their single-beta SU(6) interpretation.

Two decisions are deliberately evidence-gated:

- the exact identity key for identical mixed daughters is finalized at the
  Phase 2 gate;
- HO/FD tolerances come from the Phase 4 convergence certificate.

Neither blocks Phases 0 or 1.

## 12. Completion levels

The word "complete" has three deliberately separate meanings:

### A. API foundation complete — Phases 0--3

The domain model, legacy adapter, coherent mixing, coefficient/integral
factorization, helicity representation, and partial-wave projection exist and
are tested. The architecture is usable for the reference backend, but it does
not yet establish a solver-native strong-decay prediction.

### B. Transition project complete — Phases 0--4 plus applicable Phase 7 work

The project requested here is complete when the GI Eq. (19) operator works
end-to-end with resolved physical states, analytic SHO and native HO/FD waves,
real and supported complex momenta, all allowed partial waves, explicit
normalization, mixed-state interference, legacy cross-checks, convergence
evidence, and user documentation.

This is the principal completion milestone. It produces a validated off-shell
vertex suitable as input to a future continuum calculation. It does not claim
that `^3P0` or continuum dressing has been completed.

### C. Extended strong-decay program complete — Phases 0--7

The `^3P0` prototype and one calibrated/validated sector are also complete,
with branching fractions, numerical and parameter sensitivity studies, and a
documented comparison to the GI emission operator. This demonstrates that the
API truly supports a second operator rather than only abstracting one.

Continuum self-energies, pole searches, and scattering remain a separate
project after level C.

## 13. Definition of done for the transition project

A user can resolve independent physical states, select an operator, request an
on-shell or explicit real/complex off-shell amplitude, inspect the primitive
helicity vector and all allowed partial waves, inspect coherent terms, obtain a
normalization-aware width, reproduce legacy Table V including quasi-two-body
rows, switch between native HO and FD waves without changing the public call,
derive supported algebra without per-row production formulas, and reuse the
off-shell vertex in a separately scoped continuum calculation.

In addition, the following project-level evidence must exist:

- the complete verification script passes with byte-identical legacy reports;
- at least one channel with multiple partial waves proves shared-integral reuse;
- at least one `mixing_only` channel proves coherent composition;
- at least one node-sensitive channel has an HO/FD convergence certificate;
- analytic SHO, numerical SHO, HO, and FD routes agree in their overlapping
  domains;
- Appendix C normalization and projection conventions are traceable through the
  phase ledger;
- the public tutorial contains on-shell, real off-shell, complex off-shell,
  mixed-state, and failure-mode examples;
- unsupported states, operators, continuations, and normalization conversions
  fail explicitly rather than falling back.

## 14. Load and complexity estimate

Assumptions: one contributor already comfortable with Julia and the GI paper;
focused person-weeks; existing wave solvers and Table V data remain usable; code
review is included, but a publication and a broad phenomenological survey are
not. Calendar time will be longer if this is interleaved with other work.

| Phase | Main load | Complexity | Estimate |
|---|---|---:|---:|
| 0 | rename, units, convention types, freeze legacy reports | medium | 0.5--1 week |
| 1 | resolved states, channels, kinematics, result and reference adapter | medium-high | 1.5--3 weeks |
| 2 | coherent three-state mixing and identical daughters | high | 1.5--3 weeks |
| 3 | phase ledger, exact recoupling/flavor algebra, helicity projection | very high | 2--4 weeks |
| 4 | native Eq. (19), SHO/HO/FD paths, normalization and convergence | very high | 3--6 weeks |
| 5 | three-route `^3P0` spatial research prototype | research/high uncertainty | 3--7 weeks |
| 6 | `^3P0` sector calibration, validation, and sensitivity analysis | research/high uncertainty | 4--8 weeks |
| 7 | tutorials, API navigation, release and documentation audits | medium | 1--2 weeks |

Expected totals:

- API foundation (A): **6--11 person-weeks**.
- Transition project with native GI emission (B): **10--19 person-weeks**,
  including the relevant documentation work.
- Extended two-operator program (C): **17--34 person-weeks**.

For a single researcher/developer working part-time and allowing for review
cycles, level B is realistically a **3--6 month calendar project**. Level C is
closer to **6--12 months**. The ranges should be revised after Phase 3: by then
the phase conventions and integral-reuse architecture—the largest risks to the
Eq. (19) estimate—will have been tested.

Likely implementation volume is roughly 1,500--3,000 lines of library,
generated-table, test, and documentation changes through level B. This is only
a planning indicator: acceptance is evidence-based, not line-count based.

### Critical path and risk concentration

The critical path is `0 -> 1 -> 2 -> 3 -> 4`. Phase 3 is the conceptual gate:
if signs and Table XI projections are not independently reproduced, numerical
operator work must not proceed. Phase 4 is the numerical gate: the transition
project is not complete until node-sensitive HO/FD agreement is demonstrated.

The largest schedule uncertainties are:

1. reconciling paper, spherical-harmonic, flavor, and stored-wave phases;
2. defining identical-particle symmetry for mixed component states;
3. obtaining stable derivative/recoil integrals on FD waves;
4. establishing complex-momentum support without overstating the analytic
   continuation of numerical wave representations;
5. for the extended program only, reducing and validating the shifted `^3P0`
   convolution.

Phases 0--2 are predominantly engineering. Phases 3--4 are mixed
physics/engineering research. Phases 5--6 are research tasks and should not be
managed as fixed-scope feature work.

## 15. Readiness

The architecture is ready to implement. Phases 0 and 1 have no unresolved
design dependency. Later uncertainty is contained behind explicit gates: the
identical-mixed-state identity rule is a Phase 2 result, numerical tolerance is
a Phase 4 result, and no `^3P0` interface is frozen before its Phase 5 spike.
