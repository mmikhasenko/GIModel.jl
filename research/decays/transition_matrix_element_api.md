# Transition matrix elements: API and implementation architecture

Status: revised design, 2026-09-20. The staged delivery and acceptance gates
are in [matrix_element_api_implementation_plan.md](matrix_element_api_implementation_plan.md).

## Objective

Make a transition read like the physics,

```julia
amp = matrix_element(final, operator, initial; kinematics)
```

while separating physical-state composition, operator definition, exact
recoupling algebra, shared spatial integrals, helicity-to-partial-wave
projection, normalization, and provenance.

`matrix_element` means the complete transition. The current one-argument
method is renamed `reduced_matrix_element` at the package's `0.3.0` boundary.

## User-facing model

Resolve states once:

```julia
rho2  = physical_state(rho_spectrum, "2^3S_1")
pi    = physical_state(pi_spectrum, "1^1S_0")
omega = physical_state(omega_spectrum, "1^3S_1")

final = TwoMesonChannel(omega, pi)
```

`PhysicalState` owns a physical mass and the fully resolved basis components,
coefficients, waves, and provenance. It does not retain a `Spectrum` reference.
A `ReferenceState(label, mass)` represents legacy Table V daughters for which
no composite wave exists; native operators reject such a state.

Evaluate the operator once:

```julia
op = PseudoscalarEmission(g, h)
amp = matrix_element(final, op, rho2; kinematics=OnShell())

amp.helicity
partial_waves(amp)
amp[PartialWave(1, 1)]
decay_width(amp)
```

For an on-shell width query that may be closed,
`decay_width(final, op, rho2)` performs the threshold check first and returns
zero without constructing a nonexistent real-momentum amplitude.

The channel does not contain `(L,S)`. Those values are consequences of the
external quantum numbers and exchange symmetry. A single call returns the
complete helicity vector and every allowed partial-wave projection.

For loop work the same method is evaluated off shell:

```julia
vertex = matrix_element(final, op, rho2;
                        kinematics=CMKinematics(k))
```

`CMKinematics` accepts real or complex momentum. A closed on-shell channel has
zero width; asking for its vertex through `OnShell()` raises a targeted error
and requires an explicit analytic-continuation momentum.

## Result and normalization

The result records:

- operator, initial state, final channel, and kinematics;
- primitive helicity amplitudes;
- the complete derived partial-wave set;
- coherent component/topology terms;
- normalization and typed provenance.

Its fields and collections are concrete. It has no privileged scalar `total`,
because a physical decay may contain several correlated partial waves.

Width conversion belongs to the normalization:

```julia
decay_width(amp) = partial_width(amp.normalization, amp)
```

`RelativisticTwoBodyNormalization` is canonical. The legacy
`GITableVNormalization` remains an explicit reference convention, connected to
the canonical convention by Appendix C, Eq. (C2). There is no generic
`abs2(amplitude)` fallback.

Long-form `show(io, MIME"text/plain", amp)` prints the full decomposition.
Only if that becomes unwieldy should a narrowly named `term_table` accessor be
added. A generic `explain` verb is not part of the API.

## Calculation architecture

### 1. Resolve physical states

`physical_state(spec, label)` combines `physical_components`, the physical
mass, a stable identity, and provenance. All higher layers work with resolved
values and can be tested using synthetic states.

For

```math
|A\rangle=\sum_a c^A_a|a\rangle,\quad
|B\rangle=\sum_b c^B_b|b\rangle,\quad
|C\rangle=\sum_c c^C_c|c\rangle,
```

the internal composer computes

```math
\mathcal M_{A\to BC}
=\sum_{a,b,c}c^A_a(c^B_b)^*(c^C_c)^*
  \langle bc|T|a\rangle.
```

It retains every component triple. All terms use the same physical external
masses and kinematics.

### 2. Enumerate topology and derive pure coefficients

The operator supplies its physical interaction and allowed quark-line
topologies. Generic spin, flavor, color, and orbital algebra produces

```math
\mathcal M=\sum_{t,\mu}\xi_{t\mu} I_{t\mu}(k),
```

where `xi` depends only on quantum numbers and conventions, while `I` depends
on the operator, orbital label, waves, momentum, and numerical controls.

The internal seam is:

```julia
orbital_decomposition(op, topology, final_basis, initial_basis)
spatial_integral(op, label::OrbitalLabel, waves..., k)
```

This permits exact tests of coefficients and reuse of expensive integrals
across helicities, partial waves, multiplet members, and mixed components. A
spatial cache is never keyed by full spectroscopic identity when the wave and
orbital label are sufficient.

For GI elementary emission, Table IV's `A, A', A'', A0, S, D, P` classes are
the initial orbital-label basis. A `^3P0` operator may require a different one.

### 3. Form helicity amplitudes

Helicity amplitudes are the computational primitive, following Eq. (19)'s
recipe and Appendix C. The label type supports both daughter helicities rather
than encoding only the paper's `h_0,h_1` case.

### 4. Project all allowed partial waves

```julia
allowed_partial_waves(final, initial)
partial_wave_projection(final, initial)
```

derive the allowed `(L,S)` set and the Jacob--Wick linear map. Appendix C/Table
XI are fixtures for this generic layer. Several partial waves may share the
same helicity amplitudes and spatial-integral basis, so they are always
calculated together.

`TwoMesonChannel` preserves daughter order so an operator can assign physical
roles without a separate channel type. For GI elementary emission it is
`(surviving, emitted)`. Identical-particle exchange uses state identity, the
derived `(L,S)`, and phase `(-1)^(L+S)` to remove forbidden waves and apply the
normalization once.

### 5. Convert to observables

The normalization object owns amplitude units, phase space, symmetry factors,
and conversion to a partial width. The width sums the complete partial-wave set
once. Unsupported conversions fail through dispatch.

## Operator-specific numerics

### GI pseudoscalar emission

GI Eq. (19) is the first native operator. The pseudoscalar is elementary, so
the spatial part is a recoil- and derivative-dependent one-body operator
between parent and spectator-daughter waves. It can reuse the repository's
two-wave overlap infrastructure.

`g,h` are operator parameters. Legacy `A,S0,beta`, printed Table IV
polynomials, and `LeadingS0` are reproduction/calibration constructs. A
relationship is asserted only within the paper's single-beta SU(6) assumptions.

Native implementation activates the Appendix C helicity definitions and Table
XI projection that are folded into the current legacy tables.

### `^3P0` pair creation

All three mesons are composite. Its spatial amplitude is a recoil-shifted
momentum-space convolution. Before any public spatial interface is frozen, a
prototype compares direct multidimensional quadrature, an analytically reduced
integral, and a closed-form SHO reference. Routing, interpolation, bounds, and
convergence controls remain visible in provenance.

## Algebra and generation policy

Use generic Clebsch--Gordan and Jacob--Wick expressions, auditing
`PartialWaveFunctions.jl` for CG and Wigner-d conventions. Generate finite
6j/9j tables offline where appropriate. Flavor states are sparse sums of
ordered quark--antiquark pairs, and operator tensors plus rearrangement
topologies derive the flavor factors.

Every generated artifact has declared inputs, phase conventions, source
equations, deterministic output, and a stale-output check. Runtime symbolic
evaluation and Julia `@generated` functions are not initial dependencies.
The phase ledger incorporates `fix_outer_phase`,
`fix_annihilation_phase!`, and the oscillator-basis phase rule already present
in the repository. It verifies that the B31--B36 angular `(-1)^L` is introduced
exactly once by the Condon--Shortley spherical-harmonic convention.

## Phase conventions and mixed states

Preserve the radial-wave outer-lobe phase, mixing-eigenvector anchor, one flavor
phase convention, one helicity convention, and explicit daughter order. An
overall external-state phase may rotate an amplitude but cannot alter a width;
relative component and topology phases remain observable.

Discrete GI mixing is resolved into `PhysicalState`. Continuum-induced mixing
is not: an eventual coupled-channel layer retains the vertex vector and builds

```math
D^{-1}(E)=E\mathbf 1-M^{(0)}-\Sigma(E).
```

Its complex poles, left/right vectors, and residues require a separate
`ResonancePole` or `PoleState` representation and normalization semantics.

## Legacy compatibility

The Table V implementation is an explicit reference backend. It keeps the
canonical CSV, including quasi-two-body rows, and its current numerical path.
The one-argument `matrix_element(::StrongDecayAmplitude)` becomes
`reduced_matrix_element` in one breaking `0.3.0` commit without a deprecated
alias. `convention::Symbol` becomes singleton dispatch via
`TableIVPolynomial()` and `LeadingS0()`.

## What remains hardcoded

Keep explicit:

- the physical definition and parameters of each operator;
- calibration provenance;
- normalization and phase-space conventions;
- included phenomenological channels;
- approximations such as stable daughters or omitted topologies.

Derive selection rules, recoupling coefficients, flavor factors, partial-wave
projections, exchange symmetry, and SHO reductions.

## First vertical slice

1. Resolve `rho`, `pi`, `pi` into typed states and a daughter-only channel.
2. Derive the spin/flavor coefficient for `rho -> pi pi`.
3. Produce its helicity amplitude and derive the allowed partial wave.
4. Reproduce the analytic equal-beta legacy result.
5. Evaluate the same Eq. (19) operator on native HO and FD waves.
6. Display the complete term and convention decomposition.
7. Convert through both overlapping normalization conventions.
8. Vary real off-shell `k`, then test an explicitly supported complex `k`.

The next slice must have multiple partial waves from one helicity vector. The
third must be unequal-mass or `mixing_only`, so recoil and coherent composition
are exercised before generalization.
