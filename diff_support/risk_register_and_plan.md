# Risk register and staged plan

The transformation should preserve GIModel as a physics package. AD is a narrow
capability of the prediction kernel, not the package's primary architecture.

## Risk register

| ID | Risk | Likelihood / impact | Early detection | Planned response |
|---|---|---|---|---|
| R1 | `central` abstract field causes dynamic dispatch and `Any` return paths | mitigated / medium | type-stability probe | aggregate now carries concrete component types |
| R2 | forced `Float64` storage drops Dual/tracked values | mitigated on central assembly / high for operator overloading | non-`Float64` assembly tests | active leaves and central kernel storage are generic; presentation remains Float64 |
| R3 | constituent masses are rounded before physics evaluation | mitigated / high for small perturbations | exact-value/cache-key tests | exact promoted physics values; rounding only in cache keys |
| R4 | LAPACK eigendecomposition is unsupported or fragile in a backend | likely / high | smallest central spike | precompute `p²`; explicit symmetric eigenvalue rule |
| R5 | eigenvector derivatives blow up at small gaps | intrinsic / high for spin corrections | gap monitor and synthetic avoided crossing | central eigenvalues first; separate wave-dependent contract |
| R6 | eigenvalue/state ordering changes at crossings | intrinsic / high | track gaps and overlaps across perturbations | state matching outside AD; refuse/flag ambiguous points |
| R7 | HO beta grid selection changes under perturbation | likely in broad scans / high | perturb parameters and record selected beta | fixed-beta derivative context plus best/second-best selection margin |
| R8 | adaptive HO quadrature changes iteration count | possible / medium | log quadrature order under perturbations | keep outside initial trace or freeze prepared quadrature |
| R9 | caches mix values from different parameter points | possible / high | cache-key/property tests | no active-value result cache inside derivative kernel |
| R10 | local mutation is incorrectly removed for one backend | design risk / medium | benchmark allocations/performance | prefer Enzyme/Mooncake or isolate rules; do not optimize for Zygote prematurely |
| R11 | parameter-dependent guards introduce kinks | low in physical domain / medium | domain-boundary tests | validate positive-domain parameters before kernel |
| R12 | full package genericization causes compile-time/API cost | moderate / medium | method-instance and timing baseline | parametric leaves; concrete aggregate component types; Float64 reporting remains allowed |
| R13 | derivative is numerically precise for an under-resolved representation | possible / high | FD grid and HO basis/quadrature derivative convergence | convergence-check values and derivatives independently on each route |
| R14 | a partial derivative accidentally ignores wavefunction response | likely if spin stage is rushed / high | compare with end-to-end finite differences | do not call stopped-wave result a total derivative |
| R15 | backend/API maintenance becomes larger than physics work | moderate / high | dependency and extension review | one backend contract first, with separate FD/HO implementations; optional integrations only |

## Design decision: not one global numeric type

A single `GIParameters{T}` is sufficient for a simple ForwardDiff pass if every
active value shares `T`, but it couples unrelated groups and makes mixed numeric
contexts awkward. It is also unnecessary for Enzyme.

Prefer parametric leaf groups and an aggregate concrete in its component types:

```julia
struct ConfinementPotential{T<:Real}
    b::T
    c::T
end

struct RelativisticSmearing{T<:Real}
    sigma0::T
    s::T
end

struct GIParameters{P,C,S,R,F,A}
    potential::P
    central::C
    smearing::S
    factors::R
    fine_structure::F
    annihilation::A
end
```

The exact parameter list can be simplified during implementation, but the
principle is important:

- component fields are concrete;
- leaf constructors promote their own related scalars;
- users still write `GIParameters(...)` and need not spell type parameters;
- the numerical route is selected by `FiniteDifferenceSolver` or
  `OscillatorSolver`, not encoded in physics parameters;
- boolean/discrete configuration is not forced into an active scalar type.

## Stage 0: research package (this directory)

Status: complete for the initial audit. See [`audit_report.md`](audit_report.md).

Deliverables:

- mutation/control-flow/type-barrier inventory;
- active Hamiltonian equations;
- independent FD and fixed-beta HO Hellmann--Feynman probes;
- backend-specific research and source links;
- bounded transformation plan.

No production code or dependency changes are part of this stage.

## Stage 1: canonical numerical cleanup

Status: complete in the current implementation.

This stage must stand on its own without mentioning AD in its public API.

### 1A. Baselines

- Add inference checks for parameter lookup and both Hamiltonian/mass-vector
  routes.
- Record allocations and steady-state time for representative light, charm,
  heavy-light, and bottom sectors.
- Record FD grid convergence and HO basis/quadrature convergence separately.
- Keep the root and GIPaper regression suites as the physics oracle.

### 1B. Parameter ownership and constructors

- Build on the completed split into `FiniteDifferenceSolver` and
  `OscillatorSolver`, and removal of the phantom basis parameter from
  `GIParameters`; do not introduce a competing route marker.
- Make the central method a concrete field of `GIParameters`.
- Add copy-with-overrides constructors to each parameter group and the aggregate.
- Keep parameter-file parsing explicitly `Float64`; parsing is not a generic
  numerical kernel.
- Document which fields are continuous, discrete, or solver configuration.

### 1C. Exact masses and cache normalization

- Stop rounding inside `ConstituentMasses`.
- If solve coalescing still needs normalization, create it in
  `RadialChannelKey` or a dedicated cache-key constructor.
- Test separately that physics values remain exact and intended cache keys
  coalesce.

### 1D. Eigensolver contracts

- Replace runtime-`Symbol` branching across full/Krylov implementations with
  dispatch on concrete solver/algorithm context.
- Normalize each implementation to a concrete `(values, vectors)` return shape.
- Keep Krylov differentiation outside the first derivative contract.

Acceptance gate:

- no AD dependencies;
- model and report regressions unchanged within existing tolerances;
- `potential_diagonal` and both route-specific central mass-vector return types
  inferred;
- no material steady-state regression.

Estimated effort: 3--5 developer days.

## Stage 2: generic numerical representation

Status: complete for the bounded central-Hamiltonian target. Report and wave
presentation types remain deliberately `Float64`.

This is the real representation transformation, but not yet an AD promise.

### 2A. Parametric leaves

Parameterize:

- confinement and smearing groups;
- relativistic/fine-structure scalar groups where useful;
- `ConstituentMasses`;
- only the wave/result types needed by the target kernel.

Use promotion inside each leaf. Do not force all groups to one `T`.

### 2B. Generic central kernels

Remove active-path barriers from:

- central/smeared scalar functions;
- Hamiltonian diagonal and dense-matrix storage;
- relativistic kinetic and momentum-factor assembly;
- selected mass-vector output.

Keep report and CSV types `Float64` until a concrete derivative use requires
otherwise.

### 2C. Numeric-generic tests

- `Float32` scalar kernels and a small Hamiltonian;
- optionally `BigFloat` scalar/Hamiltonian checks where external LAPACK is not
  required;
- mixed input promotion tests;
- the existing `Float64` regression suite.

Acceptance gate:

- small FD and fixed-beta HO central Hamiltonians can be assembled with a
  non-`Float64` real type;
- `Float64` performance and results remain stable;
- no AD runtime dependency is required.

Estimated effort: 4--8 developer days.

## Stage 3: two prepared central prediction kernels

FD and HO get distinct prepared types and implementations constructed from
`FiniteDifferenceSolver` and `OscillatorSolver`, respectively. A common public
name may dispatch on them, but there should be no runtime route switch inside
the differentiated function and no route marker in `GIParameters`.

### 3A. Finite-difference route

```julia
struct PreparedFDProblem
    r
    h
    p2_values
    p2_vectors
    L
    nlevels
end
```

Preparation fixes the grid and validates its resolution. Repeated assembly uses
the prepared `p²` eigensystem.

FD acceptance gate:

- reproduce today's FD central solver on the same grid;
- no `p²` eigendecomposition inside active evaluation;
- full symmetric eigensolver initially; Krylov differentiation is separate;
- values and derivatives converge under FD grid refinement.

See [`routes/finite_difference.md`](routes/finite_difference.md).

Estimated effort: 2--4 developer days.

### 3B. Harmonic-oscillator route

```julia
struct PreparedHOBranch
    beta
    L
    nbasis
    nq
    p2_values
    p2_vectors
    quadrature_r
    quadrature_projector
    nlevels
end
```

Preparation fixes one beta branch and a converged quadrature order. A separate
`HOSelection` value records the selected beta and its margin to the second-best
candidate.

HO acceptance gate:

- reproduce today's fixed-beta HO matrices/eigenvalues;
- no `p²` or quadrature eigendecomposition inside active evaluation;
- fixed quadrature remains converged under tested perturbations;
- selected-branch derivatives agree with end-to-end differences whenever beta
  selection stays unchanged;
- beta switches are detected and reported as nonsmooth points;
- values and derivatives converge under basis/quadrature refinement.

See [`routes/harmonic_oscillator.md`](routes/harmonic_oscillator.md).

Estimated effort: 3--6 developer days.

Both routes' outputs should be numeric masses only, not `Spectrum`, dictionaries,
strings, warnings, or reconstructed waves.

## Stage 4: spectral primitive and backend experiments

### 4A. Explicit derivative oracle

- implement/test `E_dot[i] = dot(v_i, H_dot * v_i)` on both route matrices;
- implement the eigenvalue VJP `H_bar = V * Diagonal(E_bar) * V'`;
- test `dE/dc = 1`, the analytic `b` derivative, and finite differences;
- add gap diagnostics and define behavior below a documented gap tolerance.

### 4B. Backend spikes

In separate optional/test environments:

1. Enzyme native forward and reverse APIs on FD and fixed-beta HO assembly;
2. ForwardDiff on both generic assembly paths;
3. Mooncake reverse mode if readily compatible;
4. Zygote only through the isolated spectral rule, if cheap.

Measure compilation, steady-state time, allocations, and full-Jacobian cost
separately. A backend may pass one route and fail the other; record that rather
than treating support as automatically shared.

### 4C. Dependency/API decision

Possible outcomes:

- manual spectral derivatives plus Enzyme extension;
- lightweight `ChainRulesCore` rules plus one tested engine;
- a `DifferentiationInterface` facade if multiple backends pass;
- no runtime AD dependency, exposing only explicit JVP/VJP operations.

The last option is entirely acceptable for a project where differentiation is a
secondary feature.

Acceptance gate:

- all supported derivatives agree with route-specific finite differences and
  analytic oracles;
- behavior at small gaps is tested and documented;
- unsupported backends fail clearly rather than returning silent zeros;
- the chosen interface is smaller than the internal experiment matrix.

Estimated effort: 5--10 developer days, with backend compatibility the largest
uncertainty.

## Stage 5: optional spin/mixing expansion

Only start this if central derivatives are demonstrably useful.

Ordered difficulty:

1. nonperturbative contact Hamiltonian eigenvalues;
2. small mixing-block eigenvalues with fixed basis identity;
3. explicit parameter dependence of first-order spin operators;
4. central-wave response in first-order expectations;
5. phase-sensitive annihilation amplitudes and flavor/radial assignment.

Each item should extend the supported prediction contract separately. “Full
spectrum differentiable” should not be one milestone. Every item also needs an
FD implementation/test and an HO implementation/test: for example, the FD
nonperturbative contact solve acts on the radial grid, while the HO version
projects and diagonalizes in the finite oscillator space.

## Recommended pull-request boundaries

1. type-stability cleanup and copy constructors;
2. exact constituent masses plus cache-local normalization;
3. parametric leaves and generic central scalar kernels;
4. generic Hamiltonian storage plus numeric-type tests;
5. prepared FD workspace and mass-vector API;
6. prepared HO branch, quadrature preparation, and beta-selection margin;
7. shared spectral derivative oracle with separate route tests;
8. one backend experiment on both routes and a decision record.

Small boundaries make physics regressions attributable and allow stopping after
any stage without leaving a half-adopted AD framework.

## Explicit non-goals

- fitting experimental data;
- differentiating TOML parsing, reports, or caches;
- differentiating discrete solver/method choices;
- promising derivatives of labelled states through exact crossings;
- supporting every Julia AD package;
- second derivatives before first-order support is stable;
- rewriting performant local mutation merely for stylistic purity.
