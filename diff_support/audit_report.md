# General numerical-code audit before differentiation

Audit snapshot: commit `7aef9a8` on 2026-08-06, plus the uncommitted
`gauss_laguerre_dvr` memoization work in `src/harmonic_oscillator_basis.jl`.
That file and its tests were treated as read-only concurrent work. Results that
involve the HO route therefore describe this worktree snapshot, not the last
commit.

Implementation update: Gates A--C were subsequently implemented and reverified
after the cache cleanup at commit `fd281b1`. See
[`implementation_status.md`](implementation_status.md).
The findings below are retained as the before-state and rationale.

## Executive conclusion

The transformation is realistic, but AD should not be the next production-code
change. The code first needs a small, useful cleanup that is justified even if
differentiation is never shipped.

The present parameter grouping is conceptually good: confinement, smearing,
relativistic factors, fine structure, and annihilation are separate objects.
Three representation details are not yet a good foundation for repeated
parameter evaluation:

1. the physics parameters carry the numerical basis as a phantom type;
2. the `central` field is abstract, so otherwise concrete parameter objects lose
   inference at the potential boundary;
3. continuous leaves and masses force `Float64`, and masses are rounded in the
   physics-value constructor.

The completed solver cleanup made the right first structural change: separate
`FiniteDifferenceSolver` and `OscillatorSolver`, then remove the basis marker
from `GIParameters`. Any future differentiation work should build on that
result rather than creating a competing basis-dispatch design.

Local mutation is not a general blocker. Most writes fill newly allocated
matrices or normalize new result arrays and can remain. Persistent caches,
adaptive choices, state assignment, and phase conventions must stay outside a
narrow derivative kernel.

## Findings at the audit snapshot, ordered by action priority

| Priority | Finding | Evidence | Required response |
|---|---|---|---|
| P0 | Solver identity belongs to numerical context, not physics parameters | `GIParameters{Basis}` and `with_basis`; the subsequently completed solver split | Complete the solver split and remove the basis parameter before diff-specific APIs |
| P1 | Central-potential dispatch is not inferred | `GIParameters{FiniteDifferenceBasis}` is concrete, but `central::CentralPotentialMethod`; `potential_diagonal` infers as `Any` | Parameterize the aggregate by the concrete central method (and, preferably, concrete component types) |
| P1 | The eigensolver wrapper has a union-shaped contract | `lowest_eigenpairs` branches on a runtime `Symbol`; inferred output contains real/complex and matrix/`Vector{Any}` alternatives | Separate full and Krylov implementations by dispatch and normalize each return type |
| P1 | Physics masses are rounded before evaluation | `ConstituentMasses` converts to `Float64` and rounds to 12 significant digits | Store exact promoted physics values; normalize only a dedicated cache key |
| P2 | Active scalar/storage types are fixed to `Float64` | parameter leaves, masses, Hamiltonian scratch arrays, waves and results | Genericize only continuous leaves and the selected central kernels after cleanup |
| P2 | HO quadrature memoization is global mutable state | worktree `_GAUSS_LAGUERRE_DVR::Dict` populated by `get!` | Treat it as preparation; define synchronization/ownership before concurrent use; never make it active AD state |
| P2 | High-level APIs combine numerical work with orchestration | caches, warnings, beta selection, wave reconstruction, state assignment | Introduce separate prepared FD and fixed-beta HO central-mass kernels |
| P3 | Several branches are intrinsically nonsmooth | beta `argmin`, eigenvalue ordering/crossings, rank-based state assignment, eigenvector phase fixes | Expose margins/gaps and differentiate only an identified local branch |

## Type-stability evidence before cleanup

The repeatable probe in `probes/type_stability_probe.jl` reports on the current
worktree:

```text
GIParameters{FiniteDifferenceBasis} concrete type: true
potential_diagonal          => Any
relativistic_hamiltonian    => Tuple{Symmetric, Vector{Float64}}
FD channel_solution         => Tuple{Any, Any, Vector{Float64}}
HO fixed-beta Hamiltonian   => Tuple{Symmetric, Matrix{Float64}}
HO channel_solution         => Tuple{Any, Any, Vector{Float64}}
```

The concrete name of `GIParameters` is therefore misleading as an inference
signal: its `central` field remains abstract. Calling
`central_potential_values(AppendixAMomentumSandwich(), ...)` directly is
inferred as `Vector{Float64}`, confirming the dispatch boundary as the cause of
`potential_diagonal`'s `Any` result.

This is not the sole issue. `lowest_eigenpairs` is not concretely inferred even
for `Symmetric{Float64,Matrix{Float64}}`. The full branch uses a dense symmetric
factorization with concrete output, while the Krylov branch introduces a broader
contract. Passing `eigensolver=:full` as a keyword does not eliminate that union
at the current API boundary. The solver split is a good opportunity to replace
runtime-symbol selection with method dispatch or concrete solver-algorithm
markers.

Inference gates should target stable public/internal boundaries, not demand that
every large orchestration function be fully inferred. The useful initial gates
are:

- central potential vector;
- FD central Hamiltonian and ordered mass vector;
- fixed-beta HO central Hamiltonian and ordered mass vector;
- full-eigensolver output contract.

## Parameter and ownership assessment

### Keep

- The physical grouping into small immutable structs.
- Discrete booleans and method choices as non-active configuration.
- TOML parsing as an explicitly `Float64` boundary.
- A separate constituent-mass object, because masses are sector input rather
  than global model parameters.

### Change during general cleanup

- Target an aggregate such as `GIParameters{P,C,S,R,F,A}`, with no numerical
  basis parameter and with concrete component fields.
- Add consistent copy-with-overrides constructors for leaf groups and the
  aggregate. These support ordinary parameter studies as well as later fitting.
- Let each related continuous leaf promote its own values. A single global
  numeric type `T` is unnecessary, and Enzyme does not require one.
- Make `ConstituentMasses` exact and parametric; move 12-significant-digit
  normalization to `RadialChannelKey` or another explicit cache-key constructor.
- Put route knobs in `FiniteDifferenceSolver` and `OscillatorSolver`. In
  particular, the HO operator should not carry irrelevant FD `ngrid`/`rmax`;
  a reporting reconstruction mesh is a separate concern.

### Defer

- Generic report, CSV, catalog, and presentation types.
- Spin mixing, annihilation calibration, and phase-dependent wave derivatives.
- A flat fitting vector or optimizer API. Parameter flattening belongs to the
  later fitting project and can be layered over stable model objects.

## Mutation and state audit

### Safe to retain inside numerical kernels

- filling fresh finite-difference diagonals and dense Hamiltonian storage;
- local quadrature work arrays;
- small local matrix assembly;
- normalizing a newly allocated wave array, where waves are actually requested.

These mutations neither change the parameter objects nor persistent model state.
They are compatible with compiler-level differentiation and can be made
ForwardDiff-compatible by choosing promoted element types.

### Keep outside the derivative boundary

- `SectorComputation` channel and wave dictionaries;
- contact-resummation result caches;
- the process-global HO DVR dictionary in the concurrent worktree;
- beta-grid selection and adaptive quadrature-order selection;
- warnings, report construction, and reference-state assignment;
- eigenvector sign fixes and mixed-state rank assignment.

The HO DVR cache stores parameter-independent quadrature data, which is a good
thing to reuse. Its current `Dict` plus `get!` form should nevertheless be
regarded as provisional orchestration: concurrent first writes need an explicit
policy, test isolation may need a reset/prewarm policy, and an AD call should
receive already prepared immutable arrays rather than populate the cache.

## Route-specific baseline

The script `probes/audit_baseline.jl` uses charmonium, `L=0`, the active
Appendix-A central path, FD `ngrid=450, rmax=24`, and HO
`beta=0.65, nbasis=24`. Built-in `@timed` measurements are intentionally
lightweight baselines rather than a benchmarking suite.

| Operation | Central values (GeV) | Warm median | Minimum allocation |
|---|---|---:|---:|
| FD assemble `H` (450×450) | — | 13.9 ms | 21.50 MB |
| FD assemble `H` + three eigenvalues | 3.0643437698, 3.6659025045, 4.0905614773 | 25.5 ms | 23.35 MB |
| HO fixed-beta assemble `H` (24×24) | — | 1.1 ms | 1.51 MB |
| HO fixed-beta `H` + three eigenvalues | 3.0646517538, 3.6663194062, 4.0910938089 | 1.1 ms | 1.52 MB |
| HO selected 22-beta route | 3.0645353473, 3.6662037469, 4.0909843271 | 23.5 ms | 30.47 MB |

For the current full 22-candidate HO beta selection, a separate clean-process
first call was 6.53 s and allocated 2.57 GB. The cold number combines compilation
and preparation and must not be read as steady-state physics cost. The large
separation from the warm row does show why preparation, cache ownership, and
cold/warm measurements must be explicit. Timings are machine-load-sensitive;
the allocations and returned values are the more useful regression anchors.

The FD matrix allocation is substantial for repeated evaluation and is a later
prepared-workspace opportunity. It is not evidence that mutation should be
removed.

## Feasible transformation sequence

### Gate A — integrate the active engineering work — complete

1. Let the HO memoization and B3 solver split land without overlapping edits.
2. Rebase/re-read the affected dispatch and test sites.
3. Rerun the inference and numerical baselines; run the complete root and
   GIPaper regression gates.

Exit: a stable baseline with distinct FD/HO solver ownership and no basis marker
in physics parameters.

### Gate B — canonical cleanup, no AD dependency — complete

1. Make the central method and parameter components concrete in the aggregate.
2. Add uniform copy-with-overrides constructors.
3. Preserve exact constituent masses and isolate cache-key normalization.
4. Split full/Krylov eigensolver implementations into concrete contracts.
5. Add inference regression tests at the four boundaries listed above.

Exit: identical physics reports within existing tolerances, inferred central
boundaries, exact mass perturbations, and no material warm-performance loss.
Estimated effort after B3: roughly 3--5 focused developer days.

### Gate C — narrow numeric genericity, still no AD promise — complete

1. Parameterize continuous leaf groups and `ConstituentMasses` independently.
2. Remove `Float64` barriers only from scalar physics and central matrix/mass
   output paths.
3. Add `Float32`, promotion, and small non-LAPACK generic assembly tests.
4. Keep report and wave-result genericity out unless a chosen derivative target
   requires it.

Exit: small FD and fixed-beta HO central matrices assemble with a non-`Float64`
real type while Float64 results and performance remain stable. Estimated effort:
4--8 developer days, with the HO route likely requiring more care around
quadrature storage.

### Gate D — two prepared differentiable experiments

Build separate prepared problems from the two solver types:

```julia
central_masses(theta, prepared::PreparedFDProblem)
central_masses(theta, prepared::PreparedHOBranch)
```

The HO problem represents one fixed beta and quadrature branch; beta selection
returns separate margin metadata. Both use an explicit symmetric-eigenvalue
JVP/VJP as the reference oracle. Try ForwardDiff and Enzyme against that oracle
before choosing any package dependency or public AD promise.

Exit criteria are route-specific derivative convergence, backend agreement with
finite differences and the spectral oracle, gap diagnostics, and measured
compile/warm costs. Failure of one route or backend does not invalidate the
other.

## Overall realism

- General cleanup: **low-to-moderate risk** and independently worthwhile.
- Numeric genericity: **moderate effort**, mainly mechanical but broad enough to
  demand regression discipline.
- Central eigenvalue derivatives on FD: **highly realistic**.
- Central eigenvalue derivatives on fixed-beta HO: **realistic**, with more
  preparation semantics.
- Differentiating beta selection, state identity, phase-sensitive waves, or the
  complete spin/annihilation pipeline: **not a reasonable first scope**.

This keeps differentiation subordinate to the package's main purpose. If the AD
backend experiment is later abandoned, Gates A--C still leave a cleaner,
better-inferred, and easier-to-study numerical model.

## Validation performed for this audit

- `probes/type_stability_probe.jl`: completed and reproduced the inference
  failures recorded above.
- `probes/central_hamiltonian_probe.jl`: completed; at relative step `1e-4`,
  direct and Hellmann--Feynman level derivatives agree to `9.6e-8` or better,
  with `max(abs(dH/dc - I)) = 1.4e-10`.
- `probes/ho_fixed_beta_probe.jl`: completed; the same comparison agrees to
  `5.5e-8` or better, with `max(abs(dH/dc - I)) = 3.8e-11`.
- `probes/audit_baseline.jl`: completed with the values and warm measurements
  above.

After implementation, the root and GIPaper package suites and the complete
`scripts/verify_project.sh` gate passed. See
[`implementation_status.md`](implementation_status.md) for the current state.
