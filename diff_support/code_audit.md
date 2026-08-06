# Code audit for differentiation

This audit asks two separate questions:

1. Is the present code canonical and type-stable enough to be a good numerical
   foundation?
2. Which operations belong inside a future differentiated kernel?

Mutation is not automatically a defect. Local writes into newly allocated
storage are useful numerical implementation details and are supported well by
compiler-level AD. Mutation of caches, output objects, or state identities is a
different concern and should remain outside the derivative boundary.

## Parameter-to-output path

For the active finite-difference central calculation, the path is:

```text
GIParameters + ConstituentMasses + finite-difference solver context
    -> radial_grid
    -> p2_operator
    -> eigen(p2)
    -> relativistic kinetic matrix
    -> smeared G and S
    -> optional A(p) * G * A(p) momentum sandwich
    -> symmetric H
    -> lowest_eigenpairs
    -> physical wave normalization
    -> central_spectrum cache and CentralState values
```

For the active harmonic-oscillator central calculation, the path is:

```text
GIParameters + ConstituentMasses + oscillator solver/beta context
    -> exact ho_p2_matrix
    -> eigen(p2)
    -> relativistic kinetic and A(p) matrices
    -> adaptive Gauss-Laguerre matrices for smeared G and S
    -> symmetric finite-basis H(beta)
    -> eigen(H(beta))
    -> discrete beta selection by the nlevels-th state
    -> optional mesh wave reconstruction
    -> central_spectrum cache and CentralState values
```

The narrow differentiable boundary should start after parsing and solver setup,
and stop at an ordered vector of eigenvalues. For HO it is a fixed-beta boundary;
quadrature and beta selection are prepared orchestration. Spectrum construction,
caching, and state presentation can consume either route's vector outside AD.

See the separate route audits:

- [`routes/finite_difference.md`](routes/finite_difference.md)
- [`routes/harmonic_oscillator.md`](routes/harmonic_oscillator.md)

## Continuous and discrete inputs

| Input | Current owner | First affected stage | Initial derivative support |
|---|---|---|---|
| `b`, `c` | `ConfinementPotential` | central Hamiltonian | yes |
| `sigma0`, `s` | `RelativisticSmearing` | smeared central potential | yes |
| constituent quark masses | `ConstituentMasses` / mass table | kinetic, smearing, all spin terms | yes |
| `epsilon_c` | `RelativisticFactors` | contact hyperfine | later |
| `epsilon_t`, `epsilon_so_vector`, `epsilon_so_scalar` | `RelativisticFactors` | fine structure | later |
| `k_spin_orbit`, `k_tensor` | `FineStructure` | fine structure | later |
| annihilation amplitudes and scales | `AnnihilationAmplitudes` | flavor/radial mixing | later/separate |
| numerical route and central-potential method | solver / singleton | algorithm selection | no; discrete |
| booleans in `RelativisticFactors`, `FineStructure`, `SpinTerms` | configuration | algorithm selection | no; discrete |
| `ngrid`, `rmax`, eigensolver, level count | `FiniteDifferenceSolver` | numerical method | no; fixed context |
| HO `beta`, basis size, quadrature order | HO solver | numerical method | fixed context; local selected-branch derivative |

The running-coupling coefficients and Gaussian scales are constants in
`src/constants.jl`, not fields of `GIParameters`. If they are eventually made
fit parameters, that is a model/API expansion rather than only an AD change.

## Pre-transformation type-stability baseline

This section records the audit snapshot. The issues were subsequently resolved;
the current inference results are in
[`implementation_status.md`](implementation_status.md).

In the audited snapshot, `GIParameters{FiniteDifferenceBasis}` is a concrete Julia type, but its
`central::CentralPotentialMethod` field is abstract. Therefore the aggregate's
name being concrete does not imply that field access is inferable.

The repeatable probe currently reports:

```text
fieldtypes(GIParameters{FiniteDifferenceBasis}) =
    (ConfinementPotential, CentralPotentialMethod, RelativisticSmearing,
     RelativisticFactors, FineStructure, AnnihilationAmplitudes)

potential_diagonal       => Any
relativistic_hamiltonian => Tuple{LinearAlgebra.Symmetric, Vector{Float64}}
FD channel_solution      => Tuple{Any, Any, Vector{Float64}}
HO fixed-beta Hamiltonian=> Tuple{LinearAlgebra.Symmetric, Matrix{Float64}}
HO channel_solution      => Tuple{Any, Any, Vector{Float64}}
```

This is a general cleanup issue, independent of AD. The central method should be
represented in the aggregate type, for example as a concrete field type
parameter. The scheduled solver split should also remove the unrelated basis
parameter from this physics aggregate. The probe lives at
`probes/type_stability_probe.jl`.

### Concrete numeric barriers

At the audit snapshot, the most relevant forced-`Float64` sites were:

- every continuous leaf parameter in `src/parameters.jl`;
- `ConstituentMasses`, including a `Float64` conversion and rounding in its
  physics-value constructor;
- `MeshWave` and `physically_normalized_waves`;
- `ChannelRadialSolution` and spectrum state/result fields;
- `MixingBlock` and `MixingResult`;
- explicit `Matrix{Float64}`, `Vector{Float64}`, `zeros(Float64, ...)`,
  `similar(..., Float64)`, and `collect(Float64, ...)` allocations;
- scalar accumulators initialized as `0.0` in kernels that may later receive a
  non-`Float64` element type;
- calls to `float`, which preserve `Float32` but are not a substitute for
  deliberate promotion and generic storage.

These barriers prevent an operator-overloading backend such as ForwardDiff from
propagating `Dual` values. They do not by themselves prevent Enzyme from
differentiating a `Float64` primal computation.

Do not genericize every report structure pre-emptively. The initial target only
needs generic parameter leaves, Hamiltonian storage, and mass output storage.

## Mutation inventory

### A. Local array construction: keep, but make storage element types correct

| Location | Mutation | Necessary? | Proposed treatment |
|---|---|---:|---|
| `radial_grid.jl:p2_operator` | fills tridiagonal diagonal | reasonable | keep; allocate from promoted element type or precompute |
| `hamiltonian.jl:nonrelativistic_hamiltonian` | fills diagonal | reasonable | keep; generic local storage |
| `smearing_appendix_a.jl:smear_3d_radial` | fills quadrature weights/output | reasonable | keep for Enzyme/Mooncake; genericize if ForwardDiff support is wanted |
| `radial_1d_coulomb_smear.jl` | fills output | reasonable | comparator path; not initial target |
| `appendix_a_derivative_potential.jl` | finite-difference output/boundaries | reasonable | comparator path; not initial target |
| `harmonic_oscillator_basis.jl` | fills basis/reconstructed waves | reasonable | exclude reconstruction from central-mass kernel; prepare fixed-beta HO operators separately |
| `pseudoscalar_annihilation.jl` | assembles small matrices | reasonable | later target with local generic storage |
| `physically_normalized_waves` | normalizes a new copy in place | reasonable | not needed for eigenvalue-only kernel |

These mutations do not alter model inputs or persistent global state. Rewriting
them as comprehensions solely for AD would be unnecessary unless Zygote is made
a primary backend.

### B. Eigensolver phase mutation: mathematically conventional, not mass-relevant

The code flips eigenvector columns in:

- `state_mixing.jl:diagonalize_mixing_block`;
- `pseudoscalar_annihilation.jl` phase helpers;
- `harmonic_oscillator_basis.jl:orthonormalize_physical_basis`;
- `spectrum.jl` via `fix_annihilation_phase!`.

An eigenvector sign is arbitrary. A branch that changes the sign to satisfy a
phase convention is discontinuous when its anchor component crosses zero.
Eigenvalues and phase-invariant expectations are unaffected. Therefore:

- no phase fix belongs in the central-mass derivative;
- wavefunction outputs should only promise derivatives for phase-invariant
  observables;
- phase-dependent annihilation amplitudes need a separately documented local
  convention and tests away from anchor zeros.

### C. Caches and orchestration mutation: exclude from AD

`central_spectrum` fills `Dict{RadialChannelKey,ChannelRadialSolution}` caches,
and `add_spin_corrections` uses a `get!` cache for contact-resummed levels. These
are useful at the orchestration layer but unsuitable as active derivative state.

The derivative kernel should instead receive an immutable/prepared workspace
containing constant arrays. A cache may store workspaces keyed by discrete
problem identity, but cache lookup and insertion should happen before AD.

### D. State assignment mutation: exclude from AD

`_assign_block_members!`, `_apply_same_j_spin_orbit_mixing!`, and
`_apply_tensor_mixing!` mutate vectors of `MixedState`. More importantly, they
sort masses and assign eigenvalues by rank. The concern is the discontinuous
identity decision, not the array write.

An eventual fit predictor must establish state correspondence outside the
derivative kernel. It can then differentiate the selected local branch while
monitoring for crossings.

### E. Unrelated public mutation

`MesonMasses.setindex!` and strong-decay channel-list construction do not lie on
the spectrum parameter-to-mass path. They need no change for this project.

## Control-flow and nondifferentiability inventory

| Construct | Dependence | Risk | Decision |
|---|---|---|---|
| basis, kinetic, eigensolver, and physics-term branches | discrete configuration | low if fixed | specialize or pass as constant context |
| `params.central isa ...` | discrete method stored abstractly | inference issue | make concrete in type; fixed for one trace |
| `max.(p2_eigenvalues, 0)` | numerical guard | kink at zero | verify spectrum is safely positive; prepared `p²` should make guard unnecessary in derivative path |
| `max(sigma, 1e-12)` | domain guard | kink near zero | enforce positive physical parameters before kernel |
| small-`r` analytic branches | fixed grid values | branch independent of active parameters | safe when grid is constant |
| naive 3D-smearing tail length `ceil(8/(sigma*h))` | active sigma | changes array length | exclude comparator path initially or fix integration domain |
| adaptive HO quadrature order | convergence test depends on values | piecewise algorithm | choose and validate fixed `nq` during HO preparation |
| HO beta grid `argmin` | active masses/potential | discontinuous selection | differentiate the unique selected branch; expose best/second-best margin and flag switches |
| full/Krylov eigenvalue sorting | eigenvalues | nonsmooth at crossings | full symmetric solve plus explicit state policy initially |
| mixed-state sorting/assignment | masses | discontinuous identity | outside kernel; monitor crossings |
| under-resolution warning | masses/waves | side effect only | omit from kernel; validate workspace beforehand |

## What can be prepared once

For fixed FD `ngrid`, `rmax`, and `L`:

- radial grid `r` and spacing `h`;
- finite-difference `p²` matrix;
- eigenvalues and eigenvectors of `p²`;
- constant running-coupling coefficient arrays;
- requested eigenvalue indices and quantum-number metadata;
- scratch matrices of the correct element type, if the selected backend supports
  explicit caches.

The current `p2_operator` accepts a mass but does not use it. The FD numerical probe
checks `p2(m) == p2(2m)` exactly. Removing the unused argument or moving `p²`
into a prepared workspace would clarify both the physics and the derivative
boundary.

For fixed HO `beta`, `L`, `nbasis`, and `nq`, prepare separately:

- the exact HO `p²` matrix and eigensystem;
- generalized Gauss--Laguerre radii and projector;
- the selected level count and beta-branch metadata.

The fixed-beta probe confirms that branch's spectral derivatives. Beta selection
must remain outside AD and carry a margin to the second-best grid point.

## General-cleanup findings and implementation state

The following cleanup was identified. Items 1--5 and 7--8 are implemented;
prepared route problems in item 6 remain Gate D:

1. Complete the scheduled separate FD/HO solver types and remove the basis marker
   from `GIParameters`.
2. Make `GIParameters.central` concrete through an aggregate field type.
3. Split the full and Krylov eigensolver paths into concrete return contracts.
4. Add consistent copy-with-overrides constructors for parameter groups.
5. Preserve exact constituent masses as physics values; round only in a cache
   key if cache coalescing remains desired.
6. Separate `PreparedFDProblem` and `PreparedHOBranch` from
   `SectorComputation`'s result cache.
7. Establish inference tests for both Hamiltonian and central mass-vector APIs.
8. Record separate FD and HO numerical/allocation baselines before
   genericization.

These changes are worthwhile even if the AD experiment is later abandoned.
