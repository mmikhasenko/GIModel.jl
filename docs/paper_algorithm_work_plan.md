# Paper-Order Spectrum Work Plan

Status: live code-stream queue, 2026-08-07.

This document turns the findings in
[`original_1985_algorithm_audit.md`](original_1985_algorithm_audit.md) into
small, dependency-ordered work units. It is the source-controlled follow-up
board for restoring the paper's HO calculation. Physics reproduction status
still belongs in `docs/paper_manifest/*.toml`; historical reproduction work
remains in [`work_plan.md`](work_plan.md); cross-cutting engineering invariants
remain in [`engineering_work_plan.md`](engineering_work_plan.md).

## Goal in one sentence

Once a radial problem has been solved, all later code should receive a solved
state and ask it for matrix elements or representations through dispatch. It
must not know whether that state came from HO coefficients or FD samples, and
it must not read a mesh from the state.

The qualification is important: an FD implementation necessarily owns grid
samples internally. The grid is an implementation detail of `MeshWave`, not an
argument or field used by physics consumers. An HO implementation owns `beta`,
basis labels, and coefficients internally and never needs to create a mesh for
production calculations.

## What exists and what is incomplete

The intended boundary is already described in `src/model_objects.jl`:

- `RadialWave` is the abstract representation of one radial eigenlevel;
- `MeshWave <: RadialWave` implements mesh quadrature;
- `MomentumWave` and `MeshMomentumWave` represent the transformed wave;
- `wave_norm`, `radial_expect`, `radial_overlap`, `momentum_wave`,
  `momentum_expect`, and `momentum_functional` are the beginning of a generic
  operation set.

The first implementation slice now meets the radial representation boundary:

- `OscillatorWave` retains `L`, `beta`, and normalized HO coefficients;
- `ChannelRadialSolution` stores native waves; the temporary tuple/virtual-mesh
  compatibility adapter has been removed;
- exact oscillator coordinate and momentum operations implement the same
  interface as `MeshWave`;
- spectrum, spin shifts, spectroscopic mixing, annihilation, and decay
  observables consume `RadialWave` operations; and
- `radial_wave(spectrum, ...)` returns the solver-native wave without exposing a
  mesh to the caller.

What remains is the larger paper-order problem: native HO *full fixed-sector*
Hamiltonians, complete mixing blocks, and final physical-state components. The
current spectrum still applies spin corrections in stages to a central wave,
and `radial_wave` on a mixed spectrum still resolves that central component
rather than the physical linear combination.

## Forensic trace of the unfinished abstraction

The repository history makes the intended continuation unusually explicit:

| Commit | Intended step | What remained at that commit |
| --- | --- | --- |
| `30dcbfb` - "A1: one interface for a radial wave" | Introduced `RadialWave`, renamed the concrete grid object to `MeshWave`, retained `RadialWaveOnUniformMesh` as an alias, and described a future `OscillatorWave(beta, coefficients)`. | Only the mesh implementation was added. The commit message counted 23 consuming functions but did not migrate their mesh-specific signatures and field access. |
| `522dcd1` - "A1b: momentum space is its own noun" | Added `MomentumWave`, `MeshMomentumWave`, and distinct quadratic/linear momentum operations. Removed `origin_amplitude` from the interface for the correct physics reason. | The older `RadialWave` docstring was not corrected and no oscillator momentum-wave implementation followed. |
| `1ef42ee`, `3a92ae8` - A2 slices | Migrated charge-radius and Appendix-D radial-moment calculations to `radial_expect`/`radial_overlap`; added scale-invariance checks. | These were two contained consumers, not a complete consumer migration. Spin, annihilation, and most transition code still require `MeshWave`. |
| `9c4a319` through `8655d68` - A3/A17 | Implemented exact HO `p^2`, position matrices, phase convention, and the mesh-free central A17 Hamiltonian. | `oscillator_channel_solution` keeps `(beta, coeffs)` only in a temporary `best` named tuple, reconstructs `best.basis * best.coeffs` on a reporting mesh, and returns only `(values, waves, r)`. The native solution is discarded at the solver boundary. |
| `de07702` - B3 solver types | Correctly moved the method choice out of `GIParameters` into `RadialSolver`, `FiniteDifferenceSolver`, and `OscillatorSolver`. | This means a new `SpectrumRepresentation`, `GIBasis`, `HORepresentation`, or `FDRepresentation` hierarchy would duplicate a distinction already owned by solver dispatch. |
| `fd281b1` - remove second wave cache | Removed a separate HO cache that had compensated for a normalization bug and established one solver per spectrum. | The one remaining `ChannelRadialSolution` is still a mesh-shaped container, so even an HO spectrum loses its native wave representation. |

That historical break has now been repaired: native HO information crosses the
solver boundary and downstream radial consumers dispatch on it. The next break
is one level higher, between a central radial component and the paper's full
fixed-`j,l,s` plus mixed physical state.

## Reuse map: entities to keep, extend, or retire

| Responsibility | Existing entity | Decision |
| --- | --- | --- |
| Numerical method selection | `RadialSolver`, `FiniteDifferenceSolver`, `OscillatorSolver` | Keep. Do not recreate a basis/representation hierarchy. |
| One radial eigenlevel | `RadialWave`, `MeshWave`, `OscillatorWave` | Keep exactly. Both native representations now exist; do not add `FDWave` or `HOWave` synonyms. |
| Backward-compatible mesh names | Removed | Do not restore aliases; this pre-release code has no compatibility obligation. |
| Momentum representation | `MomentumWave`, `MeshMomentumWave`, `OscillatorMomentumWave` | Keep without synonyms or another generic momentum abstraction. |
| Multi-level eigensolution | `ChannelRadialSolution` | Keep its native `Vector{<:RadialWave}`. Add convergence metadata here only when PA-12 defines it; do not add a parallel `FixedSectorSolution`. |
| Calculation cache/provenance | `SectorComputation` | Keep and evolve its cache key/value types. It already owns parameters, solver, and solutions. |
| Spectroscopic identity | `BasisState` | Keep as the canonical `n,L,S,J,label` identity. `FineStructureMultiplet` overlaps its `L,S,J` subset; retain compatibility constructors initially, then route physics methods through `BasisState` rather than inventing another identity type. |
| Spectrum stages and physical mass | `CentralState`, `CorrectedState`, `MixedState`, `Spectrum` | Extend these to retain or resolve their stage-correct radial components. Do not add `SolvedMesonState`. |
| Mixing algebra | `MixingBlock`, `MixingResult`, `diagonalize_mixing_block` | Keep. Annihilation and `same_j_mixing` should return or wrap `MixingResult` instead of returning named tuples that repeat `block`, `masses`, and `vectors`. |
| Per-state mixing projection | `StateMixing` | Keep the public role, but stop copying a block's masses and eigenvectors into every member. It should reference one `MixingResult` plus the selected eigenstate/column, with compatibility properties during migration. |
| Annihilation input metadata | `PseudoscalarAnnihilationBasisInput` | Keep the physics-specific mass/flavor metadata with `radial::RadialWave`. Do not replace it with a generic duplicate state wrapper. |
| Matrix-element name | exported `matrix_element` generic | Extend it with radial/operator methods if a general operator protocol is adopted. Do not create a second synonymous public function. |
| HO operator construction | `ho_p2_matrix`, `ho_operator_matrix`, `oscillator_momentum_factor_matrix`, ordinary matrices | Reuse as implementation primitives. Add small operator descriptors only where dispatch must preserve an unevaluated A17 sandwich; do not wrap every existing matrix in a second object graph. |

There are also two existing result duplications to remove during the migration:

1. `StateMixing` repeats `partner_masses_GeV` and one eigenvector column for
   each mixed spectrum member even though `MixingResult` already owns the full
   eigensystem.
2. Four annihilation paths construct a `MixingBlock` and then return named
   tuples repeating the block, masses, and vectors instead of returning a
   `MixingResult` plus only their genuinely additional diagnostics.

These are compatibility migrations, not reasons to add another state/result
hierarchy.

## Basis, radial solution, and physical state are different objects

"The basis has been computed" can mean two different stages:

1. **Prepared representation.** HO functions and operator matrices, or an FD
   grid and difference operators, are ready. The Hamiltonian still has to be
   assembled and diagonalized.
2. **Solved radial state.** Eigenvalue and normalized eigenvector are known. HO
   stores coefficients; FD stores samples. Downstream radial computations can
   now use the common interface.

A single solved radial function is enough only before state mixing. The final
meson can be a linear combination of several radial/spectroscopic/flavor
components: singlet-triplet, S-D, different radial levels, or different flavor
sectors. The public result therefore needs two levels of abstraction:

```mermaid
flowchart LR
    I["Model parameters and quantum numbers"] --> P["Prepared sector"]
    P -->|"HO backend"| HB["HO matrices; beta; basis"]
    P -->|"FD backend"| FB["FD operators; private grid"]
    HB --> D["Diagonalize full fixed-sector H"]
    FB --> D
    D --> R["RadialWave components"]
    R --> M["Mix compatible sectors and flavors"]
    M --> S["Existing Spectrum / MixedState: mass plus component amplitudes"]
    S --> O["Generic matrix_element / expectation / observable methods"]
    O -. "optional view only" .-> V["Sample on a plotting or export grid"]
```

`RadialWave` hides a representation of one component. The existing `MixedState`
inside `Spectrum` should own or resolve the physical mass and normalized
composition of such components. Observables should normally consume a spectrum
state; low-level radial kernels can consume `RadialWave`.

## Target contract

This is the reuse-first target. Names without a comment already exist.

```julia
abstract type RadialWave end                       # existing

struct OscillatorWave <: RadialWave                # implemented promised type
    # private implementation data: L, beta, normalized HO coefficients
end

struct MeshWave <: RadialWave                      # existing
    # private to mesh implementation: normalized samples and grid definition
end

struct ChannelRadialSolution{W<:RadialWave}        # current type
    eigenvalues_GeV::Vector{Float64}
    waves::Vector{W}
end

radial_wave(solution::ChannelRadialSolution, n)
radial_wave(spectrum::Spectrum, state::BasisState) # return stage-correct wave/components
state_components(spectrum::Spectrum, state::MixedState)

wave_norm(w::RadialWave)
radial_expect(w::RadialWave, f)
radial_overlap(left::RadialWave, right::RadialWave, f)
momentum_wave(w::RadialWave, L)
momentum_expect(w::MomentumWave, f)
momentum_functional(w::MomentumWave, f)
```

The existing small operation set covers local observables and transforms. It
does not yet cover the relativized Hamiltonian's nonlocal sandwiches such as
`B(p) f(r) B(p)`. PA-01 must decide the smallest extension: either one general
method on the already-exported `matrix_element` generic, or explicitly named
sandwich methods. If an abstract `RadialOperator` is introduced, it is the only
new operator abstraction and must describe unevaluated physics, not duplicate
the existing HO/FD matrices.

Useful operator families are:

- a local position kernel `f(r)`;
- a spectral momentum kernel `g(p^2)`;
- a sandwich `B(p^2) f(r) B(p^2)`;
- sums and scalar multiples of operators; and
- spin-angular factors separated from their radial operators.

For an existing `MixedState`, a matrix element is the component sum
`sum(conj(c_i) c_j <u_i|O_ij|u_j>)`. That rule belongs in one generic method,
not in each decay or spectrum report.

## What is needed after solving `u`

If `u` means a normalized eigenstate with its native representation retained,
then no further basis setup is needed. Later calculations require only:

1. its eigenvalue and quantum-number identity;
2. normalization and a deterministic phase convention;
3. generic diagonal and transition matrix elements in coordinate and momentum
   space;
4. composition coefficients if spectroscopic or flavor mixing has occurred;
5. convergence/provenance metadata so a result says how accurately it was
   solved; and
6. optional sampling for plots, exports, or explicit FD/HO comparisons.

The mesh is not one of these public requirements. The basis is also not a
public requirement. Both remain private data needed by the relevant dispatch
implementation.

## Work-unit rules

States used below:

- **done** - merged and covered by the normal verification gate;
- **ready** - dependencies are satisfied and the unit can be taken next;
- **pending** - waits for the listed dependency;
- **blocked** - needs a recorded physics or API decision.

Each unit should be one reviewable change with one primary responsibility. A
unit is complete only when its acceptance tests and documentation land in the
same change. During migration, FD, the current HO-hybrid path, and paper-order
HO remain explicitly named; no silent fallback is allowed.

## Dependency board

| ID | Status | Depends on | Deliverable | Acceptance |
| --- | --- | --- | --- | --- |
| PA-00 | done | - | Architecture audit of the paper and current paths | [`original_1985_algorithm_audit.md`](original_1985_algorithm_audit.md) identifies native, hybrid, partial, and missing stages. |
| PA-00R | done | PA-00 | Trace existing abstractions and history before naming new entities | The forensic trace and reuse map above identify the exact solver boundary where HO data is lost and prohibit duplicate wave, representation, solution, and physical-state hierarchies. |
| PA-01 | done | PA-00R | Freeze the reuse-first wave, solution, state, and operator contracts | Existing wave/solution/state types were retained; no parallel representation or solution hierarchy was added. |
| PA-02 | done | PA-01 | Implement native `OscillatorWave` | `L`, `beta`, normalized coefficients, analytic coordinate evaluation, and exact HO momentum transform are retained without a mesh. |
| PA-03 | done | PA-02 | Implement local position matrix elements for `OscillatorWave` | Diagonal/cross expectations use HO matrices or independent infinite-domain quadrature. |
| PA-04 | done | PA-02 | Implement momentum matrix elements and linear functionals for `OscillatorWave` | Momentum expectations, overlaps, and annihilation functionals use the exact oscillator representation. |
| PA-05 | done | PA-01 | Make `MeshWave` satisfy the same finalized contract | FD remains native mesh internally and shares the public normalization/expectation/overlap operations. |
| PA-06 | done | PA-03, PA-04, PA-05 | Complete the needed nonlocal radial primitive | Explicit two-sided momentum-sandwich methods dispatch for both backends; no general operator object graph was added because current consumers need only this primitive. |
| PA-07 | done | PA-06 | Migrate leaf observables from mesh fields | Spectrum spin shifts/mixing, annihilation, mock-meson overlaps, leptonic/two-photon observables, and charge radii accept `RadialWave`. |
| PA-08 | pending | PA-06 | Implement native HO contact matrix | Fixed-beta matrix elements reproduce the intended A15 kernel and independent refined quadrature, including the origin limit. |
| PA-09 | pending | PA-06 | Implement native HO vector and scalar spin-orbit matrices | A15/A16 momentum sandwiches use exact HO `p2` spectral factors and position matrices; no mesh operator is constructed. |
| PA-10 | pending | PA-06 | Implement native HO tensor matrices | Diagonal and cross-sector matrix elements pass refined-quadrature and symmetry tests. |
| PA-11 | pending | PA-08, PA-09, PA-10 | Add the full stage-1 builder to the generalized `ChannelRadialSolution` | Native-wave storage is done; the remaining acceptance is full fixed-`L,S,J` diagonalization so states sharing `L` may have different waves. |
| PA-12 | pending | PA-11 | Add beta refinement and basis convergence | Each sector records beta, basis size, tolerance, and convergence; requested radial range controls the variational objective rather than cache depth. |
| PA-13 | pending | PA-11, PA-12 | Route paper-order stage 1 through `compute_spectrum` | Contact, diagonal tensor, and symmetric spin-orbit terms are in the same diagonalization; paper mode has no fitted bridge factors or hidden mesh dependency. |
| PA-14 | pending | PA-10, PA-13 | Assemble complete stage-2 spectroscopic blocks | Tensor and antisymmetric spin-orbit blocks use stage-1 eigenvectors and all compatible requested radial states; pairwise results are recovered when other couplings vanish. |
| PA-15 | pending | PA-14 | Complete `MixedState`/`Spectrum` physical components | Existing state objects retain or resolve the full transformation; mass and observables use the same component waves rather than the central precursor. |
| PA-16 | pending | PA-04, PA-15 | Integrate stage-3 annihilation | Literal paper annihilation acts on stage-2 states in GIModel; calibrated comparison modes remain visibly separate in GIPaper. |
| PA-17 | pending | PA-07, PA-16 | Finish physical-state consumer migration | Consumers use the PA-15 mixed-state components rather than a central precursor; mesh-specific code remains only in the FD implementation and explicit diagnostics. |
| PA-18 | pending | PA-17 | Switch and certify the headline path | Full verification and convergence reports pass; HO paper mode, HO legacy-hybrid, and FD comparator are named; default changes only by an explicit reviewed decision. |

## Completed first follow-up: PA-01 through PA-07

The reuse-first slice completed the following checklist:

1. Retained `MeshWave`, `MeshMomentumWave`, and the already documented
   `OscillatorWave`; removed the old long-name aliases.
2. Added the smallest missing nonlocal primitive as explicitly named
   momentum-sandwich dispatch. No abstract operator holder was needed for the
   current consumers.
3. Specified how `ChannelRadialSolution`, `CentralState`/`CorrectedState`,
   `MixedState`, `StateMixing`, and `MixingResult` share rather than copy native
   waves and transformations. Mass and mixing coefficients do not belong to a
   bare radial wave.
4. Corrected the `RadialWave` docstring and implemented its promised HO form.
5. Added reusable contract tests covering normalization, phase,
   diagonal expectation, transition overlap, and coordinate/momentum duality.
6. Reduced direct mesh access to mesh implementations, explicit plot/export
   sampling, and the explicitly legacy `FDOriginP2Smearing` diagnostic.
7. Added migration regression tests before changing the `StateMixing` and
   `ChannelRadialSolution` storage layouts.

The next implementation target is PA-08 through PA-11: construct the native HO
contact, spin-orbit, and tensor matrices and diagonalize their sum with the
central Hamiltonian in one fixed-sector solve.

## Migration boundary

Direct access to `.u`, `.r`, and `.h` is legitimate only inside:

- the FD/mesh implementation of the wave and operator interfaces;
- constructors and validation of that representation; and
- an explicitly named diagnostic/export adapter.

It is not legitimate inside contact, fine-structure, mixing, annihilation,
decay, or report logic. Likewise, access to HO coefficients and `beta` is
legitimate inside the HO implementation, but not in consumers.

The practical review question for every later change is therefore simple:

> Could this method run on an HO state without constructing a mesh, and on an
> FD state without knowing the grid?

If the answer is no, the method is either a representation implementation or
the abstraction is still incomplete.
