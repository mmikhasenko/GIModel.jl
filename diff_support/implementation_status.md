# Differentiation-enabling cleanup: implementation status

Implementation snapshot: integrated and reverified on top of commit `fd281b1`,
2026-08-06.

Gates A--C from the audit are implemented. No AD package or runtime dependency
has been added. The next bounded step is Gate D: separate prepared FD and
fixed-beta HO kernels, an explicit symmetric-eigenvalue derivative oracle, and
backend experiments.

## Implemented

### Concrete model and solver dispatch

- `GIParameters{P,C,S,R,F,A}` carries each concrete component type, including
  the central-potential method.
- `FiniteDifferenceSolver{Kinetic,Eigensolver}` retains the public Symbol-based
  constructor and fields while encoding both choices in its concrete type.
- Relativistic/nonrelativistic FD channel solves dispatch separately.
- Full and Krylov eigenpair implementations dispatch separately and normalize
  their result to concrete `(Vector, Matrix)` contracts.
- `Meson`, `SectorComputation`, and `Spectrum` carry their concrete mass,
  parameter, and solver ownership types rather than reintroducing abstract
  fields at the orchestration layer.

### Parameter updates and exact physics inputs

- All parameter groups and `GIParameters` have copy-with-overrides constructors.
- Continuous leaf groups promote their own related scalars independently:
  `ConfinementPotential{T}`, `RelativisticSmearing{T}`,
  `RelativisticFactors{T}`, `FineStructure{T}`, and
  `AnnihilationAmplitudes{T}`.
- `ConstituentMasses{T}` promotes and preserves exact physics values.
- Twelve-significant-digit normalization now occurs only in
  `RadialChannelKey`; tests distinguish exact physics values from normalized
  cache identity.

### Narrow numeric genericity

- FD relativistic and nonrelativistic central matrix assembly no longer forces
  active storage back to `Float64`.
- HO Gauss--Laguerre operator values retain the scalar type returned by the
  operator function.
- Small FD and fixed-beta HO central Hamiltonians assemble as `BigFloat`, while
  ordinary loaded parameters still produce the original `Float64` results.
- Wave/result/report structures remain `Float64`; they are outside the initial
  eigenvalue-only derivative boundary.

## Inference gate

`probes/type_stability_probe.jl` now reports concrete returns for:

```text
potential_diagonal              Vector{Float64}
relativistic_hamiltonian        Tuple{Symmetric{Float64,Matrix{Float64}},Vector{Float64}}
full lowest_eigenpairs          Tuple{Vector{Float64},Matrix{Float64}}
HO fixed-beta Hamiltonian       Tuple{Symmetric{Float64,Matrix{Float64}},Matrix{Float64}}
FD channel implementation       Tuple{Vector{Float64},Matrix{Float64},Vector{Float64}}
HO channel implementation       Tuple{Vector{Float64},Matrix{Float64},Vector{Float64}}
```

The root test suite additionally applies `@inferred` to these boundaries and to
both full and Krylov eigenpair paths.

## Verification

- FD Hellmann--Feynman probe: passed, including `dH/dc = I` and analytic
  `dH/db` checks.
- Fixed-beta HO Hellmann--Feynman probe: passed, including `dH/dc = I`.
- GIModel `Pkg.test()`: passed.
- GIPaper `Pkg.test()`: passed.
- `scripts/verify_project.sh`: passed, including report regeneration, FD/HO
  audits, data checks, and manifest validation.
- Generated numerical report content had zero drift; three generation-time-only
  changes were restored to keep the worktree focused.

## Deliberately not implemented yet

- No `PreparedFDProblem` or `PreparedHOBranch` production API.
- No differentiation through LAPACK or KrylovKit.
- No eigenvalue JVP/VJP implementation or ChainRules rule.
- No ForwardDiff, Enzyme, Mooncake, or DifferentiationInterface dependency.
- No derivative promise for beta selection, adaptive quadrature selection,
  eigenvectors, spin mixing, annihilation, state assignment, or reports.

These are Gate D decisions. The representation cleanup no longer blocks those
experiments, and it remains independently useful if no backend is adopted.
