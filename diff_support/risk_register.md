# Differentiation Risk Register

This is a design record, not an active work plan. It preserves the constraints
found by the differentiation audit so that a future implementation does not
distort GIModel's physics architecture. Differentiation remains an optional
beyond-paper capability.

| ID | Risk | Likelihood / impact | Early detection | Mitigation |
|---|---|---|---|---|
| R1 | An abstract central-potential field causes dynamic dispatch and `Any` return paths | mitigated / medium | type-stability probe | Keep aggregate component types concrete. |
| R2 | Forced `Float64` storage drops dual or tracked values | mitigated on central assembly / high for operator overloading | non-`Float64` assembly tests | Keep active leaves and central-kernel storage generic; presentation may remain `Float64`. |
| R3 | Constituent masses are rounded before physics evaluation | mitigated / high for small perturbations | exact-value/cache-key tests | Preserve exact promoted physics values and round only cache keys. |
| R4 | General eigendecomposition is unsupported or fragile in an AD backend | likely / high | smallest central derivative probe | Use an explicit symmetric-eigenvalue derivative rule. |
| R5 | Eigenvector derivatives become singular at small gaps | intrinsic / high | gap monitor and synthetic avoided crossing | Separate wave-dependent derivatives and report unresolved gaps. |
| R6 | Eigenvalue/state ordering changes at crossings | intrinsic / high | track gaps and overlaps across perturbations | Match states outside AD and reject ambiguous points. |
| R7 | HO beta selection changes under perturbation | likely in broad scans / high | record the selected beta and runner-up margin | Differentiate a fixed-beta local branch and expose branch switches. |
| R8 | Adaptive HO quadrature changes iteration count | possible / medium | record quadrature order under perturbations | Freeze a converged prepared quadrature outside the derivative trace. |
| R9 | Caches mix values from different parameter points | possible / high | cache-key/property tests | Keep active-value result caches outside the derivative kernel. |
| R10 | Mutation is removed merely for one backend | design risk / medium | allocation/performance benchmark | Select a compatible backend or isolate derivative rules; preserve useful local mutation. |
| R11 | Parameter-dependent guards introduce kinks | low in physical domain / medium | domain-boundary tests | Validate the physical parameter domain before entering the kernel. |
| R12 | Whole-package genericization creates compile-time/API cost | moderate / medium | method-instance and timing baseline | Genericize only active leaves and bounded numerical kernels. |
| R13 | A precise derivative is taken through an under-resolved representation | possible / high | independent FD-grid and HO-basis convergence | Converge both values and derivatives on each numerical route. |
| R14 | A partial derivative omits wavefunction response | likely for spin corrections / high | end-to-end finite-difference comparison | Label stopped-wave results explicitly; do not call them total derivatives. |
| R15 | Backend maintenance exceeds the physics value | moderate / high | dependency and extension review | Keep one small derivative contract and separate FD/HO implementations. |

## Preserved design decision

Do not force the entire model into one global numeric type. Parameterize related
continuous leaf groups and keep `GIParameters` concrete in its component types.
The numerical route remains a property of `FiniteDifferenceSolver` or
`OscillatorSolver`, never of the physics parameters.

A future derivative-facing kernel should return only numeric predictions from a
prepared, locally smooth FD or fixed-beta HO problem. Parsing, caches, adaptive
beta selection, report generation, state assignment, and experimental fitting
stay outside that kernel.

## Non-goals retained from the audit

- differentiating TOML parsing, reports, caches, or discrete solver choices;
- promising derivatives of labelled states through exact crossings;
- supporting every Julia differentiation backend;
- fitting experimental data as part of the differentiation interface; and
- attempting second derivatives before first-order behavior is stable.
