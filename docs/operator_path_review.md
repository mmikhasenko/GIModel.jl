# Orbital execution-path review

Date: 2026-09-27. Scope: `src/` Hamiltonian construction, contact diagnostics,
spin contributions, mixing, and orbital branches in wave operations.

## Defect and correction

The original baseline (`6d71a35`) suppressed contact for every L>0 despite a
finite-width kernel. `6519f63` asserted that behavior in a test, and `b9d3474`
carried it into native HO and fixed-sector FD assembly. The physical correction
was committed in `4c84a9e`; this follow-up removes the architecture that let
separate solving/reporting paths repeat the assumption.

The remaining contact-only diagnostic also excluded L>0 and hardcoded zero
orbital momentum in its explicit-grid solve. Both overloads now delegate to the
production fixed-sector solver for every L. Invalid multiplicities and meshes
throw; empty results no longer represent unsupported physics.

`ContactHyperfine` is the single definition for contact strength, kernel and
prescription. Both native matrix adapters serve solving and reporting. Explicit
legacy local and momentum-sandwich entry points remain as thin adapters. The
local approximation can be represented in either basis; production parameters
continue to select the Appendix-A sandwich.

## Retained branches and their reasons

| Location | Branch | Mathematical or API reason |
|---|---|---|
| `spin_fine_structure.jl` | diagonal L>0, multiplicity 3 | Diagonal L·S and tensor angular factors vanish for L=0 or S=0. Centralized for FD, HO and reporting. |
| `spin_fine_structure.jl`, `spectrum.jl` | antisymmetric spin-orbit L>0, J=L | The singlet/triplet matrix element has factor sqrt(L(L+1)); only J=L has both spins. |
| `spectrum.jl` | tensor triplet, J>0, L=J±1 | Off-diagonal tensor connects the same-J pair separated by ΔL=2. Includes S-D. |
| `spectrum.jl` | S-wave triplet J=1 | Angular-momentum addition; avoids enumerating invalid J=-1,0 triplets at L=0. |
| `flavor_mixing.jl`, `pseudoscalar_annihilation.jl` | prescribed pseudoscalar labels | Explicit scope of the phenomenological annihilation model; other states are rejected, not silently assigned an annihilation calculation. |
| `harmonic_oscillator_basis.jl` | L=0 at r=0 | Analytic origin limit of the radial basis; positive-r higher waves still use the general formula. |
| `momentum_waves.jl` | spherical Bessel L=0 | Recurrence seed for j_L; subsequent partial waves use the recurrence. |
| wave overlap and momentum routines | matching L/beta | Representation compatibility or an exact same-basis shortcut, not a suppressed interaction. |
| `spectrum.jl`, `flavor_mixing.jl` | equal quantum labels | State lookup, grouping, and orthogonality. |

The small-r guards in derivative kernels are numerical limit handling, not
orbital exclusions. This review does not independently rederive those kernels
or certify the remaining published mixing discrepancies as solved.

## Regression design

- Prescribed n=0 Gaussian/HO waves for S, P, D, F and G: independent closed-form
  smeared-Gaussian integral, nonzero contact, singlet/triplet ratio -3.
- Dressed contact in both native representations for the same orbital range.
- Contact-only diagnostic versus production solve for every tested L and both
  solvers.
- Hellmann-Feynman finite difference of the assembled energy versus its reported
  contact expectation, on fixed grids/bases.
- Existing independent P-wave contact reference and full package suites.

These checks address operator completeness and consistency in addition to
numerical convergence. No test count by itself establishes physical correctness.

## Validation results

- Full suites passed: GIModel 1107/1107, QuarkModelTransitions 598/598,
  GIPaper 2652/2652. The final expanded contact test file was then run
  separately: 88/88, including G waves and local-prescription/grid checks.
- Compared the pre-refactor source at `85e3bdb` with the refactor on 112 fixed
  states: q-qbar, c-ubar, c-cbar and b-bbar, n=1 S/P/D/F, FD and HO.
  The HO comparison used beta bracket 0.15:0.1:4.15 GeV because the old
  bottomonium n=1 solve reached the default upper bracket; both versions used
  the same expanded bracket. Maximum mass, spin-orbit and tensor differences
  were exactly zero; maximum contact-expectation difference was
  5.083e-11 GeV. The latter is from using the same native HO matrix for
  reporting instead of a separate adaptive integral.
- Regenerated Tables V–VII input traces; all calculated outputs were byte
  unchanged. Source fingerprints and 220/79/61 canonical row coverage passed.
- Regenerated all three density caches; values were unchanged, with only the
  four edited source hashes updated. Report figures therefore need no change.
- Existing diffuse-HO quadrature and deliberately under-resolved-grid warnings
  remain visible. This refactor does not change their tolerances.

No parameters were fitted, and no merge into main was performed.
