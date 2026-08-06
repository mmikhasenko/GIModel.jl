# Harmonic-oscillator differentiation route

The harmonic-oscillator (HO) route is the paper-style Appendix-A Eq. (A17)
calculation and is a peer of the finite-difference route. Its central operator is
not an FD matrix projected onto an oscillator basis: for the active Appendix-A
path, both momentum- and position-space matrix elements are constructed in the
finite HO space.

## Primal contract

At a fixed oscillator scale `beta` and basis size `nbasis`, the route is:

```text
exact HO p² matrix
    -> prepared p² eigensystem
generalized Gauss-Laguerre nodes/projector
    -> G and S operator matrices
    -> dense H_HO(theta; beta)
    -> symmetric eigenpairs
    -> central masses for this beta branch
```

The uniform spatial mesh in `oscillator_hamiltonian_for_beta` is only used to
reconstruct wavefunctions for reporting. It is not part of the active
Appendix-A central Hamiltonian.

## Exact momentum-space matrix

For radial basis index `n=0,...,nbasis-1`,

```math
\langle n|p^2|n\rangle
= \beta^2(2n+L+3/2),
```

```math
\langle n|p^2|n+1\rangle
= \beta^2\sqrt{(n+1)(n+L+3/2)}.
```

At fixed beta, this matrix and its eigensystem are constant with respect to the
physical model parameters. Kinetic and momentum-factor matrices are spectral
functions of it, just as on the FD route.

## Position-space operator matrix

For generalized Gauss--Laguerre/DVR nodes, the code diagonalizes the Jacobi
matrix represented by `beta² * ho_r2_matrix`. Let its eigenvalues be `x_i` and
the relevant eigenvector rows be `Z_ai`. Then

```math
r_i=\frac{\sqrt{x_i}}{\beta},
\qquad
\langle a|g(r)|b\rangle
=\sum_i Z_{ai}g(r_i)Z_{bi}.
```

For the active route, `g` is `G_tilde` or `S_tilde`. The central matrix has the
same abstract structure as the FD route,

```math
H_{\mathrm{HO}}(\theta;\beta)
=T + A\,\widetilde G_{\mathrm{HO}}A
+\widetilde S_{\mathrm{HO}},
```

but all matrices have size `nbasis x nbasis` and the position operators are
quadrature matrices, not mesh diagonals.

For physical-parameter derivatives at fixed beta, the Jacobi eigensystem and
nodes can be prepared. Only `g(r_i;theta)` and the mass-dependent spectral
functions remain active.

## Two distinct derivative semantics

### Fixed-beta branch

```math
E_i(\theta;\beta)
```

is smooth away from eigenvalue degeneracies. This is the fundamental HO
derivative kernel and is what `ho_fixed_beta_probe.jl` checks.

### Beta-selected channel result

The current channel solver scans

```text
beta = 0.25, 0.35, ..., 2.35
```

and chooses the beta minimizing the `nlevels`-th eigenvalue. Define

```math
\beta_*(\theta)=\arg\min_{\beta\in\mathcal B}
E_{n_\mathrm{levels}}(\theta;\beta).
```

The reported result is `E_i(theta; beta_*(theta))`. On an open region where one
grid point is the unique winner, `beta_*` is constant and the derivative equals
the fixed-beta branch derivative. At a switch between grid points the result is
generally nondifferentiable. Moreover, beta is selected using only the
`nlevels`-th state: at a tie in that selection score, the other reported levels
on the two finite-basis branches need not agree and can jump when the selected
branch changes.

Therefore beta selection should be evaluated before AD and return:

- selected beta;
- best and second-best selection scores;
- the score margin;
- optionally, stability of the selected beta under the finite-difference steps
  used by derivative tests.

The local derivative is supported only while the selection margin is safely
nonzero. This is not the same contract as treating beta as a continuously
optimized or fitted parameter.

The selected beta also depends on `nlevels`, because the last requested level is
the selection score. `nlevels` is therefore part of the discrete prepared
problem identity.

## Adaptive quadrature

`ho_operator_matrix` doubles `nq` until successive matrices meet a tolerance.
The mathematical integral is smooth, but the implementation changes work and
matrix nodes when the accepted `nq` changes.

For a robust derivative kernel:

1. determine a converged `nq` during preparation for the selected beta, `L`, and
   basis size;
2. prepare the Jacobi eigensystem/nodes at that fixed order;
3. evaluate active scalar functions on those fixed nodes during AD;
4. verify after perturbation that the chosen order remains converged.

This removes both repeated constant LAPACK work and value-dependent adaptive
control flow from the differentiated region.

## Wave reconstruction and phase

Central masses do not require `ho_basis_matrix`, mesh QR orthonormalization, or
wave reconstruction. Those operations should not be in the fixed-beta
eigenvalue kernel.

Wavefunction-dependent observables need a later HO-specific contract. The QR and
eigenvector phase fixes are discontinuous when an anchor crosses zero, although
phase-invariant expectations remain well-defined. First-order distorted states
also contain explicit energy denominators and therefore require gap handling.

## Probe evidence

Run:

```bash
julia --project=. diff_support/probes/ho_fixed_beta_probe.jl
```

The probe uses charmonium, `L=0`, `beta=0.65`, and `nbasis=8`. The small basis is
chosen for a quick structural test, not production accuracy. It independently
checks `b`, `c`, `sigma0`, `s`, and the common charm mass.

At relative step `1e-4`, fixed-beta direct eigenvalue differences and
Hellmann--Feynman derivatives agree to about `5.5e-8` or better for the first
three levels. It also finds `max(abs(dH/dc-I)) = 3.8e-11`.

This proves the smooth fixed-beta branch is numerically straightforward. It does
not yet validate the beta-selection margin or an AD backend. The probe calls the
current adaptive `ho_operator_matrix`; it also does not replace the adaptive loop
with the proposed prepared fixed-order quadrature.

## HO-specific implementation shape

A possible split is:

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

struct HOSelection
    selected
    best_score
    second_score
    margin
end
```

`HOSelection` is discrete orchestration. `PreparedHOBranch` is constant context
for the differentiable fixed-beta matrix assembly.

## HO acceptance gate

The HO route is supported only when:

- fixed-beta prepared assembly reproduces today's HO matrix/eigenvalues;
- `p²` and quadrature eigensystems are not recomputed in active evaluations;
- fixed-order quadrature remains converged over tested perturbations;
- analytic, finite-difference, JVP, and VJP checks agree on a fixed branch;
- beta selection exposes a margin and flags switch boundaries;
- the selected-branch derivative agrees with end-to-end finite differences when
  those differences keep the same beta;
- several beta values, orbital channels, basis sizes, and mass sectors are tested;
- wave/phase differentiation is explicitly excluded until separately supported.

Passing this gate says nothing by itself about FD support.
