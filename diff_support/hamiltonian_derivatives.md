# Shared physics and spectral derivative equations

This note records physics functions and spectral rules shared by two separately
implemented computation routes:

- [`routes/finite_difference.md`](routes/finite_difference.md)
- [`routes/harmonic_oscillator.md`](routes/harmonic_oscillator.md)

The matrix representation below is explicitly labelled by route. Symbols follow
the source code rather than introducing a separate mathematical API.

## FD radial representation

For `ngrid = N`, `rmax = R`, and orbital angular momentum `L`,

```math
h = \frac{R}{N+1}, \qquad r_i = ih, \quad i=1,\ldots,N.
```

The finite-difference representation of radial `p²` is symmetric tridiagonal:

```math
(P^2)_{ii} = \frac{2}{h^2} + \frac{L(L+1)}{r_i^2},
\qquad
(P^2)_{i,i+1}=(P^2)_{i+1,i}=-\frac{1}{h^2}.
```

Write its eigendecomposition as

```math
P^2 = Q\,\mathrm{diag}(\lambda_a)Q^T.
```

For a fixed numerical problem, `r`, `h`, `P²`, `Q`, and `lambda` are independent
of `b`, `c`, `sigma0`, `s`, and the quark masses. They should be constant members
of a prepared workspace.

## Smearing width and central functions

The code implements Appendix A's mass-dependent width

```math
\sigma^2 = \sigma_0^2\left[\frac12 + \frac12
\left(\frac{4m_1m_2}{(m_1+m_2)^2}\right)^4\right]
+ s^2\left(\frac{2m_1m_2}{m_1+m_2}\right)^2.
```

For the Gaussian components of the running coupling,

```math
\tau_k^{-2}=\sigma^{-2}+\gamma_k^{-2}.
```

The closed-form smeared Coulomb function is

```math
\widetilde G(r) = -\sum_k \frac{4\alpha_k}{3r}
\operatorname{erf}(\tau_k r),
```

with its analytic `r -> 0` limit used in the implementation. The smeared linear
plus constant confinement is

```math
\widetilde S(r) = br\left[
\frac{e^{-\sigma^2r^2}}{\sqrt{\pi}\,\sigma r}
+ \left(1+\frac{1}{2\sigma^2r^2}\right)\operatorname{erf}(\sigma r)
\right] + c,
```

again with an analytic origin limit. On the strictly positive FD interior grid,
the special origin branch is not normally taken.

## FD active relativistic Hamiltonian

Define spectral one-particle energies

```math
E_{ja}=\sqrt{\lambda_a+m_j^2}, \qquad j\in\{1,2\}.
```

The kinetic matrix is

```math
T = Q\,\mathrm{diag}(E_{1a}+E_{2a})Q^T.
```

The central relativistic momentum factor is

```math
A = Q\,\mathrm{diag}\left(
\sqrt{1+\frac{\lambda_a}{E_{1a}E_{2a}}}
\right)Q^T.
```

With `G = diag(G_tilde(r_i))` and `S = diag(S_tilde(r_i))`, the active
Hamiltonian is

```math
H(\theta)=T + AGA + S.
```

This is exactly the FD structure assembled by
`appendix_a_momentum_sandwich_matrix` and `relativistic_hamiltonian`.

All parameter-dependent operations in this expression are scalar elementary
functions, diagonal construction, and dense matrix multiplication. This makes
FD Hamiltonian assembly a plausible AD target once the constant `p²`
eigendecomposition is moved out of the trace.

## HO fixed-beta Hamiltonian

At fixed beta, `L`, basis size, and quadrature order, the HO route constructs an
exact finite-basis `p²`, Gauss--Laguerre matrices for `G_tilde` and `S_tilde`,
and then

```math
H_{\mathrm{HO}}(\theta;\beta)
=T(P^2_{\mathrm{HO}};m_1,m_2)
+A(P^2_{\mathrm{HO}};m_1,m_2)\widetilde G_{\mathrm{HO}}A(P^2_{\mathrm{HO}};m_1,m_2)
+\widetilde S_{\mathrm{HO}}.
```

This has the same parameter-dependent spectral and matrix-product structure as
FD, but it is a different finite representation. Its `p²` eigensystem and
quadrature projector can be prepared at fixed beta. The detailed HO equations,
beta-selection semantics, and acceptance gate are in
[`routes/harmonic_oscillator.md`](routes/harmonic_oscillator.md).

## Useful analytic checks

### Confinement offset

On both FD and a correctly normalized HO branch, the offset appears only as
`cI`, hence

```math
\frac{\partial H}{\partial c}=I,
\qquad
\frac{\partial E_i}{\partial c}=1.
```

The probes find `max(abs(dH/dc - I)) = 1.4e-10` for FD and `3.8e-11` for the
fixed-beta HO branch, limited by the finite-difference calculations used in the
probes.

### Confinement slope

`b` enters linearly only in `S_tilde`, while `sigma` is independent of `b`. On
the FD route:

```math
\frac{\partial H}{\partial b}
= \operatorname{diag}\left(r_i B(\sigma r_i)\right),
```

where the bracket `B` is the bracketed expression in `S_tilde`. This is another
simple rule/oracle for testing AD. The probe compares this analytic matrix to a
finite-difference `dH/db` directly.

### Kinetic mass derivative

Because `Q` and `lambda` are fixed,

```math
\frac{\partial}{\partial m_j}
Q\,\mathrm{diag}(E_{ja})Q^T
= Q\,\mathrm{diag}\left(\frac{m_j}{E_{ja}}\right)Q^T.
```

Masses also affect `A`, `sigma`, `G_tilde`, and `S_tilde`, so this is a component
check rather than the complete mass derivative.

## Symmetric spectral derivatives

Let

```math
Hv_i=E_iv_i, \qquad v_i^Tv_i=1.
```

For a simple eigenvalue and perturbation `Hdot`, the
Hellmann--Feynman relation gives

```math
\dot E_i=v_i^T\dot H v_i.
```

Equivalently, the matrix cotangent for a scalar objective `ell(E)` is

```math
\overline H = V\,\mathrm{diag}(\overline E) V^T,
\qquad
\overline E_i=\frac{\partial \ell}{\partial E_i}.
```

This eigenvalue-only VJP contains no division by an eigenvalue gap. It is the
preferred primitive for central masses.

For eigenvectors, one common gauge gives

```math
\dot v_i = \sum_{j\ne i} v_j
\frac{v_j^T\dot H v_i}{E_i-E_j}.
```

This exposes the main difficulty of differentiating wavefunction-dependent spin
corrections: small gaps amplify the derivative, and exact degeneracy makes an
individual eigenvector undefined. A phase/sign convention does not solve the
degenerate-subspace ambiguity.

At an exact `k`-fold degeneracy, the first-order eigenvalue shifts are the
eigenvalues of the perturbation projected into the degenerate subspace. An API
that promises derivatives for individually labelled states must therefore state
how it handles degeneracy and crossings.

## Consequences for later spectrum stages

### First-order spin corrections

The present fine-structure code evaluates expectations on central eigenvectors.
For an operator `O(theta)`,

```math
\frac{d}{d\theta}\langle v|O|v\rangle
= \langle v|\dot O|v\rangle
+ 2\,\operatorname{Re}\langle \dot v|O|v\rangle.
```

Thus differentiating corrected masses requires eigenvector response, not only
the simple eigenvalue rule. Options to investigate later are:

1. use the explicit gap formula/projected linear solve away from degeneracy;
2. formulate the correction as an eigenvalue of a larger effective Hamiltonian;
3. initially treat central waves as stopped/constant and clearly label the
   result as a partial derivative (not recommended as the default contract).

### Nonperturbative contact term and small mixing blocks

Where the code already forms `H + V` and diagonalizes it, the same eigenvalue
rule applies directly. The `2x2` and small annihilation/mixing matrices are also
straightforward spectral problems. Their difficulty is state identity at
crossings, not matrix size.

### HO basis

For fixed beta, the HO matrix uses the shared spectral derivative. The current
solver chooses beta by minimizing an eigenvalue over a fixed grid, making the
reported result a locally smooth selected branch and nonsmooth at branch
switches. The HO contract therefore differentiates the selected fixed-beta
branch only when it has a safe margin to the second-best beta. See
[`routes/harmonic_oscillator.md`](routes/harmonic_oscillator.md).

## FD numerical probe result

The table below uses the active central path, charmonium, `L=0`, `ngrid=160`, and
a relative central-difference step of `1e-4`. `FD(E)` differentiates eigenvalues
directly; `HF` evaluates `v' * FD(H) * v` at the unperturbed eigenvector.

| Direction | Level | `FD(E)` | `HF` | absolute difference |
|---|---:|---:|---:|---:|
| `b` | 1 | 1.4871604914 | 1.4871604769 | 1.46e-8 |
| `c` | 1 | 1.0000000001 | 1.0000000000 | 6.85e-11 |
| `sigma0` | 1 | -0.0127150439 | -0.0127150440 | 4.30e-11 |
| `s` | 1 | -0.0290192285 | -0.0290192285 | 7.74e-11 |
| common charm mass | 1 | 1.8405312271 | 1.8405312267 | 3.64e-10 |

The first three central eigenvalues in this probe are `3.06310471`,
`3.66390418`, and `4.08782729` GeV; their smallest adjacent gap is about
`0.424` GeV. The probe is therefore well away from a degeneracy.

This validates the proposed spectral decomposition. It does not yet validate a
specific AD backend or behavior near a crossing.

## HO fixed-beta numerical probe result

The independent HO probe uses charmonium, `L=0`, `beta=0.65`, `nbasis=8`, and a
relative central-difference step of `1e-4`. The first three masses are
`3.06767557`, `3.66947560`, and `4.09831118` GeV, with a smallest adjacent gap of
about `0.429` GeV.

| Direction | Level | `FD(E)` | `HF` | absolute difference |
|---|---:|---:|---:|---:|
| `b` | 1 | 1.5047830522 | 1.5047830375 | 1.47e-8 |
| `c` | 1 | 1.0000000000 | 1.0000000000 | 4.35e-12 |
| `sigma0` | 1 | -0.0112588661 | -0.0112588661 | 5.43e-12 |
| `s` | 1 | -0.0256958299 | -0.0256958299 | 2.26e-11 |
| common charm mass | 1 | 1.8455541517 | 1.8455541513 | 4.18e-10 |

This validates only a smooth fixed-beta HO branch. The selected-beta route also
needs the branch-margin checks defined in
[`routes/harmonic_oscillator.md`](routes/harmonic_oscillator.md).
