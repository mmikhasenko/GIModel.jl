# Finite-difference differentiation route

The finite-difference (FD) route is one of two peer computation routes. It has
its own preparation object, numerical risks, derivative probe, and acceptance
gate. It must not be used as an implementation surrogate for the harmonic-
oscillator route.

## Primal contract

For fixed orbital `L` and solver settings `(ngrid, rmax, nlevels)`, the route is:

```text
uniform interior grid
    -> FD p² tridiagonal matrix
    -> prepared p² eigensystem
    -> dense relativistic H_FD(theta)
    -> lowest symmetric eigenpairs
    -> selected central masses
```

The active Appendix-A Hamiltonian is

```math
H_{\mathrm{FD}} = T(P^2;m_1,m_2)
+ A(P^2;m_1,m_2)\,\operatorname{diag}(\widetilde G(r_i))\,
  A(P^2;m_1,m_2)
+ \operatorname{diag}(\widetilde S(r_i)).
```

The detailed scalar and spectral equations are in
[`../hamiltonian_derivatives.md`](../hamiltonian_derivatives.md).

## What is fixed and what is active

Fixed context:

- `L`, `ngrid`, `rmax`, `nlevels`;
- the radial grid and spacing;
- the FD `p²` matrix and its eigensystem;
- the central-potential method and eigensolver choice.

Initially active:

- `b`, `c`, `sigma0`, `s`;
- the two constituent masses.

Later spin parameters require a separate extension of the contract.

## FD-specific preparation

A possible internal shape is:

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

The exact fields should be concretely typed. `p2_operator` currently accepts a
mass argument but does not use it. The probe verifies exactly that
`p2(m) == p2(2m)`, so the eigensystem can be moved out of the active trace.

The existing under-resolution diagnostic should run during preparation or
validation, not during AD. A derivative of a converged discrete result is still
unhelpful if the discrete result itself is not converged.

## FD-specific mutation

The route uses local writes to assemble tridiagonal diagonals and normalize
newly allocated wave matrices. These are not persistent state mutations.

For central eigenvalues:

- keep local assembly mutation for Enzyme/Mooncake;
- make storage element types generic for ForwardDiff;
- omit wave normalization entirely from an eigenvalue-only kernel;
- keep `SectorComputation` dictionaries outside the route.

The dense `eigen(H)` call is the spectral boundary. The project should own its
eigenvalue JVP/VJP rather than requiring a backend to differentiate LAPACK.

## FD-specific numerical issues

1. **Grid convergence.** Values and derivatives need separate `ngrid/rmax`
   convergence tests.
2. **Dense cost.** The default `N=450` Hamiltonian is much larger than the HO
   matrix. Forward-mode chunks may repeat expensive dense work.
3. **Krylov path.** The differentiated route should initially use the full
   symmetric solver. Differentiating `eigsolve` is a separate algorithmic task.
4. **Small `p²` clamps.** Prepared `p²` eigenvalues should be checked once; a
   `max(lambda,0)` kink should not silently become part of the derivative model.
5. **Mass rounding.** Rounding in `ConstituentMasses` must be removed from the
   physics value before AD support.

There is no beta selection or adaptive position-space quadrature in the active
FD central route.

## Probe evidence

Run:

```bash
julia --project=. diff_support/probes/central_hamiltonian_probe.jl
```

The probe uses charmonium, `L=0`, and `ngrid=160` for speed. It verifies:

- `p²` independence from the mass argument;
- Hellmann--Feynman derivatives for `b`, `c`, `sigma0`, `s`, and common charm
  mass against direct eigenvalue finite differences;
- `dH/dc = I`;
- the analytic matrix `dH/db`.

At relative step `1e-4`, the level-1 spectral comparisons agree to `1.5e-8` or
better. This is a route-structure check, not a production spectrum benchmark.

## FD acceptance gate

The FD route is supported only when:

- a prepared FD problem reproduces today's FD central masses on the same grid;
- repeated evaluations do not recompute the `p²` eigensystem;
- analytic, finite-difference, JVP, and VJP checks agree;
- derivatives converge under grid refinement;
- the API reports/monitors small eigenvalue gaps;
- no spectrum cache, report object, or state-assignment mutation is inside AD;
- supported backend versions and performance are recorded.

Passing this gate says nothing by itself about HO support.
