# Differentiable state-energy readiness

## Question

How close is GIModel to providing the reverse-mode derivative

```text
d(state_energy) / d(parameters)
```

through its radial solution?

This is an architecture and numerical-readiness investigation. It does not add
an automatic-differentiation dependency or claim that the current
`compute_spectrum` workflow is differentiable.

## Conclusion

An isolated energy from a fixed `(L,S,J)` harmonic-oscillator calculation is
close to differentiable. The correct reverse rule differentiates the converged
eigenvalue problem, not the iterations used to find it. For

```math
H(\theta)\psi_n = E_n\psi_n, \qquad \psi_n^\dagger\psi_n=1,
```

the pullback of a scalar energy is

```math
\overline H = \overline E_n\,\psi_n\psi_n^\dagger,
```

equivalently the Hellmann–Feynman derivative

```math
\frac{\partial E_n}{\partial\theta_a}
=
\psi_n^\dagger\frac{\partial H}{\partial\theta_a}\psi_n.
```

This makes reverse-mode cost essentially independent of the number of fitted
parameters after one state has been solved.

GIModel already assembles central, contact, symmetric spin-orbit, and diagonal
tensor terms into one fixed-sector Hamiltonian before diagonalization. The same
eigenvalue rule therefore applies to the complete fixed-sector energy, not only
to its central part.

## Numerical probe

The companion script evaluates the native HO charmonium `1S` central
Hamiltonian at fixed `beta = 0.65 GeV` and `nbasis = 12`. It differentiates the
ground-state energy with respect to the confinement slope `b` in two independent
orders:

1. finite-difference the final eigenvalue;
2. finite-difference the Hamiltonian and contract it with the normalized
   eigenvector, `psi' * (dH/db) * psi`.

The observed values were

| quantity | value |
|---|---:|
| `E_1S` | `3.0656706422771114 GeV` |
| finite-difference `dE/db` | `1.4965695061208704 GeV^-1` |
| Hellmann–Feynman `dE/db` | `1.4965695059344297 GeV^-1` |
| relative difference | `1.25e-10` |

This checks the proposed eigensolver boundary using GIModel's actual
Hamiltonian. It is not yet a reverse-AD test of Hamiltonian assembly.

Regenerate it from the repository root with:

```sh
julia GIPaper/scripts/investigate_state_energy_differentiability.jl
```

## Boundary of the differentiable calculation

The adaptive numerical workflow should prepare a calculation, then the reverse
pass should operate on a frozen numerical problem:

```text
parameters
    -> converged numerical preparation (beta, basis size, quadrature order)
    -> pure H(parameters; prepared channel)
    -> selected isolated eigenvalue
    -> state energy
```

The preparation stage may use searches, caches, convergence tests, and integer
choices. These are solver decisions rather than model parameters. Differentiating
through them would make a derivative depend on branch choices such as which beta
grid cell won or which basis size first passed a tolerance.

The beta search is especially important. At infinite basis size a physical
energy is beta-independent. At finite size, beta dependence is numerical error;
the gradient must therefore receive its own convergence check under nearby beta
and increasing basis size. The current search minimizes the highest requested
level, so the envelope theorem does not by itself justify differentiating every
lower state through that selected beta.

## Existing blockers

- `GIParameters` and `ConstituentMasses` preserve generic real-number types and
  are suitable inputs for a differentiable kernel.
- The HO matrix construction is compact and is the natural first backend.
- Adaptive beta refinement, adaptive quadrature order, basis convergence,
  caching, sorting, and state assignment are intentionally non-smooth and should
  remain in the preparation layer.
- `ChannelRadialSolution`, `OscillatorWave`, `MixingBlock`, `MixingResult`, and
  several spectrum result types store or explicitly collect `Float64`. A scalar
  AD path must return the energy before crossing those reporting containers, or
  those containers must later become parametric.
- A named eigenvalue is differentiable only while it is isolated. At a level
  crossing, state identity must be tracked by overlap and the derivative should
  be formulated for a subspace/projector when degeneracy is physical.
- Final same-`J`, tensor, flavor, and radial mixing adds another eigenproblem.
  Its eigenvalue pullback is straightforward, but its matrix elements depend on
  radial waves and therefore require eigenvector-response derivatives from the
  fixed-sector solves.
- Mass-dependent annihilation poles require a separate implicit derivative of
  the fixed-point equation.

## Proposed implementation sequence

1. Introduce a prepared HO channel containing fixed beta, basis size,
   quadrature data, and parameter-independent momentum data.
2. Add a pure scalar `state_energy(parameter_vector, prepared_channel)` path
   which does not construct reporting objects.
3. Add a backend-neutral `ChainRulesCore.rrule` for the selected symmetric
   eigenvalue, returning `psi * psi'` as the Hamiltonian cotangent.
4. Validate reverse gradients for `b`, `c`, `sigma0`, `s`, the relativistic
   epsilon factors, and constituent masses against finite differences.
5. Require gradient convergence under basis, beta, and quadrature refinement.
6. Extend the same boundary to small mixing matrices, with explicit overlap
   tracking and degeneracy tests.

The first three steps are a focused GIModel feature, not a solver rewrite. A
fully differentiable `compute_spectrum` is neither required nor desirable for
the first useful parameter gradient.
