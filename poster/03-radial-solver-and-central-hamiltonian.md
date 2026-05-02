# Poster 3: Radial Solver And Central Hamiltonian

## Purpose

Explain the computational backbone: the finite-difference radial Hamiltonian that
produces spin-independent eigenvalues and eigenvectors. Later layers only make
sense once this base machine is visible.

## Headline

Layer 2: the radial Hamiltonian builds the spin-independent spectrum.

## One-Sentence Takeaway

The solver turns quark masses and central potentials into radial wavefunctions
`u(r)`, and those same wavefunctions become the stage for contact, tensor, and
spin-orbit corrections.

## Main Visual

Draw a radial line `r = 0 ... rmax` as a mountain trail. A llama walks along the
grid points carrying a kinetic-energy compass. Under the trail, an orca shadow
marks the momentum-space operator `p^2`.

Show the Hamiltonian as a machine:

```text
quark masses + L + radial grid
        |
        v
T = sqrt(p^2 + m1^2) + sqrt(p^2 + m2^2)
Vcentral(r)
        |
        v
H u_nL(r) = E_nL u_nL(r)
```

## Core Equation

Use a compact callout:

```text
H = sqrt(p^2 + m1^2) + sqrt(p^2 + m2^2) + Vcentral
```

For nonrelativistic diagnostic paths:

```text
H_NR = m1 + m2 + p^2/(2 mu) + Vcentral
```

But the 1.0 reproduction path uses the relativistic kinetic operator.

## Code Anchors

- `radial_grid(...)` creates a uniform radial mesh.
- `p2_operator(...)` builds the radial momentum-squared operator.
- `relativistic_hamiltonian(...)` builds the active heavy-quarkonium Hamiltonian.
- `solve_sector(...)` tabulates equal-mass eigenvalues per `(n, L)` for one flavor.
- `compute_sector(...)` (**`src/sector_comparison.jl`**) caches spin-independent radial solves per `RadialChannelKey`
  (`ConstituentMasses` + orbital letter `L`) built from **`ReferenceStateWithMasses`** rows.
- `compare(computed, annotated)` (**`src/sector_comparison.jl`**) uses the same rows, maps levels to reference labels,
  applies contact / fine-structure shifts, and returns residual **`NamedTuple`** rows; **`write_residual_report`**
  is markdown-only on those rows. See **`docs/code_architecture.md`**.

## Important Convention

The radial eigenvector is treated as a reduced radial wavefunction:

```text
integral |u(r)|^2 dr = 1
```

There is no extra `4 pi` factor in radial expectation values because the angular
spherical harmonics are already normalized.

Make this visual:

```text
full wavefunction = u(r)/r * Y_LM(theta, phi)
angular integral = 1
radial integral = integral |u(r)|^2 dr
```

## Central Potential Before GI Smearing

Introduce the simple central form:

```text
V(r) = G(r) + S(r)
G(r) = -4 alpha_s(r) / (3r)
S(r) = b r + c
```

This is the intuitive backbone, but not the final 1.0 central operator. Poster 4
explains how `G(r)` and `S(r)` become smeared and momentum-dependent.

## Visual Layer Stack

Show the base layer as:

```text
radial grid
  + p^2 operator
  + relativistic kinetic energy
  + central potential
  = spin-independent eigenstates
```

Then draw faint placeholders above it:

- contact shift;
- spin-orbit shift;
- tensor shift.

Caption: "Spin corrections do not replace the radial solver; they sit on top of
its wavefunctions."

## Diagnostic Outputs

The central solver supplies the baseline pieces later shown in reports:

- `central_GeV`
- `contact_shift_GeV`
- `spin_orbit_shift_GeV`
- `tensor_shift_GeV`
- `predicted_GeV`
- `residual_MeV`

Show one example report row as an exploded bar:

```text
central mass
 + contact
 + spin-orbit
 + tensor
 = predicted mass
```

## What This Layer Explains

- How a meson state receives a base mass.
- Why `n` and `L` are solver labels.
- Why the same wavefunction must be normalized consistently for all later
  expectation values.
- Why central spin-independent physics dominates radial/orbital spacings.

## What This Layer Does Not Explain Yet

- It does not include the full Appendix-A central relativization by itself.
- It does not split spin multiplets.
- It does not solve unequal-mass spin mixing.

## Active 1.0 Badge

Place this on the poster:

```text
Central Hamiltonian active path:
relativistic kinetic + Appendix-A momentum-sandwiched central potential
```

Then point to Poster 4 for the detailed `A(p) G~ A(p)` explanation.

## Designer Notes

- Use deep teal for the Hamiltonian machine.
- Use evenly spaced tick marks for the radial grid.
- Use the llama as a guide along the mesh, not as decoration over equations.
- Use the orca subtly to suggest "momentum-space depth" under the radial line.

