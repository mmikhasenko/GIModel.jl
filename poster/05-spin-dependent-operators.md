# Poster 5: Spin-Dependent Operators

## Purpose

Explain the final physical layer before results: contact hyperfine splitting,
spin-orbit splitting, and tensor splitting. This is the poster that shows how
the solver moves from spin-averaged levels to the detailed multiplet pattern.

## Headline

Layer 4: spin operators turn central levels into individual meson states.

## One-Sentence Takeaway

After the central Hamiltonian sets the multiplet centers, contact hyperfine,
spin-orbit, and tensor operators add state-dependent shifts controlled by spin,
orbital angular momentum, total `J`, smearing, and momentum relativization.

## Main Visual

Draw one central `nL` level splitting into multiple states:

```text
central 1P level
     |
     +-- 1^1P_1
     +-- 1^3P_0
     +-- 1^3P_1
     +-- 1^3P_2
```

Use coral arrows for spin shifts. An orca circles the tensor/spin-orbit arrows;
a llama holds the angular-momentum labels.

## Spin Layers

### Contact Hyperfine

Role:

- affects S-waves;
- splits singlets and triplets;
- uses smeared short-distance contact operator.

Visual:

```text
1S central level
   -> 1^1S_0 shifted down
   -> 1^3S_1 shifted up
```

Active path:

```toml
contact_momentum_sandwich = true
```

Meaning:

```text
contact operator -> B(p) contact B(p)
```

### Spin-Orbit

Role:

- affects states with `L > 0`;
- depends on `L dot S`;
- combines color-magnetic vector term and Thomas-precession term.

Active radial kernels:

```text
vector spin-orbit kernel ~ (1/r) dG~/dr
Thomas kernel           ~ (1/2r) d(G~ + S~)/dr
```

Active path:

```toml
fine_structure_momentum_sandwich = true
fine_structure_smeared_kernels = true
```

### Tensor

Role:

- affects triplet states with `L > 0`;
- redistributes mass inside P-, D-, and F-wave triplets;
- depends on angular tensor factors for each `L, J`.

Active radial kernel:

```text
K_tensor(r) = (1/r) dG~/dr - d^2G~/dr^2
```

This replaced the older diagnostic pointwise running-Coulomb tensor kernel in
the active 1.0 path.

## Code Anchors

- `contact_smearing_sigma(...)`
- `contact_hyperfine_shift_active(...)`
- `contact_hyperfine_shift_momentum_sandwich(...)`
- `fine_structure_components(...)`
- `radial_expect_momentum_sandwich(...)`
- `smeared_coulomb_G_prime_closed(...)`
- `smeared_coulomb_G_second_closed(...)`
- `smeared_confinement_S_prime_closed(...)`
- `tensor_kernel_smeared_coulomb(...)`
- `LdotS(...)`
- `tensor_triplet_LJ(...)`

Primary file: `src/GIModel/spin_fine_structure.jl`.

## Formula Panel

Show the predicted mass as a stacked sum:

```text
M_predicted =
  M_central
  + Delta_contact
  + Delta_spin-orbit
  + Delta_tensor
```

Then map each term:

| shift | depends on | mainly visible in |
|---|---|---|
| contact | spin singlet/triplet, S-wave wavefunction at short distance | S-wave hyperfine gaps |
| vector spin-orbit | `L dot S`, `G~'(r)` | P/D/F triplet ordering |
| Thomas spin-orbit | `L dot S`, `(G~+S~)'(r)` | cancellation/rebalance of spin-orbit |
| tensor | tensor angular factor, `G~'(r)`, `G~''(r)` | detailed triplet splittings |

## Angular Momentum Mini-Panel

Explain with a compact visual:

```text
S = 0 singlet: no triplet fine-structure tensor/L dot S split
S = 1 triplet: J = L-1, L, L+1 can split
```

For P-waves:

```text
L=1, S=1 -> J=0, 1, 2
```

This is why the model predicts separate `1^3P_0`, `1^3P_1`, and `1^3P_2`
masses.

## Why Momentum Sandwiches Matter Here Too

The GI paper does not only smear coordinate-space singularities; it also
relativizes spin operators with momentum-dependent factors. In the code this is
implemented as a Hermitian expectation:

```text
<u | B(p) K(r) B(p) | u>
```

where `K(r)` is the radial kernel for contact, tensor, or spin-orbit, and the
exponent depends on the relevant `epsilon` parameter.

Visual: same sandwich motif as Poster 4, but coral instead of blue.

## What This Layer Explains

- Why S-wave singlets and triplets split.
- Why P-, D-, and F-wave triplets separate by `J`.
- Why smeared derivatives of `G~` and `S~` were the natural next step after the
  central Appendix-A path.
- Why the charmonium P-wave residuals became much better after implementing
  smeared derivative kernels.

## What This Layer Does Not Explain Yet

- Unequal-mass antisymmetric spin-orbit terms.
- Same-`J` tensor mixing between basis states.
- Light isoscalar annihilation and flavor mixing.
- Coupled-channel mass shifts.

## Suggested Infographic Panels

### Panel A: From One Central Level To Many States

Visual: one thick teal line splits into coral fine-structure lines.

### Panel B: Contact Hyperfine

Visual: S-wave pair moving apart.

### Panel C: Spin-Orbit And Tensor

Visual: P-wave triplet fan with arrows labelled vector, Thomas, tensor.

### Panel D: Smeared Derivative Kernels

Visual: smooth `G~` curve with tangent and curvature markers:

- first derivative for spin-orbit;
- first and second derivative for tensor.

## Designer Notes

- Use coral for spin shifts.
- Use blue outlines when a spin kernel comes from smeared `G~`/`S~` to show
  continuity with Poster 4.
- Use one orca dive loop around `G~'` and `G~''` to suggest "derivatives below
  the surface."

