# Poster 4: Appendix-A Smearing And Momentum Sandwiches

## Purpose

Explain the layer that moved the project past the stuck point: GI-style smearing
and momentum-dependent relativization. This poster should make clear why using a
pointwise potential alone was not enough.

## Headline

Layer 3: singular potentials become smeared, momentum-aware operators.

## One-Sentence Takeaway

The 1.0 central breakthrough was replacing the naive pointwise Coulomb path with
closed-form smeared `G~(r)`, smeared `S~(r)`, and a Hermitian momentum sandwich
around the Coulomb operator.

## Main Visual

Draw two paths from the same paper formula:

### Old Diagnostic Path

```text
G(r), S(r) sampled pointwise on the radial grid
```

Visual: jagged sharp potential curve.

### Active 1.0 Path

```text
G(r), S(r)
  -> Gaussian smearing
  -> G~(r), S~(r)
  -> A(p) G~ A(p) + S~
```

Visual: an orca dives through a water layer labelled "smearing", then surfaces
with a smooth curve and a sandwich diagram.

## Core Concepts

### Smearing Width

GI smearing depends on the quark masses. The code uses:

```text
sigma = sigma(m1, m2; sigma0, s)
```

Physical role:

- heavier quarks smear differently than lighter quarks;
- singular short-distance behavior becomes finite;
- contact and fine-structure operators become numerically meaningful.

### Smeared Coulomb And Confinement

Show compact symbolic forms:

```text
G(r)  ->  G~(r)
S(r)  ->  S~(r)
```

Use a visual caption:

"The tilde means the paper's short-distance relativistic smearing has been
applied."

### Momentum Sandwich

The active central Coulomb operator is:

```text
G' = A(p) G~ A(p)
```

where `A(p)` is a momentum-dependent relativization factor built on the
finite-difference `p^2` eigenbasis.

Visual: a literal sandwich:

```text
A(p)  |  G~(r)  |  A(p)
```

Caption: "Hermitian by construction: the same factor appears on both sides."

## Code Anchors

- `smeared_coulomb_G_closed(...)`
- `smeared_confinement_S_closed(...)`
- `appendix_a_closed_central_values(...)`
- `momentum_relativization_matrix(...)`
- `appendix_a_momentum_sandwich_matrix(...)`
- `central_potential_mode(...)`
- `central_potential_path(...)`

These are primarily in `src/GIModel.jl` and
`src/appendix_a_status.jl`.

## Why This Layer Was Necessary

Before this layer, residuals were dominated by a central-path mismatch. Spin
corrections could not repair a wrong spin-independent backbone. The central
operator determines:

- spin-averaged multiplet centers;
- radial spacings like `2S - 1S`;
- orbital spacings like `1P - 1S`;
- wavefunctions used by later spin corrections.

Once this layer was active, heavy-quarkonium centers and spacings moved to the
few-MeV regime, making fine-structure details meaningful.

## Visual Comparison: Pointwise vs Active

Use a three-column table:

| central path | visual metaphor | purpose |
|---|---|---|
| pointwise `G + S` | sharp cliff | baseline intuition only |
| diagonal `G~ + S~` | smoothed slope | isolate smearing effect |
| `A(p) G~ A(p) + S~` | smooth slope with momentum bridge | active 1.0 path |

## Active Switches

Show:

```toml
[potential]
appendix_a_momentum_sandwich = true
```

Also show the precedence:

```text
appendix_a_momentum_sandwich
  > appendix_a_closed_form
  > appendix_a_derivative_g
  > appendix_a_smearing
  > coulomb_1d_smear
  > pointwise
```

Caption: "The dispatcher prevents ambiguous central physics."

## What This Layer Explains

- Why the central potential had to be more than pointwise `G(r)+S(r)`.
- Why external GI-family formulas were useful for moving beyond the paper-only
  stuck point.
- Why the heavy-quarkonium radial/orbital spacings became accurate.
- Why operator ordering and Hermiticity matter.

## What This Layer Does Not Explain Yet

- It does not split S-wave singlet/triplet states.
- It does not split P-wave triplets by `J`.
- It does not implement light-sector annihilation or flavor mixing.
- It does not prove identity with the original HO-basis implementation order;
  it implements the corresponding operator path on the finite-difference basis.

## Suggested Infographic Panels

### Panel A: The Problem With Sharp Potentials

Visual: a llama stops at a cliff labelled "singular short-distance behavior".

### Panel B: Smearing

Visual: an orca passes the cliff through a smoothing wave, producing `G~` and
`S~`.

### Panel C: Momentum Sandwich

Visual: `A(p)` slices around `G~`, with a "Hermitian" stamp.

### Panel D: Why Results Improve

Visual: multiplet centers snapping into alignment before spin splittings are
added.

## Designer Notes

- Use ocean blue for smearing/momentum effects.
- Make the sandwich diagram instantly recognizable but not silly.
- Keep `G~` and `S~` visually distinct from unsmeared `G` and `S`.

