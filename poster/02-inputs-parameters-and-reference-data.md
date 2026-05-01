# Poster 2: Inputs, Parameters, And Reference Data

## Purpose

Explain the first layer: before solving anything, the model needs trustworthy
inputs. This poster should make the data provenance visible and show why
parameter/reference hygiene is a physics requirement, not bookkeeping.

## Headline

Layer 1: make the paper numerically reproducible before touching the solver.

## One-Sentence Takeaway

The model can only reproduce GI if the same masses, potential constants,
relativization exponents, reference spectra, and state labels enter the code
with documented provenance.

## Main Visual

Draw three input streams merging into one "validated model input" reservoir:

1. Table II parameters.
2. Reference spectrum rows.
3. Formula/research audit notes.

The llama mascot carries labelled crates:

- `m_u, m_s, m_c, m_b`
- `b, c`
- `sigma0, s`
- `epsilon_c, epsilon_t, epsilon_so`

An orca fin appears under the "relativization exponents" crate to foreshadow
operator physics.

## Input Objects

### Parameter File

Code anchor: `data/parameters.provisional.toml`

The active 1.0 model reads:

- quark masses;
- linear confinement parameters `b` and `c`;
- relativistic smearing parameters `sigma0` and `s`;
- Coulomb running parameters from the model implementation;
- Appendix-A central path switches;
- contact and fine-structure momentum switches;
- diagnostic global spin scales.

### Reference Spectra

Code anchors:

- `data/reference_spectrum_charmonium.csv`
- `data/reference_spectrum_bottomonium.csv`
- `data/reference_spectrum_charmed.csv`
- `data/reference_spectrum_b_flavored.csv`
- `data/reference_spectrum_strange.csv`
- `data/reference_spectrum_isovector.csv`
- `data/reference_spectrum_isoscalar.csv`

Every row must identify:

- sector;
- quark content;
- radial quantum number `n`;
- spin multiplicity;
- orbital label `L`;
- total angular momentum `J`;
- reference mass;
- confidence/provenance.

### Research And Formula Audit

Documentation anchors:

- `docs/formula_map.md`
- `docs/research_spin_independent_gi.md`
- `docs/appendix_a_from_paper.md`

These files explain which formulas are active, which are diagnostic branches,
and which paper/web references justified the implementation path.

## Visual Table: What Each Input Controls

| input | physical meaning | downstream effect |
|---|---|---|
| quark masses | constituent mass scale | kinetic energy, smearing width, thresholds |
| `b` | linear confinement slope | radial/orbital level spacing |
| `c` | constant energy offset | absolute mass scale |
| `sigma0`, `s` | GI smearing width | softened singularities and Appendix-A kernels |
| `epsilon_c` | contact relativization | S-wave hyperfine splitting |
| `epsilon_t` | tensor relativization | triplet fine splittings |
| `epsilon_so_vector` | color-magnetic spin-orbit relativization | vector L dot S |
| `epsilon_so_scalar` | Thomas/scalar spin-orbit relativization | scalar L dot S |

## Key Validation Gates

Show these as stamps on the input streams:

- Table II digitization matches TOML.
- Reference CSVs contain required columns.
- State labels are parsed consistently.
- Formula map identifies code location and equation role.
- Reports are regenerated from scripts, not hand-edited.

## Suggested Infographic Panels

### Panel A: Paper Values Become Machine Values

Visual: a paper table passes through a scanner into `parameters.provisional.toml`.

Caption: "Digitization errors masquerade as physics errors."

### Panel B: State Labels Are Contracts

Visual: one meson row exploded into fields:

```text
2^3P_1 -> n=2, multiplicity=3, L=P, J=1
```

Explain that the solver groups states by `n` and `L`, while spin corrections
depend on multiplicity and `J`.

### Panel C: Active Switches Define The Experiment

Show the current 1.0 active switches:

```toml
appendix_a_momentum_sandwich = true
contact_momentum_sandwich = true
fine_structure_momentum_sandwich = true
fine_structure_smeared_kernels = true
```

Caption: "The scorecard only means something if the active path is explicit."

### Panel D: Provenance Prevents Overfitting

Visual: two paths:

- bad path: residual -> retune knobs;
- good path: residual -> formula audit -> layer implementation -> tests.

## What This Layer Explains

- Why the model starts from specific constants and state rows.
- Why reports include confidence/provenance information.
- Why the formula map is part of the scientific result.
- Why we can say "GI reproduction" rather than "some fitted spectrum".

## What This Layer Does Not Explain Yet

- It does not solve the radial equation.
- It does not decide the correct operator ordering.
- It does not fix missing unequal-mass mixing or annihilation.

## Designer Notes

- Use warm yellow for input/provenance.
- Keep the data stream clean and rectangular.
- Let the llama mascot carry the parameter crates, but keep all numbers readable.
- Add one tiny speech bubble: "Same inputs, or it is not the same test."

