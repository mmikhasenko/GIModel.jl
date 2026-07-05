# Poster 1: Introduction And Global Workflow

## Purpose

Explain the whole project in one glance: we are reproducing the
Godfrey-Isgur meson spectrum by translating the paper's physical ingredients
into a transparent computational pipeline. This poster is the map that makes the
other five posters readable.

## Headline

From paper formulas to reproduced meson masses: a layered Godfrey-Isgur
reproduction.

## One-Sentence Takeaway

Version 1.0 succeeds for heavy quarkonium because the solver now combines the
right central Hamiltonian, Appendix-A smearing, momentum-dependent operator
sandwiches, and spin-dependent splittings in the same auditable pipeline.

## Main Visual

Create a left-to-right flow diagram with six stations:

1. Paper and reference data.
2. Parameter loading and provenance checks.
3. Radial Hamiltonian and eigenstates.
4. Appendix-A central smearing and momentum sandwich.
5. Contact, tensor, and spin-orbit shifts.
6. Residual reports and scorecard.

Use the cartoon theme lightly: a llama starts at the paper stack, an orca dives
under the smearing/spin stations, and both appear near the final scorecard.

## Workflow Diagram Content

```text
GI paper + web/audit references
        |
        v
Table II parameters + reference spectra
        |
        v
radial grid + sqrt(p^2 + m^2) kinetic energy
        |
        v
closed-form G~(r), S~(r) + A(p) G~ A(p)
        |
        v
contact hyperfine + tensor + spin-orbit
        |
        v
predicted masses -> residuals -> scorecard
```

## Key Idea Blocks

### What Is Being Reproduced

- The target is the Godfrey-Isgur relativized quark model spectrum.
- The first completed validation target is heavy-heavy mesons:
  - charmonium, `c cbar`;
  - bottomonium, `b bbar`.
- Heavy-light and light sectors are tracked, but still need mixing and
  annihilation machinery before they should be judged as final.

### Why Layering Matters

Each physical correction has a different job:

| layer | job | visible symptom when missing |
|---|---|---|
| parameters/reference data | define the same problem as the paper | impossible to compare honestly |
| central Hamiltonian | set spin-averaged masses and radial/orbital spacings | global offsets and wrong level spacings |
| smearing/momentum sandwich | implement GI relativization of singular operators | heavy-heavy masses miss by tens to hundreds of MeV |
| contact hyperfine | split S-wave singlet and triplet states | `1^1S_0` / `1^3S_1` gap wrong |
| tensor/spin-orbit | split P-, D-, F-wave triplets | fine/hyperfine pattern distorted |
| diagnostics/tests | keep claims reproducible | no way to know which layer changed what |

### The 1.0 Active Path

Show this as four stacked toggle cards:

- `appendix_a_momentum_sandwich = true`
- `contact_momentum_sandwich = true`
- `fine_structure_momentum_sandwich = true`
- `fine_structure_smeared_kernels = true`

Caption: these are not arbitrary switches; they encode the physical path that
made heavy-quarkonium reproduction work.

## Code Anchors

- `data/parameters.provisional.toml`: active switches, Table II masses (`[masses]`), and paper parameters;
  load with `load_parameters_and_quark_masses` (see `docs/code_architecture.md`).
- `src/GIModel.jl`: module includes (types/setup/numerics), Hamiltonian, central operator,
  comparison reports.
- `src/spin_fine_structure.jl`: tensor and spin-orbit implementation.
- `docs/residual_reports/scorecard.md`: final sector-level performance.

## Infographic Panels

### Panel A: The Question

"Can we reproduce the GI paper numbers from formula-level ingredients, not by
blind tuning?"

Visual: paper mountain, llama carrying parameter packs, orca below the surface
labelled "operator details".

### Panel B: The Method

"Build one auditable layer at a time."

Visual: stacked transparent layers:

1. Data.
2. Central solver.
3. Smearing.
4. Momentum factors.
5. Spin shifts.
6. Reports/tests.

### Panel C: The Check

"Every layer leaves a diagnostic trail."

Visual: arrows from model layers to:

- tests;
- formula map;
- residual reports;
- scorecard.

### Panel D: The Outcome

Show the final 1.0 numbers:

| sector | mean abs residual | max abs residual |
|---|---:|---:|
| bottomonium | 4.0 MeV | 12.6 MeV |
| charmonium | 6.7 MeV | 18.0 MeV |

## What This Poster Should Make Clear

- The project is not one formula; it is a layered reproduction pipeline.
- The "stuck point" was not solved by retuning one number, but by expanding from
  the paper alone to web-supported GI-family formulas and then encoding those
  formulas in testable code.
- The final result is a heavy-quarkonium reproduction milestone, not yet a full
  light-meson model.

## What This Poster Should Not Claim

- Do not claim all meson sectors are solved.
- Do not claim the finite-difference workflow is identical to the original
  paper's harmonic-oscillator workflow.
- Do not hide the remaining unequal-mass and mixing gaps.

