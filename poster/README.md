# Poster Packet: Godfrey-Isgur Reproduction 1.0

This folder contains six coherent markdown briefs for visual explanation posters.
They are written to become detailed infographic panels, not prose-only notes.
Each poster should stand alone on a wall, while also reading as one continuous
story from left to right.

The current visual theme is defined separately in `theme.md`. The scientific
content in the six poster files should not depend on llamas or orcas; those are
presentation mascots and visual metaphors only. To change style, replace
`theme.md` and keep the poster files intact.

## Narrative Spine

The 1.0 story is:

1. We start with the Godfrey-Isgur model and ask whether the published meson
   spectrum can be reproduced from the paper's ingredients.
2. We organize the reproduction as a transparent pipeline: paper parameters and
   reference rows enter a radial solver, then successive physical layers are
   added.
3. The central spin-independent Hamiltonian is the backbone. The decisive
   improvement is the closed-form smeared Appendix-A potential plus the Coulomb
   momentum sandwich.
4. Contact hyperfine splitting fixes S-wave singlet/triplet separation.
5. Tensor and spin-orbit fine structure clean up P-, D-, and F-wave multiplets.
6. The 1.0 scorecard shows heavy-heavy reproduction at few-MeV accuracy, and it
   clearly separates completed physics from known missing light-sector machinery.

## File Map

| file | poster | role |
|---|---|---|
| `01-introduction-global-workflow.md` | 1 | Big picture: what we are reproducing and the full pipeline. |
| `02-inputs-parameters-and-reference-data.md` | 2 | Layer 1: data inputs, Table II parameters, reference spectra, provenance checks. |
| `03-radial-solver-and-central-hamiltonian.md` | 3 | Layer 2: radial grid, kinetic energy, central potential, eigenstates. |
| `04-appendix-a-smearing-and-momentum-sandwiches.md` | 4 | Layer 3: GI smearing, closed-form `G~`, `S~`, and momentum-dependent central operator. |
| `05-spin-dependent-operators.md` | 5 | Layer 4: contact hyperfine, tensor, spin-orbit, smeared derivative kernels. |
| `06-results-and-remaining-physics.md` | 6 | Final result: scorecard, diagnostics, interpretation, next physics. |
| `theme.md` | all | Replaceable visual style guide. Current style: cartoon llamas and orcas. |

## Coherence Rules For All Posters

- Use the same color coding for physical layers:
  - Input/provenance: warm yellow.
  - Central Hamiltonian: deep teal.
  - Smearing/momentum relativization: blue.
  - Spin-dependent shifts: coral.
  - Validation/results: green.
  - Missing/future physics: gray-violet.
- Use the same three-line legend:
  - "Paper object" means a formula, parameter, or reference value from the GI
    paper or a GI-family reference.
  - "Code object" means a function, data file, or script in this repository.
  - "Diagnostic object" means a report, test, residual, or scorecard row.
- Put the current 1.0 active switches in a small badge on Posters 3-6:
  - `appendix_a_momentum_sandwich = true`
  - `contact_momentum_sandwich = true`
  - `fine_structure_momentum_sandwich = true`
  - `fine_structure_smeared_kernels = true`
- Show formulas as compact visual callouts. Do not turn the posters into a
  textbook page; each formula should be tied to one visual or one code object.
- Each poster should include a "What this layer explains" box and a "What it
  does not explain yet" box.

## Suggested Wall Layout

Use a 2 x 3 grid:

```text
[1 Workflow] [2 Inputs] [3 Central Solver]
[4 Smearing] [5 Spin]   [6 Results]
```

When arranged this way, the top row describes the base machine and the bottom row
shows the corrections that made version 1.0 work.

## Source Anchors

Use these repository anchors when building the visuals:

- Parameters: `data/parameters.provisional.toml`
- Julia package: **`GIModel`** (`Project.toml` at repo root; module entry `src/GIModel.jl`).
- Formula audit: `docs/formula_map.md`
- Research notes: `docs/research_spin_independent_gi.md`
- Core implementation: `src/GIModel.jl`
- Fine structure implementation: `src/spin_fine_structure.jl`
- Central path status: `src/appendix_a_status.jl`
- Tests: `test/runtests.jl`
- Scorecard: `docs/residual_reports/scorecard.md`
- Heavy-quarkonium diagnostics:
  `docs/residual_reports/heavy_quarkonium_diagnostics.md`

## Version 1.0 Claim

The version tagged `1.0` reproduces heavy quarkonium at:

| sector | rows | mean absolute residual | maximum absolute residual |
|---|---:|---:|---:|
| bottomonium | 30 | 4.0 MeV | 12.6 MeV |
| charmonium | 28 | 6.7 MeV | 18.0 MeV |

The posters should treat those numbers as the final arrival point of the story,
not as an isolated table.

