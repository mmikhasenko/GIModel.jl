# Poster 6: Results And Remaining Physics

## Purpose

Show what version 1.0 achieves, why the result is meaningful, and what remains
outside the current implementation. This should be the satisfying endpoint of
the six-poster sequence without overselling the scope.

## Headline

Version 1.0: heavy quarkonium reproduced at few-MeV accuracy.

## One-Sentence Takeaway

With Appendix-A central smearing, momentum sandwiches, contact hyperfine, and
smeared-derivative fine structure active, the model reproduces bottomonium and
charmonium well enough that remaining discrepancies are now detailed convention
and missing-physics questions rather than gross solver failures.

## Main Visual

Create a large scorecard centerpiece:

| sector | rows | mean abs residual | max abs residual | interpretation |
|---|---:|---:|---:|---|
| bottomonium | 30 | 4.0 MeV | 12.6 MeV | reproduced |
| charmonium | 28 | 6.7 MeV | 18.0 MeV | reproduced |
| b-flavored | 21 | 17.0 MeV | 78.0 MeV | useful diagnostic |
| charmed | 22 | 33.0 MeV | 79.5 MeV | useful diagnostic |
| strange | 30 | 60.8 MeV | 258.9 MeV | incomplete physics |
| isovector | 30 | 96.4 MeV | 464.9 MeV | incomplete physics |
| isoscalar | 49 | 191.9 MeV | 1191.0 MeV | missing annihilation/mixing |

Use green check styling only for bottomonium and charmonium. Use neutral
diagnostic styling for heavy-light. Use gray-violet "future physics" styling for
light and isoscalar sectors.

## Result Sources

- `docs/residual_reports/scorecard.md`
- `docs/residual_reports/charmonium_residuals.md`
- `docs/residual_reports/bottomonium_residuals.md`
- `docs/residual_reports/ccbar_baseline.md`
- `docs/residual_reports/bbbar_baseline.md`
- `docs/residual_reports/heavy_quarkonium_diagnostics.md`

## Heavy-Quarkonium Diagnostics

Show the deeper interpretation:

### Charmonium

- rows: 28;
- mean residual/common offset: about `-0.6 MeV`;
- mean absolute residual: `6.7 MeV`;
- radial/orbital spacing mean absolute error: about `2.2 MeV`.

Interpretation:

- multiplet centers and radial/orbital spacings are now excellent;
- remaining visible tension is mostly fine splitting, especially P-wave triplet
  patterns and some S-wave hyperfine details;
- this points toward extraction/rounding audit and paper workflow ordering, not
  a broken central Hamiltonian.

### Bottomonium

- rows: 30;
- mean residual/common offset: about `-0.8 MeV`;
- mean absolute residual: `4.0 MeV`;
- maximum absolute residual: `12.6 MeV`.

Interpretation:

- the heavy, compact system is reproduced extremely well;
- spin-dependent corrections are small and controlled;
- this is the cleanest validation of the 1.0 physical path.

## Before/After Story

Create a visual staircase:

```text
pointwise central only
    -> closed-form smeared central
    -> central momentum sandwich
    -> contact momentum sandwich
    -> fine-structure momentum sandwich
    -> smeared G~/S~ derivative kernels
    -> 1.0 heavy-quarkonium scorecard
```

Caption: "The result came from putting the missing operator layers in place, not
from hiding residuals with a fit."

## What Is Solved At 1.0

Use green cards:

- Heavy-heavy central masses and radial/orbital spacings.
- S-wave contact hyperfine mechanism.
- Equal-mass tensor and spin-orbit fine-structure mechanism.
- Closed-form smeared `G~`, `S~` central path.
- Momentum-dependent Hermitian sandwiches for central, contact, and
  fine-structure operators.
- Tests and reports that regenerate the numbers.

## What Is Not Solved Yet

Use gray-violet cards:

- Unequal-mass antisymmetric spin-orbit mixing.
- Same-`J` tensor mixing between basis states.
- Isoscalar annihilation and `n nbar` / `s sbar` flavor mixing.
- Coupled-channel effects.
- Exact comparison to the original paper's harmonic-oscillator basis ordering.
- Final extraction audit for every plotted/table reference value.

## Visual Interpretation Guide

Create a "read the residuals correctly" panel:

| sector type | how to read residuals |
|---|---|
| heavy-heavy | meaningful 1.0 reproduction target |
| heavy-light | useful diagnostic, missing unequal-mass mixing |
| light non-isoscalar | qualitatively informative, incomplete relativistic/mixing physics |
| isoscalar | not final until annihilation and flavor mixing are implemented |

## Validation Commands

Show as a reproducibility strip:

```bash
julia test/runtests.jl
julia scripts/analyze_heavy_quarkonium.jl
julia scripts/run_baseline_solver.jl
julia scripts/run_all_spectrum_checks.jl
julia scripts/compare_central_paths.jl
```

Caption: "The posters explain a pipeline that can be rerun."

## Final Message

Use a big concluding callout:

"1.0 is the point where the GI heavy-heavy reproduction becomes a controlled
physics problem instead of a debugging problem."

The final visual can show the llama and orca together at the scorecard:

- llama holds the reproducibility checklist;
- orca surfaces with the `G~`, `S~`, and spin-kernel badges;
- the scorecard sits in the center, clean and readable.

## Designer Notes

- Use validation green for the heavy-heavy result cards.
- Do not make all sectors green. The distinction between reproduced and
  diagnostic is central to scientific honesty.
- Put the exact 1.0 numbers in large readable type.
- Keep the future-physics list prominent enough that the poster feels accurate,
  not triumphalist.

