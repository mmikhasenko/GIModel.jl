# GIPaper investigations

This directory contains reproducible, hypothesis-driven studies of mechanisms
inside the Godfrey–Isgur calculation.

It is deliberately separate from `../residual_reports/`:

- `residual_reports/` is the authoritative paper-reproduction ledger. It asks
  which published values differ and records every comparison row.
- `investigations/` asks why a localized discrepancy occurs. A study may
  separate operators, compare numerical realizations, or test a physical
  hypothesis, but it must link back to the canonical residual report rather
  than copy its paper data or implement a second comparison path.

## Standards for an investigation

Every numerical study here should:

1. name its executable generator and numerical settings;
2. use GIModel's public computation API wherever possible;
3. distinguish coordinate-kernel plots from complete momentum-sandwiched
   matrix elements;
4. include an independent numerical control when the conclusion depends on a
   solver;
5. state what was ruled out, what remains unresolved, and the next calculation
   that could discriminate between explanations;
6. keep generated figures under `figures/` and avoid embedding paper reference
   values already owned by a residual report.

## Studies

- [Light and charmed P-wave fine structure](pwave_fine_structure.md) — separates
  vector spin-orbit, scalar/Thomas, and tensor kernels; compares central-wave
  first order with full fixed-sector diagonalization; and checks native HO
  against finite differences.

Regenerate the P-wave study from the repository root with:

```sh
julia GIPaper/scripts/investigate_pwave_fine_structure.jl
```
