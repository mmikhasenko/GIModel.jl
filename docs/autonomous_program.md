# Autonomous Program: Godfrey-Isgur Reproduction

This file is the local equivalent of an `autoresearch` program file. It tells an
autonomous coding agent what progress means in this repository.

## Mission

Turn the current diagnostic Godfrey-Isgur solver into a paper-auditable,
reproducible implementation of the 1985 relativized quark model.

Do not chase residual improvements by tuning undocumented factors. Prefer small
auditable changes that make one discrepancy class clearer.

## Current Strategic Priorities

Work in this order unless the code state makes a later item an obvious
prerequisite:

1. Formula and convention audit.
   - Keep `docs/formula_map.md` synchronized with active code.
   - Check signs, factors, color factors, units, and radial implementation.
   - Add targeted regression tests for any formula that is fixed or clarified.

2. Radial normalization audit.
   - Decide whether finite-difference eigenvectors represent reduced radial
     functions or full radial wavefunctions.
   - Make contact, tensor, and spin-orbit expectation values use one documented
     normalization convention.
   - Remove or justify any suspicious `4π` factors.

3. Heavy-quarkonium validation.
   - Treat `ccbar` and `bbbar` as the main validation sectors.
   - Preserve diagnostics that separate common offsets, spin-averaged spacings,
     and spin splittings.
   - Explain residual classes before adding new model freedom.

4. Appendix A central potential.
   - Replace or clearly bracket the current pointwise central potential with a
     paper-faithful Appendix A implementation.
   - Keep before/after diagnostics so the effect is visible.

5. Heavy-light mixing.
   - Implement `^1L_J`/`^3L_J` and same-`J` tensor mixing as explicit mass-matrix
     diagonalization only after the central/fine-structure base is trustworthy.
   - Reports should show unmixed and mixed predictions side by side.

## Guardrails

- Do not silently refit parameters.
- Do not promote `k_spin_orbit`, `k_tensor`, or other bridge factors to physics
  parameters unless they are explicitly paper-derived.
- Do not expand claims in light or isoscalar sectors until missing mixing and
  annihilation machinery are implemented.
- Do not overwrite raw extraction/provenance data without preserving source,
  page/table/figure, extraction method, and confidence.
- Do not use later quoted GI tables as primary authority when the original
  paper can answer the question.
- If a discrepancy appears, classify it before changing code:
  extraction error, convention error, solver/numerics error, physics-model
  implementation error, or later-source mismatch.

## Standard Verification

Run the relevant subset during development, and expect the autonomous runner to
run the full gate before accepting an iteration:

```bash
python3 scripts/validate_seed.py
julia test/runtests.jl
julia scripts/analyze_heavy_quarkonium.jl
julia scripts/run_all_spectrum_checks.jl
```

## Preferred Work Unit

When work is *clear a priori* (same file family, one physics goal, one
verification), **batch it in a single session** and ship one meaningful commit
or a short chain (code + docs + test), not many micro-commits. Reserve small
patches for truly isolated fixes.

Each reviewable unit should still be **auditable**:

- one formula fix plus a regression test;
- one normalization clarification plus updated diagnostics;
- one report improvement that exposes a residual class;
- one documentation update that removes stale or contradictory project state;
- or a **coherent slice** of the completion path in `docs/orchestrator_task.md`
  (e.g. Table II lock + `verify_table_ii_toml` + doc update) in one go.

The final response for each iteration should say:

- what changed;
- why it advances the reproduction goal;
- which verification commands passed;
- which risk remains.

