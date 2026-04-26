# Midterm Review Brief

This document is for an external reviewer. Please be strict. Treat numerical
agreement, formula provenance, units, signs, normalization, and sector
conventions as possible failure points until they are demonstrated. Do not
reward progress that advances while blind to bugs, unclassified discrepancies,
or missing physics.

## Project Goal

Reproduce the Godfrey-Isgur relativized quark model meson spectra locally,
using the original 1985 paper as the primary authority and the digitized Figure
3-9 spectra as current reference targets.

The implementation is not complete. The current code is a diagnostic solver
that exposes which missing terms matter next.

## Current Working Commands

Run these from the repository root:

```bash
python3 scripts/validate_seed.py
julia --project=. test/runtests.jl
julia --project=. scripts/run_all_spectrum_checks.jl
julia --project=. scripts/analyze_heavy_quarkonium.jl
```

Expected generated/relevant reports:

- `docs/residual_reports/scorecard.md`
- `docs/residual_reports/heavy_quarkonium_diagnostics.md`
- `docs/residual_reports/*_residuals.md`
- `docs/residual_reports/ccbar_baseline.md`
- `docs/residual_reports/bbbar_baseline.md`

## Current Scorecard

From `docs/residual_reports/scorecard.md` after commit `fa52c02`:

| sector | rows | mean abs MeV | max abs MeV |
|---|---:|---:|---:|
| bottomonium | 30 | 38.8 | 77.2 |
| b_flavored | 21 | 58.7 | 303.5 |
| charmonium | 28 | 79.3 | 129.0 |
| charmed | 22 | 93.8 | 269.6 |
| strange | 30 | 205.7 | 503.2 |
| isoscalar | 49 | 306.2 | 1408.8 |
| isovector | 30 | 328.0 | 1157.4 |

Heavy-heavy and heavy-light sectors are currently useful diagnostics. Light and
isoscalar sectors are not complete because annihilation and mixing machinery is
still missing.

## Encouraging Result, Not Proof

`docs/residual_reports/heavy_quarkonium_diagnostics.md` decomposes the heavy
quarkonia residuals:

- charmonium common offset: `+79.3 MeV`
- charmonium mean absolute residual after removing offset: `15.6 MeV`
- bottomonium common offset: `+38.8 MeV`
- bottomonium mean absolute residual after removing offset: `10.0 MeV`

This suggests the central spectrum is coherent, but it does not prove the
Hamiltonian is correct. The offset and fine-structure defects still require
explanation from paper-derived terms, not tuning.

## Implemented Pieces

- Top-level working reference data in `data/reference_spectrum_*.csv`.
- Table II parameters in `data/parameters.provisional.toml`.
- Semirelativistic radial finite-difference kinetic energy:
  `sqrt(p^2 + m_1^2) + sqrt(p^2 + m_2^2)`.
- Pointwise central potential:
  `b r - 4 alpha_s(r)/(3r) + c`.
- Fig. 2 running-coupling ansatz.
- Smeared S-wave contact hyperfine diagnostic.
- First-order diagnostic spin-orbit and tensor terms.
- Unequal-mass constituent selection from `quark_content`.
- All-sector report generation.
- Heavy-quarkonium residual decomposition.

## Known Provisional Or Dangerous Parts

These are not final GI reproduction pieces:

- The central potential is not the full Appendix A GI effective potential.
- `appendix_a_smearing` exists but is off by default. The current naive 3D
  convolution is explicitly experimental and should not be promoted.
- Tensor radial term currently uses an unsmeared `alpha_s(r)/r^3` diagnostic
  proxy for `L>0`. This was chosen because using the broad contact width
  overdamped splittings. It is not the final Appendix A tensor operator.
- Fine-structure scale factors `k_spin_orbit` and `k_tensor` are diagnostic
  bridge factors. They are not paper-derived parameters and should be treated
  with suspicion.
- No proper off-diagonal tensor or antisymmetric spin-orbit mixing is
  implemented.
- Isoscalar annihilation and `n nbar`/`s sbar` mixing are missing.
- Light-sector residuals are large and should not be interpreted as a completed
  reproduction.
- The reference spectra are digitized from figure labels, not clean promoted
  `data/clean` products.

## Review Priorities

Please review in this order.

1. Check formulas against the paper.
   - `src/GIModel/GIModel.jl`
   - `src/GIModel/spin_fine_structure.jl`
   - `docs/formula_map.md`
   Verify signs, factors of 2/3/4, color factors, units, and whether the radial
   implementation corresponds to the formula claimed.

2. Check normalization.
   The finite-difference eigenvectors are discrete radial vectors. Verify all
   expectation values use a consistent normalization. The current code uses a
   `4π` normalization helper in spin diagnostics; this may be conceptually
   inconsistent with reduced radial functions. Be intolerant here.

3. Check the heavy-light mass parsing.
   `src/GIModel/masses_from_content.jl` previously had a bug class where
   `composition_raw` was mistaken for `quark_content`. Confirm this is truly
   fixed and covered.

4. Check whether fine-structure “improvements” are real.
   The latest changes improved several reports but worsened high-L bottom-light
   states. Decide whether this is acceptable evidence for missing mixing or a
   sign/radial-integral bug.

5. Check diagnostics for blindness.
   The project should not tune parameters or add terms without making the
   residual class clearer. Any undocumented factor or fudge should be rejected.

6. Check generated reports for reproducibility.
   Delete reports, rerun commands, and confirm they regenerate identically.

## Specific Questions For Reviewer

- Is the spin-orbit radial coefficient currently coded with the right sign and
  relative vector/scalar factor?
- Are the triplet tensor angular factors correct for the operator convention
  used by Godfrey-Isgur?
- Is the finite-difference radial basis compatible with the expectation-value
  formulas used for contact, tensor, and spin-orbit terms?
- Should the contact hyperfine expectation contain the current normalization
  and smearing factors, or is it dimensionally inconsistent?
- Are the `epsilon_*` factors being applied in a way that is defensible, or are
  they currently cosmetic?
- Is the common offset in heavy quarkonia an expected missing Appendix A effect,
  a bad constant convention, or an implementation bug?
- Should the next implementation step be Appendix A central potential or
  heavy-light mixing?

## Suggested Next Step If Review Passes

Implement heavy-light `^1L_J`/`^3L_J` mixing and same-`J` tensor mixing as
explicit mass-matrix diagonalization, with reports showing unmixed and mixed
predictions side by side.

## Suggested Next Step If Review Finds Formula Bugs

Stop sector expansion. Fix the formula, add a targeted regression test, and
regenerate:

```bash
julia --project=. test/runtests.jl
julia --project=. scripts/analyze_heavy_quarkonium.jl
julia --project=. scripts/run_all_spectrum_checks.jl
```

Do not tune around a formula bug.

## Current Git Context

Recent relevant commits:

- `fa52c02` Strengthen diagnostic spin splittings
- `a3777e5` Add heavy quarkonium mismatch diagnostics
- `7c1d2e4` Add spectrum scorecard report
- `22bd085` Correct fine-structure angular factors
- `017bda3` Extend spectrum checks across sectors

The reviewer should inspect the diff from before `017bda3` if they want to
audit the full solver evolution.
