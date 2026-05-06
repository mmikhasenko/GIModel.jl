# Isoscalar Annihilation Progress Metric

This metric is the gate for step 4: literal Godfrey-Isgur isoscalar
annihilation and pseudoscalar P1/P2 implementation. It is designed to make
progress measurable without hiding extraction or convention uncertainty.

## Score

Use a 100-point score, with each component computed from committed artifacts.

| Component | Points | Gauge |
| --- | ---: | --- |
| Formula provenance | 20 | Eq. (16), Eq. (17), Eq. (18a), Eq. (18b), and Table III each have a PDF/image-audited row in `docs/formula_map.md` or a linked audit note. Four points per item. |
| Clean targets | 20 | `data/clean/masses.csv` has the relevant isoscalar pseudoscalar masses, and `data/clean/mixings.csv` has audited Table III P1/P2 rows with provenance. Ten points for masses, ten for mixings. |
| Implementation modes | 20 | Dispatch exposes `:none`, `:calibrated_p1`, literal `:paper_p1`, and literal `:paper_p2`; each mode has at least one regression test and does not change Table II parameters. Five points per mode. |
| Spectral fidelity | 20 | For the audited `1^1S_0` and `2^1S_0` isoscalar target masses, compute mean absolute residual. Award 20 points at `<= 25 MeV`, 10 points at `<= 50 MeV`, otherwise 0. P1 and P2 should be scored separately. |
| Mixing fidelity | 20 | Compare normalized eigenvectors to audited Table III amplitudes up to overall sign. Award 20 points if RMS amplitude error is `<= 0.10`, 10 points if `<= 0.20`, otherwise 0. P1 and P2 should be scored separately. |

## Required Report

Step 4 should generate `docs/residual_reports/annihilation_model_scorecard.md`
with one row per model:

```text
model, formula_points, clean_target_points, implementation_points,
spectral_points, mixing_points, total_points, mean_abs_mass_residual_MeV,
mixing_rms_error
```

The existing `:calibrated_p1` diagnostic is allowed to score implementation and
spectral points, but it must score zero formula-provenance points for literal
P1/P2 because it is a reconstruction control rather than Eq. (18a).

## Progress Gates

- **Ready to code literal P1/P2:** at least 30/40 combined formula-provenance
  plus clean-target points.
- **Ready to compare physics:** at least 70/100 for either `:paper_p1` or
  `:paper_p2`.
- **Ready to retire calibrated diagnostic as the main report path:** a literal
  paper model scores at least 80/100 and does not worsen the non-isoscalar
  scorecards, since annihilation must remain opt-in to self-conjugate isoscalar
  sectors.
