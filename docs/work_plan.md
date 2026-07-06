# Reproduction Work Plan

Coordinator-maintained. Objective: close every gap in
[`docs/reproduction_audit.md`](reproduction_audit.md) so the project can claim
**"the original Godfrey-Isgur model is fully reproduced"**, then harvest the
demo/TIL material (direction 2). This file is the single status board; workers
do **not** edit it — they deliver into their owned file scopes and the
coordinator integrates, updates statuses here, and keeps the audit document in
sync.

## Demand assessment

| Rank | Workstream | Why demanding |
| --- | --- | --- |
| 1 | **W2 Table V completion** | Multipage digitization with crop-level OCR verification; new physics (charm classes `A_c`/`S_c`, unequal-mass recoil behind the strange √3 anomaly); cross-layer wiring (spectrum-layer K1 angle into decay rows). Split into W2a/W2b/W2c so data work, wiring, and new physics proceed independently. |
| 2 | **W4 Table VII** | New formula layer (D4-D9: decay constants, leptonic/γγ/gluonic widths, charge radii) plus digitization of a three-page table; everything runs on existing wavefunctions, so it is broad rather than deep. |
| 3 | **W3 Table VI completion** | Mostly assembly (mixing amplitudes exist; kernels exist) plus a careful PDF crop audit of the OCR-shifted column. |
| 4 | **W6 HO-order spin/mixing validation** | Research-flavored (paper-order perturbation theory in the HO basis); background item, not delegated until W1-W4 land. |
| 5 | **W1 Mixing angles**, **W5 Eqs. (20)-(21)** | Small; all machinery exists. |

## Workstreams

Each workstream owns a disjoint file scope (no worker touches `src/` solver
internals, shared ledgers, or another stream's files; new shared-code needs go
through the coordinator). Workers run in isolated git worktrees and commit
there; the coordinator merges.

### W1 — Mixing-angle audit (small, wave 1)

Computed-vs-paper table for the 13 quoted angles (strange θ₁P 34°, θ₁D 33°,
θ₂P 15°, θ₁F 32°, θ₂D 25°, θ₁G 33°; charm cū 1P −41°/1D −39°, cs̄ 1P −44°/1D
−39°; bottom bū −43°, bs̄ −45°, bc̄ −53°), with the paper's rotation
convention pinned once against our `StateMixing` convention.

- Files: `GIPaper/scripts/audit_mixing_angles.jl`,
  `GIPaper/docs/residual_reports/mixing_angles.md`.
- Accept: all 13 angles computed from `compute_spectrum`; convention mapping
  derived, not tuned; deviations explained or flagged.

### W2a — Table V digitization (large, wave 1)

Digitize the remaining Table V sections from the vision OCR (lines 593-1048 of
`paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`) with page-image
verification: light 2S and 1D, strange, charmed, charmonium sections.

- Files: `data/raw/digitized_tables/table_v_strong_decays/*.csv` (new files
  following the existing `table_v_light_1s_1p.provisional.csv` schema),
  `data/raw/digitized_tables/table_v_strong_decays/DIGITIZATION_NOTES.md`.
- Accept: every row carries coefficient expression, class, paper amplitude,
  confidence, and page/crop provenance; low-confidence cells flagged.

### W2b — K1 angle into Q1/Q2 decay rows (medium, wave 2, after W1)

Wire the model's strange 1P same-J mixing angle from the spectrum layer into
the Table V audit so the 15 excluded Q1/Q2 rows become scoreable.

- Files: `GIPaper/scripts/audit_table_v_decays.jl` (extend),
  `GIPaper/docs/residual_reports/table_v_light_decays.md` (regenerate).

### W2c — Charm classes and sections (medium-large, wave 2, after W2a)

`A_c`, `S_c`, `beta_c` classes in `src/strong_decays.jl` (+tests), audit of
the charmed/charmonium Table V sections; investigate the strange-parent √3
normalization (Appendix-B unequal-mass recoil factors).

### W3 — Table VI completion (medium, wave 1)

Isoscalar mixing rows via the existing Table III mixing amplitudes; PDF crop
audit of the OCR-shifted open-flavor column; remaining E1/M2 and
strange/charmed blocks; 2S → χ₀ momentum-convention check.

- Files: `GIPaper/scripts/audit_table_vi_photon_decays.jl` (extend),
  `GIPaper/docs/residual_reports/table_vi_photon_decays.md` (regenerate).

### W4 — Table VII (medium, wave 2)

D4-D9 formulas (f_P, f_V, f_{A₁}; leptonic, γγ, gluonic widths), charge radii
from model wavefunctions, digitization of Table VII (OCR lines 1275-1410+),
audit report.

- Files: `GIPaper/scripts/audit_table_vii.jl`,
  `data/raw/digitized_tables/table_vii/*.csv`,
  `GIPaper/docs/residual_reports/table_vii_annihilation_em.md`.

### W5 — Eqs. (20)-(21) realistic factors (small, wave 2, after W2a merges)

Radial-moment ratios on model wavefunctions; adds the parenthetical correction
column to the Table V audit.

### W6 — HO-order spin/mixing validation (background, unscheduled)

Paper-order perturbative comparison in the HO basis (A15-A17 fidelity).
Coordinator-owned until W1-W4 land.

## Coordination rules

1. Workers commit in their worktree; the coordinator reviews, merges to
   `main`, runs both test suites and the report gate
   (`run_all_spectrum_checks.jl` must stay byte-identical except where a
   workstream legitimately extends its own report).
2. Ledger updates (`reproduction_audit.md`, `observable_ledger.md`, this
   file) are coordinator-only, done at merge time.
3. A workstream is **done** when its audit asset is committed, regenerable
   from a script, and its acceptance criterion in
   `reproduction_audit.md` §3 is met.

## Status board

| Stream | Status | Notes |
| --- | --- | --- |
| W1 angles | **done** (merged `5340e62`) | `mixing_angles.md`: 1D blocks 0.3-2.5°, b-sector 1-5° after low/high label exchange (near-degenerate eigenstates); small-offdiagonal 1P blocks 13-30° off with inverted radial trend — physics finding, not a bug. Worker interrupted by session limit; coordinator completed. |
| W2a Table V digitization | relaunch pending | wave-1 worker cut off by session limit before output; relaunch after limit reset |
| W2b K1 wiring | queued | W1 done — unblocked; note the K1 (u sbar 1P) angle itself deviates (+19.5° vs paper +34°), so score Q1/Q2 rows with both angles |
| W2c charm classes | queued | blocked by W2a (charm rows) |
| W3 Table VI | relaunch pending | wave-1 worker cut off by session limit before output; relaunch after limit reset |
| W4 Table VII | queued | wave 2 |
| W5 Eqs. 20-21 | queued | blocked by W2a merge |
| W6 HO-order | unscheduled | coordinator; the small-offdiagonal 1P angle sensitivity found by W1 strengthens the case for this item |
