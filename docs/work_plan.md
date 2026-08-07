# Reproduction Work Plan

> **Live status moved.** Per-unit status is the manifest/dashboard
> (`docs/paper_manifest/*.toml`; `cd docs && make dashboard`); workstream-level
> remaining-work costing is [`reproduction_audit.md`](reproduction_audit.md) §2;
> the task tracker is the live board; the *code*-stream plan is
> [`engineering_work_plan.md`](engineering_work_plan.md), with the focused
> paper-order architecture queue in
> [`paper_algorithm_work_plan.md`](paper_algorithm_work_plan.md). What stays
> here is the **process record**:
> the objective, how work was decomposed and delegated, the coordination rules,
> the empirical worker-sizing lesson, and the historical commit trail of what
> landed. Treat the "Open units" as illustrative history, not live truth.

Coordinator-maintained. Objective: close every gap in
[`docs/reproduction_audit.md`](reproduction_audit.md) so the project can claim
**"the original Godfrey-Isgur model is fully reproduced"**, then harvest the
demo/TIL material (direction 2).

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

- Files: `GIPaper/scripts/reproduce_table_v.jl`,
  `GIPaper/docs/residual_reports/table_v_reproduction.md` (regenerate). *(Done:
  the Q1/Q2 + Q1c/Q2c rows are scored via the singlet/triplet rotation in the
  unified reproduction harness.)*

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

### W6 — HO-order spin/mixing validation (retired into PA queue)

Paper-order perturbative comparison in the HO basis (A15-A17 fidelity). The
audit is complete; implementation follow-up is decomposed as PA-01 through
PA-18 in [`paper_algorithm_work_plan.md`](paper_algorithm_work_plan.md).

## Coordination rules

1. Workers commit in their worktree; the coordinator reviews, merges to
   `main`, runs both test suites and the report gate
   (`run_all_spectrum_checks.jl` must stay byte-identical except where a
   workstream legitimately extends its own report).
2. Ledger updates (`reproduction_audit.md`, `observable_ledger.md`, this
   file) are coordinator-only, done at merge time.
3. A workstream is **done** when its audit asset is committed, regenerable
   from a script, and the corresponding manifest unit is updated (status +
   metric + report) so the dashboard reflects it. Global/fitted-scale changes
   must be reflected in `reproduction_audit.md` §0 (the honesty ledger).

## Worker-sizing constraint (observed)

Delegated workers are terminated at a hard ~11-minute wall-clock ceiling
(session-limit). Empirically: code-shaped tasks that finish inside that
window survive and merge (W1, W2b); slow many-round-trip digitization /
page-image-verification tasks (W2a, W3) die having produced little. **Rule
going forward: a delegated unit must be completable in <~10 min.** The two
image-heavy streams are therefore split into single-deliverable units below,
launched one/two at a time (not swarmed) to cap shared-budget burn and blast
radius.

### W2a decomposition (per Table V section; each ~one page, ~30 rows)
- **W2a-1** light 1D remainder (1³D₁, 1³D₂, 1¹D₂) + light 2S — 1³D₃ already merged (`8823643`).
- **W2a-2** strange section.
- **W2a-3** charmed section (unblocks W2c).
- **W2a-4** charmonium (ψ) section.

### W3 decomposition (Table VI; single owned script/report, so serial not parallel)
- **W3-1** crop audit of the OCR-shifted open-flavor M1 column + 2S→χ₀ model-mass-q test (verification-only, no new rows — smallest, launch first).
- **W3-2** isoscalar mixing rows via existing Table III amplitudes.
- **W3-3** remaining light E1/M2 + strange/charmed P-wave + hindered bottomonium blocks.

## Status board

Cost/blocking detail for the open items lives in
[`reproduction_audit.md`](reproduction_audit.md) §2; this board tracks
ownership and merge state only.

### Completed and merged (13 units)

| Unit | Commit | One-line outcome |
| --- | --- | --- |
| W1 angles | `5340e62` | 13 angles audited; large-offdiag blocks match, small-offdiag 1P off (→ W6). |
| W2b K1 wiring | `c5d7ccd` | Q1/Q2 scored 3 ways; paper matches its own +34°. |
| W2a-1³D₃ | `8823643` | light 1³D₃ digitized, page-16 verified. |
| W2a-3 charmed | `44e1d50` | 10 charmed rows, page-22 verified; A_c/S_c/β_c convention documented. |
| W2a-1 light 2S/1D | `824eaed` | 94 rows incl. strange 2S/1D, page-209 verified (sign correction + 1³F₄ exclusion). |
| W3-1 crop+2S→χ₀ | `07ad617` | open-flavor M1 image-verified; 2S→χ₀ shown not a q-artifact. |
| W3-2 isoscalar | `4a8645c` | 6 isoscalar M1 rows; signs right, η↔η′ ordering open (→ W6). |
| W3-3 blocks | `d678136` | light E1/M2 + strange E1 + hindered b-b̄; **Table VI structurally complete**. |
| W2c-1 charm classes | `f2dadf2` | A_c/S_c coded no-refit; **closes Table IV gap**; Q1c/Q2c reproduce at −41°. |
| W2d score 2S/1D | `a8066b8` | 2S/1D/1³D₃ scored; structure-independent median 1%. SD (S/D/P) rows excluded: paper column not raw HO (realistic factor folded in). |
| W2a-4 1³F₄ | `0407272` | 1³F₄ nonet digitized+scored (no charmonium in Table V); SI median stays 1% S→F. |
| W5 Eqs. 20-21 | `dfb3007` | type-A realistic factors from hyperfine waves; R_A>1, monotonic in L, ~paper trend. |
| W2c-2 √3 diag | `fb5f5b4` | K*₂→Kπ √3 characterized: non-kinematic, missing Appendix-B isospin factor (open). |

(W2a-2 strange: **N/A** — strange 1S+1P already in the light CSV, strange 2S/1D folded into W2a-1.)

**Table V is now digitized and scored across every section.**

### Open units

The live physics open-work list is the task tracker and
[`reproduction_audit.md`](reproduction_audit.md) §2 (Table VII, the √3 recoil
factor, Eq. 21 type-S ratio, and the src/-promotion hygiene). W6's code work is
now the source-controlled PA queue in
[`paper_algorithm_work_plan.md`](paper_algorithm_work_plan.md).
