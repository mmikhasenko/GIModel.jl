# Paper Gap Ledger

> **Detail layer.** The at-a-glance "what's done / what's missing" is the
> manifest/dashboard (`docs/paper_manifest/*.toml`; `cd docs && make dashboard`).
> This file is the **technical rationale** behind the `partial`/`missing`
> statuses — the specific physics of each remaining gap (paper-order staging,
> full Appendix-A relativization, off-diagonal tensor mixing, unequal-mass
> spin-orbit, literal isoscalar annihilation, eigenvector-fidelity scoring), which
> the per-unit `notes` only summarize.

This is the current short list of what is still missing relative to the
Godfrey-Isgur paper. It is meant to be the durable replacement for older
handoff/autonomous planning notes.

## Implemented Enough For Current Comparisons

- Table II solver inputs are audited into CSV/TOML, including
  `epsilon_so_scalar = +0.055`.
- The semirelativistic kinetic operator is implemented in the finite-difference
  path.
- The oscillator path is implemented as `OscillatorSolver` and has regression
  tests. `scripts/audit_nonmixing_contact.jl`
  compares FD and HO non-mixing scorecards, and
  `scripts/audit_appendix_a_ho_comparison.jl` compares the active Appendix-A
  central operator directly.
- The central potential uses the GI color-singlet sign and normalization.
- Fig. 2 running `alpha_s(r)` is implemented and regression-tested.
- S-wave contact hyperfine is smeared and diagonalized nonperturbatively in the
  finite-difference path.
- Equal-mass spin-orbit and tensor diagonal shifts are implemented as active
  sector diagnostics.
- The generic mixing layer exists: `MixingMechanism`, `MixingBlock`,
  `diagonalize_mixing_block`, and comparison-layer assignment methods for
  `AntisymmetricSpinOrbit`, `TensorMixing`, and the calibrated
  `IsoscalarAnnihilation` control. `IsoscalarAnnihilation` also has paper
  pseudoscalar P1/P2 modes driven by the `[annihilation]` TOML constants, the
  general non-pseudoscalar Eq. (16) block
  (`isoscalar_general_annihilation_solution`), and the `:table_iii` scheme
  that combines calibrated P1, general `^3S_1`/`^3P_2` blocks, and ideal
  mixing for all other isoscalar channels.
  `compare` now assigns open-flavor `^1L_J`/`^3L_J` same-`J` pairs and triplet
  tensor `L/L'` pairs when partner rows are present.
- Top-level spectrum CSVs cover Figs. 3-9, and promoted clean mass/mixing data
  now exists under `data/clean/`.

## Missing Or Still Approximate

### 1. Eq. (14) Paper-Order Staging, Not "Missing HO"

What the paper does: diagonalizes the relativized Hamiltonian in a large
harmonic-oscillator basis, first in fixed `|jm;ls>` sectors, then diagonalizes
smaller mass matrices for tensor, antisymmetric spin-orbit, and annihilation
mixing.

What the repo does: both FD and HO radial paths exist. The HO path is selected
by passing an `OscillatorSolver`, builds the Hamiltonian in a finite oscillator
subspace, scans the oscillator scale, reconstructs mesh wavefunctions, and is
tested. Current headline sector reports still use
FD; the HO path is used mainly for basis-comparison audits.

Why it matters: the claim that "HO is missing" is wrong. The central
Appendix-A FD/HO comparison is now complete at the sub-MeV level. The
open-flavor antisymmetric spin-orbit and triplet tensor blocks are wired into
the FD comparison path, and the calibrated isoscalar pseudoscalar annihilation
control plus literal P1/P2 formula modes are routed through the same mixing
layer. The remaining gap is paper-order staging after the fixed-sector solve:
broader Table III flavor/radial mixing and HO-order validation of the
spin/mixing blocks.

Acceptance check: keep the existing FD-vs-HO scorecards, then add a paper-order
comparison mode that starts from HO fixed-sector eigenvectors and applies the
post-diagonalization mixing blocks before assigning physical rows.

### 2. Full Appendix-A Relativization

What the paper does: replaces naive mass factors in spin-dependent and smeared
operators by operator-level `m/E` factors with fitted epsilon exponents and
momentum-space sandwiching.

What the repo does: more is implemented than the old shorthand implied.
Contact and fine-structure momentum-factor sandwiches are active when the
corresponding parameter switches are enabled; the side exponent is
`1/2 + epsilon_i`, so the Table II epsilon values are used as exponents in the
GI-style Hermitian sandwich. The central path also has a momentum sandwich.

Why it matters: this is now a follow-up/audit gap, not a blank missing feature.
The central A12-A14 operator has been compared between FD and HO bases. The
remaining uncertainty is whether every spin-dependent kernel and operator
ordering matches the exact Appendix-A HO matrix-element prescription and the
paper's perturbative ordering.

Acceptance check: add a report comparing FD-sandwich and HO-sandwich matrix
elements for contact, tensor, vector spin-orbit, and scalar spin-orbit terms on
`ccbar`, `ssbar`, and light `q qbar`, then wire the resulting paper-order
blocks into sector assignment.

### 3. Off-Diagonal Tensor Mixing

What the paper does: the tensor part of the hyperfine interaction mixes
`^3L_J` and `^3L'_J` states with the same `J` after the first fixed-sector
diagonalization.

What the repo does: the generic mixing level is present and tested. Tensor
diagonal shifts, the same-`J` off-diagonal angular factor, cross-radial
matrix elements, and `assign_mixed_rows(::TensorMixing, ...)` are implemented.
Partnered triplet `L=J-1`/`L=J+1` rows are assigned in the comparison reports.

Why it matters: this affects states such as triplet S/D candidates and can move
assignments, not just masses. Partnered tensor rows are now judged in the main
scorecards instead of being excluded as unresolved.

Acceptance check: complete for available partnered rows. Residual reports now
include same-`J` tensor mixing tables and the non-mixing scorecard includes
assigned tensor pairs.

### 4. Unequal-Mass Antisymmetric Spin-Orbit Mixing

What the paper does: for unequal constituent masses, the antisymmetric
spin-orbit term mixes `^1L_J` and `^3L_J` states.

What the repo does: the radial antisymmetric matrix element, angular convention,
`same_j_mixing` diagonalizer, equal-mass/off-equal-mass tests, and comparison
assignment layer now exist. `compare` detects candidate open-flavor
`^1L_J`/`^3L_J` pairs, builds the `2x2` mass matrix from the cached radial
solution and diagonal fine-structure masses, assigns the mixed eigenvalues by
reference-mass order, and records the off-diagonal element, mixing angle, and
singlet/triplet components.

Why it matters: heavy-light P-, D-, F-, and G-wave same-`J` rows are now
scoreable in the main residual reports instead of being excluded as unresolved
assignment/mixing failures.

Acceptance check: complete. Sector reports for strange, charmed, and
bottom-flavored states include a same-`J` spin-orbit mixing table, and assigned
rows are no longer excluded from the non-mixing scorecard.

### 5. Literal Isoscalar Annihilation, Eq. (16)-Eq. (18)

What the paper does: for self-conjugate isoscalars, it adds an annihilation
mass matrix after the fixed-sector diagonalization. Eq. (16)-Eq. (17) define
the general matrix element; Eq. (18a) and Eq. (18b) replace the pseudoscalar
bracket for P1 and P2.

What the repo does: the general Eq. (16) machinery now exists
(`isoscalar_general_annihilation_solution`) with the `4π(2L+1)` prefactor, the
`(α_s(M_j²)α_s(M_i²)/π²)^{n/2}` bracket with `n = 2`/`3` for `C = +`/`−`, and
the Eq. (17) `S_L` smearing evaluated by a `j_L` momentum transform on
HO-basis radial wavefunctions. The earlier "scale problem" had two causes,
both fixed: the FD wavefunction-at-origin scale (resolved by the HO wave
cache) and a wrong two-gluon bracket in the `^3S_1` channel (resolved by the
`n = 3` power for `C = −`). With Table II inputs and `A(^3S_1) = +2.5` the
omega/phi block reproduces Table III: amplitudes `(+1.000, -0.029)` vs
`(+0.999, -0.02)` and an omega shift of `+13 MeV` vs the paper's `+10 MeV`.
`A(^3P_2) = -0.8` reproduces the f2/f2' row, `(+0.997, +0.080)` vs
`(+0.997, +0.06)`. The literal pseudoscalar `:paper_p1`/`:paper_p2` modes on
HO waves reproduce the Table III sign structure in every row once the GI
annihilation phase convention `Φ(0) > 0` is enforced
(`fix_annihilation_phase!`): the radially excited `2 ns`/`2 ss` couplings are
then negative, which is what makes the published all-positive eta-prime row
possible. Mean amplitude RMS is 0.078 (P1) and 0.025 (P2) — the P2 eta_r row
`(+0.993, +0.076, +0.069)` vs the paper's `(+0.99, +0.07, +0.08)` is
essentially exact. Literal-mode masses sit 47 MeV (P2) / 76 MeV (P1) mean
from the paper-model targets, inherited mostly from the light unperturbed
diagonals (the FD pi is ~55 MeV below the paper's 0.15 GeV).

Why it matters: this was the last hundreds-of-MeV gap in the isoscalar
spectrum reproduction.

Acceptance check: complete for the non-pseudoscalar channels (eigenvectors and
splittings) and for the literal pseudoscalar eigenvectors. The calibrated P1
control still carries the headline mass scoring; retiring it requires fixing
the shared light-sector diagonal residuals, not the annihilation block.

### 6. General Table III Isoscalar Mixing

What the paper does: Table III gives approximate isoscalar compositions beyond
the pseudoscalar rows, including non-pseudoscalar sectors with predicted and
observed splittings, and assumes ideal mixing for every channel it does not
list.

What the repo does: the `:table_iii` comparison scheme implements exactly that
prescription: calibrated P1 for `^1S_0`, general Eq. (16) blocks for `^3S_1`
and `^3P_2`, and ideal `n nbar`/`s sbar` mixing for every other isoscalar
channel, with full contact + fine-structure diagonals for the `s sbar`
partner rows. The isoscalar residual report now scores all paired rows: the
sector mean dropped from 97.3 MeV to 13.5 MeV and the non-mixing scorecard
includes 40 isoscalar rows at 13.6 MeV mean — in line with the isovector
sector, whose channel residuals the isoscalar rows now track.

Acceptance check: complete. The previously unpaired `2^3D_2` digitized row at
2.26 GeV was crop-audited (400 dpi, page 8) and rejected: the Fig. 5 `2--`
column contains only `1^3D_2(1.70)` and `1^3D_2(1.91)`; the 2.26 value was
contamination from the Fig. 4 strange `2^3D_2(2.26)` label and has been
removed from the reference spectrum.

### 7. Eigenvector Fidelity Scoring

What the paper gives: Table III amplitudes are eigenvector/composition targets,
not just mass targets.

What the repo does: `docs/residual_reports/table_iii_mixing_audit.md` computes
the amplitude RMS for literal P1/P2 on HO waves against the promoted Table III
rows plus the literal-mode mass residuals against the paper-model targets
(Fig. 5 labels for P1; Table III Δm + Sec. VA pole masses for P2), and
`data_checks.py score-annihilation` folds both into the annihilation
scorecard as mixing and spectral points with shared thresholds.

Acceptance check: complete as a scoring mechanism, including literal-mode
spectral scoring. Remaining quality work is the shared light-sector
unperturbed diagonal that caps the literal-mode mass points.

### 8. Observables Outside The Mass Spectrum

What the paper includes: strong decay amplitudes (Sec. IV, Tables IV-V,
Appendices B-C), leptonic, two-photon, gluonic decay amplitudes, charge
radii, and decay-model discussion beyond the mass spectrum.

What the repo does: `docs/observable_ledger.md` (the acceptance artifact)
tracks decay conventions separately from the mass code. The Table IV/V
two-parameter strong-decay model is implemented in `src/strong_decays.jl`
(row-oriented `decay_amplitude` API; `A = 1.665` from `rho -> pi pi`,
`S0 = 3.287` from `B -> [omega pi]_S` in the leading-S0 convention,
`beta = 0.40 GeV`), all of Table V is digitized in
`data/raw/digitized_tables/table_v_strong_decays.csv`, and
`scripts/reproduce_table_v.jl` reproduces the paper's numeric amplitude
column on **160 / 178 scoreable rows**, with the remaining non-matches shown
(via implied-mass inversion + the mixing rotation) to be parent-mass /
mixing-angle input sensitivity, and the deferred conventions (quasi-two-body
lineshapes, sub-threshold modes) itemized in the observable ledger. Table VI
photon decays are now started: `scripts/audit_table_vi_photon_decays.jl`
evaluates the Appendix-D mock-meson overlaps `I_i(x,y)` / `E_n^i(x,y)` on
the model's own FD wavefunctions and reproduces 28 mixing-free M1/E1 rows
with no new fitted constants — quarkonium and open-flavor M1 at the 0.1-4%
level, E1 chi triplets at 3-7% (see the observable ledger for the two
open q-convention rows and the OCR column-shift caveat).

Why it matters: these observables test wavefunctions, not only eigenvalues.
The Table VI agreement is direct evidence that the FD wavefunctions — not
just the eigenvalues — match the paper's, including the contact-distorted
S waves (the hindered `psi' -> eta_c gamma` row only exists through the
singlet/triplet wavefunction difference).

Acceptance check: observable ledger exists; the first Table V block and the
first 28 Table VI rows are scored. Remaining: wire the model K1 mixing angle
into the decay audit, resolve the strange-parent recoil normalization via
Appendix B, extract the later Table V and Table VI sections, fold Table III
mixings into the isoscalar Table VI rows, and start Table VII.

## Next Clean Implementation Step

Items 3-7 are now implemented (general Eq. (16), `:table_iii` ideal-mixing
prescription, the `Φ(0) > 0` annihilation phase convention, eigenvector and
literal-mode spectral scoring, and the `2^3D_2` crop audit). The most direct
remaining order is:

1. Reduce the shared light-sector unperturbed residuals that now cap both the
   light scorecards and the literal pseudoscalar masses: the FD pi sits
   ~55 MeV below the paper's 0.15 GeV, and the P-wave singlets (`1^1P_1`,
   `1^3P_0`, `2^3P_0`) carry +25..40 MeV in both isovector and isoscalar.
2. Revisit Appendix-A/HO paper-order staging with the now-explicit mixing
   blocks (items 1-2 above): start from HO fixed-sector eigenvectors and apply
   the post-diagonalization mixing blocks before assigning physical rows.
3. Once 1-2 settle, retire the calibrated P1 control from the headline
   isoscalar report in favor of the literal P1 mode.
4. Observables layer (item 8, started): promote the Table VI overlap kernels
   from `scripts/audit_table_vi_photon_decays.jl` into `src/` with tests,
   fold Table III mixings into the isoscalar photon rows, wire the model K1
   mixing angle into the Table V decay audit, resolve the strange-parent
   recoil normalization from Appendix B, then extract the remaining Table V
   and Table VI sections.
