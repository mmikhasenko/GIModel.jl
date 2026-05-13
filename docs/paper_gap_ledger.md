# Paper Gap Ledger

This is the current short list of what is still missing relative to the
Godfrey-Isgur paper. It is meant to be the durable replacement for older
handoff/autonomous planning notes.

## Implemented Enough For Current Comparisons

- Table II solver inputs are audited into CSV/TOML, including
  `epsilon_so_scalar = +0.055`.
- The semirelativistic kinetic operator is implemented in the finite-difference
  path.
- `HarmonicOscillatorBasis` is implemented through `GIParameters{Basis}`
  dispatch and has regression tests. `scripts/audit_nonmixing_contact.jl`
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
  pseudoscalar P1/P2 modes driven by the `[annihilation]` TOML constants.
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

What the repo does: both FD and HO radial basis paths exist. The HO path is
selected through `GIParameters{HarmonicOscillatorBasis}`, projects operators
into a finite oscillator subspace, scans the oscillator scale, reconstructs
mesh wavefunctions, and is tested. Current headline sector reports still use
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

What the repo does: formulas and Table III rows are audited.
`:calibrated_p1` is routed through `assign_mixed_rows(::IsoscalarAnnihilation, ...)`
as a rank-one control that reproduces the four isoscalar pseudoscalar masses.
Literal `:paper_p1` and `:paper_p2` modes exist and their constants are
represented in `GIParameters` and the TOML parameter file. The formula path now
uses the paper's momentum-space `alpha_s(Q^2)`, evaluates the S-wave Eq. (17)
momentum integral on the FD radial wavefunction, and includes the coherent
`sqrt(2)` factor for the normalized `ns` flavor state. These local corrections
still do not reproduce the Table III eigenvectors; see
`docs/residual_reports/table_iii_mixing_audit.md` and
`docs/residual_reports/table_iii_suspect_investigation.md`.

Why it matters: calibrated P1 proves the missing physics is localized, but it
is not a paper implementation. P2 is especially different because the paper
expects mass-dependent, non-orthogonal poles.

Acceptance check: partially complete. The modes, parameters, momentum-coupling
convention, FD Eq. (17) integral, and `ns` flavor factor exist; the remaining
acceptance criterion is eigenvector/mass-splitting fidelity against Table III,
likely requiring paper-order HO basis inputs or a sharper reconstruction of the
pseudoscalar phenomenological prescription.

### 6. General Table III Isoscalar Mixing, Later Than 1-5

What the paper does: Table III gives approximate isoscalar compositions beyond
the pseudoscalar rows, including non-pseudoscalar sectors with predicted and
observed splittings.

What the repo does: `data/clean/mixings.csv` promotes the visually checked
pseudoscalar P1/P2 rows plus the visible `1^3S_1` and `1^3P_2` rows from
page 11.

Why it matters: once general annihilation is added, non-pseudoscalar isoscalar
rows should become validation targets instead of being treated as unexplained
duplicates in the isoscalar residual report.

Acceptance check: clean-data promotion is complete for the visible rows.
Remaining work is a validation report for the non-pseudoscalar Eq. (16)
composition and mass splitting.

### 7. Eigenvector Fidelity Scoring, Later Than 1-5

What the paper gives: Table III amplitudes are eigenvector/composition targets,
not just mass targets.

What the repo does: `docs/residual_reports/annihilation_model_scorecard.md`
scores formula provenance, clean targets, implementation modes, and
pseudoscalar mass residuals. `docs/residual_reports/table_iii_mixing_audit.md`
now computes RMS amplitude error for literal P1/P2 against the promoted
Table III pseudoscalar rows.

Why it matters: a model can fit the four pseudoscalar masses while producing
the wrong flavor/radial composition. The next score must prevent that false
positive.

Acceptance check: amplitude comparison exists. Next, fold that RMS into the
annihilation scorecard and fix the matrix-element path until the score is
acceptable.

### 8. Observables Outside The Mass Spectrum, Later Than 1-5

What the paper includes: leptonic, two-photon, gluonic decay amplitudes, charge
radii, and decay-model discussion beyond the mass spectrum.

What the repo does: current scope is masses, residual reports, and the inputs
needed for the mass Hamiltonian.

Why it matters: these observables test wavefunctions, not only eigenvalues.
They should be a later validation layer after mass-spectrum mechanics settle.

Acceptance check: create a separate observable ledger before implementing them,
so mass-reproduction code does not absorb decay/charge-radius conventions
implicitly.

## Next Clean Implementation Step

The most direct implementation order is:

1. Use `table_iii_mixing_audit.md` to debug why literal FD P1/P2 under-mix
   compared with Table III.
2. Implement the general non-pseudoscalar Eq. (16) model for `1^3S_1` and
   `1^3P_2`.
3. Fold Table III eigenvector/mass-shift scoring into
   `annihilation_model_scorecard.md`.
4. Revisit Appendix-A/HO paper-order staging with the now-explicit mixing
   blocks.
