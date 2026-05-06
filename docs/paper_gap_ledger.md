# Paper Gap Ledger

This is the current short list of what is still missing relative to the
Godfrey-Isgur paper. It is meant to be the durable replacement for older
handoff/autonomous planning notes.

## Implemented Enough For Current Comparisons

- Table II solver inputs are audited into CSV/TOML, including
  `epsilon_so_scalar = +0.055`.
- The semirelativistic kinetic operator is implemented in the finite-difference
  path.
- The central potential uses the GI color-singlet sign and normalization.
- Fig. 2 running `alpha_s(r)` is implemented and regression-tested.
- S-wave contact hyperfine is smeared and diagonalized nonperturbatively in the
  finite-difference path.
- Equal-mass spin-orbit and tensor diagonal shifts are implemented as active
  sector diagnostics.
- Top-level spectrum CSVs cover Figs. 3-9, and promoted clean mass/mixing data
  now exists under `data/clean/`.

## Missing Or Still Approximate

### 1. Eq. (14) Production Basis And Staging

What the paper does: diagonalizes the relativized Hamiltonian in a large
harmonic-oscillator basis, first in fixed `|jm;ls>` sectors, then diagonalizes
smaller mass matrices for tensor, antisymmetric spin-orbit, and annihilation
mixing.

What the repo does: production comparisons use a finite-difference radial grid.
An HO path exists for audit-style checks, but it is not the main path and it is
not yet the paper's staged workflow.

Why it matters: heavy-heavy masses are already close, so FD is a useful
numerical representation. But a literal paper reproduction should be able to
separate numerical-basis differences from physics differences.

Acceptance check: run the same sector through FD and HO production paths, with
identical model switches, and report whether residual differences are smaller
than the physical effects being studied.

### 2. Full Appendix-A Relativization

What the paper does: replaces naive mass factors in spin-dependent and smeared
operators by operator-level `m/E` factors with fitted epsilon exponents and
momentum-space sandwiching.

What the repo does: several smearing and momentum-sandwich ingredients are
implemented, and the central path has a GI-style momentum sandwich. Some active
spin-dependent pieces still use scalar `(1 + epsilon)` multipliers as a
controlled approximation.

Why it matters: scalar epsilon factors can tune signs and scales, but they do
not reproduce state-dependent momentum weighting. The gap matters most for
light states and radial excitations.

Acceptance check: for contact, tensor, vector spin-orbit, and scalar
spin-orbit terms, document the exact Appendix-A operator used and add a
regression comparing the operator-level implementation to the current scalar
approximation on at least `ccbar`, `ssbar`, and light `q qbar`.

### 3. Off-Diagonal Tensor Mixing

What the paper does: the tensor part of the hyperfine interaction mixes
`^3L_J` and `^3L'_J` states with the same `J` after the first fixed-sector
diagonalization.

What the repo does: tensor diagonal shifts and angular factors are tested, but
sector comparisons do not build the physical same-`J` tensor mass matrix for
reported rows.

Why it matters: this affects states such as triplet S/D candidates and can move
assignments, not just masses. Current reports mark these rows as mixing-prone
rather than judging them as final.

Acceptance check: add an explicit tensor mixing block, expose eigenvectors, and
exclude fewer rows from the non-mixing scorecard only after the block has tests.

### 4. Unequal-Mass Antisymmetric Spin-Orbit Mixing

What the paper does: for unequal constituent masses, the antisymmetric
spin-orbit term mixes `^1L_J` and `^3L_J` states.

What the repo does: same-`J` spin-orbit mixing diagnostics exist, but
heavy-light, strange, charmed, and bottom-flavored sector reports still use a
diagnostic equal-share fine-structure convention for physical comparisons.

Why it matters: heavy-light P-wave patterns are not meaningful final tests
until this is implemented. Without it, residuals can look like potential
failures when they are really assignment/mixing failures.

Acceptance check: implement a mass-matrix path for `^1L_J`/`^3L_J` pairs,
report mixing angles/eigenvectors, and update residual assignment logic for
open-flavor sectors.

### 5. Literal Isoscalar Annihilation, Eq. (16)-Eq. (18)

What the paper does: for self-conjugate isoscalars, it adds an annihilation
mass matrix after the fixed-sector diagonalization. Eq. (16)-Eq. (17) define
the general matrix element; Eq. (18a) and Eq. (18b) replace the pseudoscalar
bracket for P1 and P2.

What the repo does: formulas and Table III pseudoscalar rows are audited.
`:calibrated_p1` is a rank-one diagnostic that reproduces the four isoscalar
pseudoscalar masses, but literal `:paper_p1` and `:paper_p2` modes do not exist.

Why it matters: calibrated P1 proves the missing physics is localized, but it
is not a paper implementation. P2 is especially different because the paper
expects mass-dependent, non-orthogonal poles.

Acceptance check: implement `:paper_p1` and `:paper_p2` as separate dispatch
modes, keep `:calibrated_p1` as a control, and require scorecard rows for all
four modes.

### 6. General Table III Isoscalar Mixing

What the paper does: Table III gives approximate isoscalar compositions beyond
the pseudoscalar rows, including non-pseudoscalar sectors with predicted and
observed splittings.

What the repo does: `data/clean/mixings.csv` promotes only the visually checked
pseudoscalar P1/P2 rows. Other Table III rows remain raw/provisional.

Why it matters: once general annihilation is added, non-pseudoscalar isoscalar
rows should become validation targets instead of being treated as unexplained
duplicates in the isoscalar residual report.

Acceptance check: promote the remaining image-audited Table III rows into
clean data and add a validation report for composition and mass splitting.

### 7. Eigenvector Fidelity Scoring

What the paper gives: Table III amplitudes are eigenvector/composition targets,
not just mass targets.

What the repo does: `docs/residual_reports/annihilation_model_scorecard.md`
scores formula provenance, clean targets, implementation modes, and
pseudoscalar mass residuals. It does not compute RMS amplitude error yet.

Why it matters: a model can fit the four pseudoscalar masses while producing
the wrong flavor/radial composition. The next score must prevent that false
positive.

Acceptance check: compare model eigenvectors against `data/clean/mixings.csv`
up to overall sign, report RMS amplitude error separately for P1 and P2, and
award mixing-fidelity points only from that calculation.

### 8. Observables Outside The Mass Spectrum

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

Implement literal `:paper_p1` first. It is the narrowest useful paper mode:
Eq. (16) and Eq. (17) are shared with P2, Eq. (18a) is simpler than Eq. (18b),
and the promoted P1 pseudoscalar rows already provide both mass and composition
targets. After that, add eigenvector RMS scoring before touching P2.
