# Observable Ledger

> **Detail layer.** Per-unit status for the decay/EM units (Eqs. 19-22, Tables
> IV-VI, Appendices B-D) is the manifest/dashboard (`docs/paper_manifest/*.toml`;
> `cd docs && make dashboard`). This file is the **conventions** detail behind
> them: the Table IV class algebra, the `[PAPER]`/`[DERIVED]` provenance split,
> the leading-S₀ convention, and the photon-decay overlap kernels — the
> physics-convention reasoning the per-unit `notes` only summarize.

Scope ledger for GI observables beyond the mass spectrum (gap-ledger item 8).
The mass-spectrum reproduction machinery stays in `src/` untouched; decay and
electromagnetic conventions are tracked here before/while they are encoded.

## Strong Decays `M* -> M + P` (Sec. IV, Tables IV-V, Appendices B-C)

### Model structure

- Eq. (19): pseudoscalar emission from a quark/antiquark with two couplings
  `g` (sigma.q term) and `h` (sigma.p' recoil term), flavor operators
  `X^i_{q(qbar)}` from Appendix B, helicity amplitudes converted to partial
  waves via Appendix C.
- The paper evaluates amplitudes analytically with single-`beta`
  harmonic-oscillator wavefunctions (the SU(6) limit), `beta = 0.40 GeV`
  (Table VI footnote g confirms "as in Table V").
- Table IV reduced partial-wave amplitudes. Structure-independent classes
  `A = [g + h/4]beta`, `A' = [g - h/4]beta`, `A'' = [g + h/8]beta`,
  `A0 = [g]beta`, `A_c`; structure-dependent classes
  `S = [3h - (1/2)(g + h/4) q^2/beta^2]beta`,
  `D = [3h - (3/10)(...) q^2/beta^2]beta`,
  `P = [3h - (3/4)(...) q^2/beta^2]beta`, `S_c`.
  Writing `S0 = 3 h beta`, the structure-dependent amplitudes are
  `S(qbar) = S0 - (1/2) A qbar^2`, `D = S0 - (3/10) A qbar^2`,
  `P = S0 - (3/4) A qbar^2` with `qbar = q/beta`.
- Two-parameter fit: the paper sets `A' = A'' = A0 = A` numerically (rows
  display `(A0/A)`-style factors) and uses one `S`-type strength. The fit
  inputs are `rho -> pi pi = +12.4 MeV^(1/2)` (fixes `A`) and
  `B -> [omega pi]_S = -11 MeV^(1/2)` (fixes `S0` once `A` is known).
- Full amplitude for a Table V row with coefficient `c`, class `X`, and
  orbital power `L`:
  `amp = c * X(qbar) * qbar^L * sqrt(q/(2 pi)) * exp(-q^2/(16 beta^2))`
  in `MeV^(1/2)` with `q` in MeV inside the square root. The table caption
  confirms the suppressed factor `(q/2pi)^(1/2) exp(-q^2/16 beta^2)` and that
  identical-particle `1/sqrt(2)` factors are already inside `c`.
- Isoscalar rows use Table III mixings (P1 for pseudoscalars; `(B14)/(B15)`
  perfect mixing for the `1^1S_0` formula column).

### Equation provenance ([PAPER] transcribed vs [DERIVED] in-module)

An amplitude factorizes as `amp = c * X(qbar) * spatial_overlap`
(matrix element x numerical overlap):

| factor | symbol | provenance | where |
|---|---|---|---|
| flavor-spin coefficient | `c` | **[PAPER]** App. B, one per row | canonical CSV `coefficient` |
| reduced-amplitude class algebra | `X` form, `k = r,3/10,3/4` (`r = m_Q/(m_Q+m_q)`, `1/2` at equal mass) | **[PAPER]** Table IV | `reduced_decay_amplitude` |
| fitted strengths | `A`, `S0` | **[DERIVED]** solved from 2 fit rows | `calibrate_strong_decay_model` |
| spatial overlap | `qbar^L sqrt(q/2pi) e^{-(1/4) r^2 q^2/b_c^2}` (`r=1/2` is the light `e^{-q^2/16b^2}`) | **[DERIVED]** SHO integral | `spatial_overlap` |
| breakup momentum | `q` | **[DERIVED]** Kallen | `decay_momentum` |
| unequal-mass form factor + recoil | footnote d | **[PAPER]** shape, **[DERIVED]** value | `spatial_overlap(; heavy_fraction, recoil)` |
| partial width | `|amp|^2` | **[DERIVED]** GI normalization | `decay_width` |

**Leading-S0 convention.** The paper's structure-dependent (S/D/P) numeric
column uses the *leading constant* `S0 = 3 h beta` (dropping the
`-k A qbar^2` polynomial of Table IV). That is what reproduces the reported
`S0 ~ 3.29` and the D/P rows to ~1%; it is the `:leading` convention (default
of `decay_amplitude`), while `:table_iv` keeps the printed polynomial.

### Encoded so far

- `data/raw/digitized_tables/table_v_strong_decays.csv`: the canonical,
  page-image-verified transcription of all of Table V (220 rows, see
  `README_canonical_tables.md`).
- `src/strong_decays.jl`: the row-oriented API — `DecayChannel`,
  `decay_amplitude` returning the 3-way `StrongDecayAmplitude` decomposition
  (`coefficient`/`reduced`/`spatial_overlap`/`total`), `matrix_element`,
  `decay_width`, `MesonMasses`, the `:leading`/`:table_iv`
  conventions, and the charm footnote-d path. Leading calibration gives
  `A = 1.665`, `S0 = 3.287` (`beta = 0.40 GeV`).
- `scripts/reproduce_table_v.jl` writes
  `docs/residual_reports/table_v_reproduction.md`: **160 / 178 scoreable rows
  match** (+6 near, 12 off), 22 convention-deferred. Non-matches are all
  parent-mass or mixing-angle input sensitivity (each MISS inverts to a mass
  within 20-50 MeV of the input), not decay-algebra error.
- `GIPaper.load_table_v` owns the adapter from the paper's canonical CSV schema
  to GIModel's reusable `DecayChannel` type.

### Findings that make Table V reproduce

1. **Leading-S0** turns every S/D/P row from "excluded" to a ~1% match
   (`rho -> [omega pi]_P` -7.82 vs -7.8).
2. **`K*2 -> K pi` sqrt(3) was a digitization error, not physics.** The old
   provisional coefficient `+(3/20)^1/2` should be `+(1/20)^1/2` (page image);
   the canonical value gives +7.60 vs paper +7.7. The careful re-digitization
   fixed it.
3. **`h` (f4) is the nn isoscalar** (decays to pi pi, K Kbar), not ss; a flavor
   fix that reproduced its whole 1^3F_4 block.
4. **Same-J mixing** (Q1/Q2, Q1c/Q2c): rotating the pure singlet/triplet
   formulas by the paper's angle (1P +34, 1D +33, charm -41 deg) reproduces the
   footnote-j near-cancellations in sign and magnitude
   (`Q1->[K*pi]_S` -0.34 vs -0.3, `Q2->[K*pi]_S` +17.4 vs +16).

### Open conventions / next steps

1. Quasi-two-body subchannel daughters (`(pi pi)_eps`, `(K pi)_kappa`,
   `(eta pi)_delta2`) and sub-threshold modes depend on a lineshape the paper
   does not state; those rows are convention-deferred, not scored.
2. GI-predicted parent masses: the residual MISSes (isoscalar `H`/`H'`, `2S`
   `K'`/`rho_S`, near-threshold `K* Kbar`) track the parent mass; using the
   paper's own predicted masses would close them. Documented per-row via the
   implied-mass inversion in the report.
3. "Realistic factor" column (SHO -> realistic wavefunction ratios) is recorded
   in the CSV but not applied (the leading-S0 finding superseded the earlier
   "realistic factor folded in" hypothesis).
4. Applying Eq. (19) directly to the calculated physical waves would be a
   beyond-paper extension. The original numerical single-beta SHO treatment is
   reproduced by the Table IV/V path above.

## Photon Decays `M* -> M gamma` (Sec. IV B, Table VI, Appendix D)

### Model structure

- Appendix D mock-meson matrix elements, evaluated on the model's own
  wavefunctions (unlike Table V's single-`beta` SHO limit):
  `I_i(x,y)` is a momentum-space overlap with weight
  `(1/m_i)(m_i/E_i)^0.7`, `E_n^i(x,y)` is a position-space `r^n` moment with
  prefactor `|m_i / sqrt(<E_i>_x <E_i>_y)|^0.5`, and the mock mass is
  `M~ = <E_1> + <E_2>`. The exponents 0.7/0.5 are the paper's, fitted there
  to `rho -> pi gamma` and `A2 -> pi gamma` — so no new fitted constants
  enter on our side.
- M1 moments are listed in units of `e/2`, hence `mu/mu_N = coeff * I * M_N`.
- E1 amplitudes are `formula * sqrt(alpha q)` in `MeV^(1/2)` with the
  q-dependence explicit in the formula column.
- Open-flavor M1 coefficients follow `mu = e_q I_q - e_qbar I_qbar` (the
  antiquark charge enters flipped), which reproduces every printed
  quarkonium coefficient (`+4/3 I_c` for `psi`, `-2/3 I_b` for `Upsilon`).

### Encoded so far

- `GIPaper/scripts/audit_table_vi_photon_decays.jl` writes the complete
  79-row report and a machine-readable CSV: 42 M1, 35 E1, and 2 M2.
  Every state uses the shared native fixed-channel solver, including the
  spin-distorted P waves. Four-flavor P1/vector/tensor states are composed
  through the final MixedSpectrum, with pure-flavor charges applied once.
- The canonical loader rejects duplicate identities, preserves parenthetical
  and approximate predictions, and separates state/kinematics inputs from
  computed observables.
- Recoil, the footnote-g form factor, the explicit footnote-a +0.01 μN
  input, and M1 width conversion are covered. A1 -> pi is spin-flip E1;
  its q² power is not an M2 classification.
- Excited eta and near-cancelled amplitudes retain explicit quantitative
  residuals. The independent FD run and solver-comparison report distinguish
  these from representation errors. There are no omitted Table VI rows.

## Numerical method (Appendix A, Eq. A17)

Two independent algorithms solve the same radial problem, behind one interface:

| | `FiniteDifferenceSolver` | `OscillatorSolver` |
|---|---|---|
| representation | `u(r)` on a uniform mesh | adaptive oscillator expansion (starts at 24 states) |
| `⟨i\|f(p)\|n⟩` | spectral function of the FD `p²` | `ho_p2_matrix`, closed form |
| `⟨n\|g(r)\|j⟩` | mesh quadrature | `ho_operator_matrix`, Gauss–Laguerre |
| spatial mesh | yes, intrinsically | **none** (reporting only) |

The oscillator path is the paper's own: Eq. (A17) inserts a complete set between
`f(p)` and `g(r)` so each factor is evaluated in the space where it is simple,
which is cheap precisely because oscillator functions are the same
polynomial×Gaussian in both. It is now built that way rather than imitated on a
grid.

**Cross-checks that need no reference data.** `p²` and `r²` reconstruct the
oscillator Hamiltonian exactly diagonal (1e-14); the virial theorem holds per
state; `g = 1` returns the identity (4e-15); `g = r²` returns `ho_r2_matrix`
(1e-12); and the smeared Appendix-A potential matches independent adaptive
quadrature to 2.8e-16…1.2e-13. The production controller now enlarges the basis
until all requested energies pass a recorded 0.1 MeV convergence criterion.

**Agreement between the two.** Charm +0.18 MeV, bottom +0.55 MeV, oscillator
above finite-difference — the correct side for a variational calculation in a
finite basis. Light sectors sit ~1.5 MeV the other way, which reads as the FD
mesh being high where short-distance structure is hardest to resolve, not as a
disagreement about the physics.

**Against the paper.** With the oscillator path running the paper's own method,
the mean absolute residual against GI Tables I/II improves in every sector:

| sector | FD | HO | Δ |
|---|---:|---:|---:|
| isovector | 14.54 | 14.13 | −0.41 |
| isoscalar | 138.65 | 138.18 | −0.47 |
| strange | 10.85 | 10.70 | −0.15 |
| charmonium | 6.03 | 5.92 | −0.11 |

(MeV; the isoscalar absolute value is dominated by the annihilation-mixing gap
tracked separately, so read its Δ rather than its level.) The improvement is
small and systematic — four of four sectors, none worse — which is the expected
signature if part of the paper's own residual is its numerics rather than its
physics. It is not a large enough effect to change any conclusion about the
model; it is evidence that we are now running GI's method and not merely
matching its numbers by another route.

Comparator central methods (pointwise, 1D/3D-smeared, derivative-G) keep the
mesh on both sides deliberately: several smear numerically on it, so they are
not closed-form functions of `r`, and a hybrid Hamiltonian — exact kinetic,
mesh-projected potential — is not variational.

## Acceptance

A Table V section counts as reproduced when every row with an unambiguous
kinematic convention matches the paper's numeric amplitude to ~10% (the
paper itself rounds to 2 significant figures), with the two fit rows exact by
construction.

A Table VI block counts as reproduced when all 79 canonical rows are computed
or explicitly classified, with formula, phase, mass, and mixing conventions
itemized. Mixing-dependent rows must use the shared Table III physical-state
composition. Numerical disagreements remain reported residuals; coverage does
not require every cancellation-sensitive or excited-state amplitude to agree
within 10%.
