# Observable Ledger

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
| reduced-amplitude class algebra | `X` form, `k = 1/2,3/10,3/4` | **[PAPER]** Table IV | `reduced_decay_amplitude` |
| fitted strengths | `A`, `S0` | **[DERIVED]** solved from 2 fit rows | `calibrate_strong_decay_model` |
| spatial overlap | `qbar^L sqrt(q/2pi) e^{-q^2/16b^2}` | **[DERIVED]** SHO integral | `spatial_overlap` |
| breakup momentum | `q` | **[DERIVED]** Kallen | `decay_momentum` |
| charm form factor + recoil | footnote d | **[PAPER]** shape, **[DERIVED]** value | `spatial_overlap(; charm)` |
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
  `decay_width`, `MesonMasses`, `load_table_v`, the `:leading`/`:table_iv`
  conventions, and the charm footnote-d path. Leading calibration gives
  `A = 1.665`, `S0 = 3.287` (`beta = 0.40 GeV`).
- `scripts/reproduce_table_v.jl` writes
  `docs/residual_reports/table_v_reproduction.md`: **160 / 178 scoreable rows
  match** (+6 near, 12 off), 22 convention-deferred. Non-matches are all
  parent-mass or mixing-angle input sensitivity (each MISS inverts to a mass
  within 20-50 MeV of the input), not decay-algebra error.

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
4. Table VII (leptonic, two-photon, gluonic decays, charge radii) audit is not
   started; Table VI photon decays are in the section below.

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

- `scripts/audit_table_vi_photon_decays.jl` writes
  `docs/residual_reports/table_vi_photon_decays.md`: 28 mixing-free rows
  computed from the FD solver wavefunctions (contact-distorted S waves,
  central P waves). Quarkonium M1 rows land at the 0.1-2% level
  (`psi -> eta_c gamma` +0.684 vs +0.69, `psi' -> eta_c' gamma` +0.680 vs
  +0.68, `Upsilon` family -0.121/-0.121/-0.120 vs -0.13/-0.12/-0.12),
  open-flavor M1 at 1-4% (`D*+` -0.347 vs -0.35, `F*` -0.132 vs -0.13,
  `B*+` +1.360 vs +1.37, `F_b*` -0.550 vs -0.55), light rows at ~6%
  (`rho -> pi gamma` +0.650 vs the +0.69 fit target — the known
  light-sector wavefunction residual), and the hindered
  `psi' -> eta_c gamma` with the recoil term gets sign and magnitude
  (-0.067 vs -0.056). E1 `chi_c`/`chi_b` triplets reproduce at 3-7%; the
  two `2S -> chi_0` rows sit 20-30% high (largest q, node cancellation).

### Open conventions / next steps

1. The vision-OCR predicted column in the open-flavor M1 block is displaced
   by one row against the decay labels; the audit uses shift-corrected
   values (confirmed by 0.1-2% matches on four independent rows, and the
   orphaned `-0.55` landing exactly on the computed `F_b*`), but a crop
   audit of the printed PDF column should confirm the alignment.
2. `2S -> chi_0` E1 rows: check whether the paper used model masses rather
   than measured 1984 masses for the photon momentum `q`.
3. Isoscalar rows (`phi -> eta gamma`, `eta' -> rho gamma`, ...) need the
   Table III mixing amplitudes folded in — the mixing layer already
   provides them.
4. Remaining Table VI blocks: light E1/M2 rows (`A2 -> pi gamma` is the 0.5
   exponent fit row), strange/charmed P-wave rows, hindered bottomonium
   rows, and the `psi/Upsilon -> (light) gamma` order-of-magnitude rows
   (footnote d).
5. Promote the overlap kernels (`I_i`, `E_n^i`, mock mass, momentum waves)
   from the audit script into `src/` with regression tests once the
   conventions above settle.

## Acceptance

A Table V section counts as reproduced when every row with an unambiguous
kinematic convention matches the paper's numeric amplitude to ~10% (the
paper itself rounds to 2 significant figures), with the two fit rows exact by
construction.

A Table VI block counts as reproduced when every mixing-free row matches the
paper's moment/amplitude to ~10% with conventions itemized, and
mixing-dependent rows additionally use the Table III amplitudes from the
mixing layer.
