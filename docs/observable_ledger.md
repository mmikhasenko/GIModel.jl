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

### Encoded so far

- `data/raw/digitized_tables/table_v_strong_decays/table_v_light_1s_1p.provisional.csv`:
  the light `u-d-s` `1^3S_1`, `1^3P_2`, `1^3P_1`, `1^1P_1`, and `1^3P_0`
  sections of Table V (coefficient expressions, classes, paper numbers,
  experiment column), digitized from the vision OCR Markdown; image audit
  pending for the messier later pages.
- `src/strong_decays.jl`: reduced-amplitude algebra, two-point calibration,
  and the full-amplitude assembly above. Calibration on the paper's fit rows
  gives `A = 1.665`, `S0 = 3.918` (with `beta = 0.40 GeV` and 1984-era
  masses).
- `scripts/audit_table_v_decays.jl` writes
  `docs/residual_reports/table_v_light_decays.md`: of the 70 extracted rows,
  the 33 headline rows with a clean equal-mass two-parameter convention
  reproduce the paper's numeric column at 6% median deviation (e.g.
  `K* -> K pi` +7.97 vs +7.9, `A2 -> eta pi` +4.49 vs +4.5,
  `f -> pi pi` -10.9 vs -11); 8 are flagged and 27 are excluded with stated
  reasons (K1 mixing angle, quasi-two-body lineshapes, strange recoil).

### Open conventions / next steps

1. Decay momenta `q`: computed from a 1984-era mass table (experimental
   masses for established states, GI model masses for unobserved parents like
   `H`, `H'`, `delta2`, `epsilon`, `epsilon'`, `kappa`). Quasi-two-body
   subchannel daughter masses (`(pi pi)_epsilon`, `(K pi)_kappa`,
   `(eta pi)_delta2`) are the dominant ambiguity; rows whose computed
   amplitude deviates from the paper by more than 20% are flagged in the
   audit for a per-row kinematics audit.
2. `Q1`/`Q2` (strange axial) rows: the numeric column folds in the model's
   `K1` (`1^3P_1`/`1^1P_1`) mixing angle on top of the unmixed formula
   coefficients — the near-zero `Q1 -> [K* pi]_S` against the large
   `Q2 -> [K* pi]_S` is the familiar K1(1270)/K1(1400) selectivity. The
   spectrum layer's antisymmetric spin-orbit block provides this angle; wiring
   it into the decay audit is the clearest next step.
3. Strange-parent normalization: `K*2 -> K pi` computes ~sqrt(3) above the
   paper while `K*(892) -> K pi` is exact; the strange `1^3P_0`/`ss`-parent
   `S`-rows sit 20-30% off in a correlated way. Appendix B's
   emission amplitudes carry unequal-mass recoil factors not yet encoded.
4. Later Table V sections (light `2S`, `1D`, charmed, charmonium `psi`
   sectors with `A_c`, `S_c`, `beta_c`) are not yet extracted.
5. "Realistic factor" column (SHO -> realistic wavefunction correction
   ratios) is recorded but not modeled.
6. Table VI/VII electromagnetic amplitudes and charge radii: Table VI is
   started, see the section below. Table VII (leptonic, two-photon, gluonic
   decays, charge radii) is not started.

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
