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
6. Table VI/VII electromagnetic amplitudes and charge radii: not started.

## Acceptance

A Table V section counts as reproduced when every row with an unambiguous
kinematic convention matches the paper's numeric amplitude to ~10% (the
paper itself rounds to 2 significant figures), with the two fit rows exact by
construction.
