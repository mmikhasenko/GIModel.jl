# Canonical digitized tables (Godfrey–Isgur 1985)

These CSVs are a **from-the-page-image** transcription of the paper's tables,
verified row-by-row against `paper/vision_ocr/page_images/page-0NN.png` (the
`page` column is the PDF page = image number; `journal_pg` is the printed page).
They supersede the earlier piecemeal `table_v_*.provisional.csv` files, which
conflated the two "realistic factor" sub-columns and carried a few OCR sign/value
errors (see below).

## `table_v_strong_decays.csv` — Table V, strong decays M* → M + P (pages 14–22, 220 rows)

Column meaning (transcribed faithfully; interpretation columns flagged):

| column | meaning |
|---|---|
| `page`, `journal_pg` | PDF page (= image) and printed journal page |
| `section` | nonet header printed in the table (`n^(2S+1)L_J` + strange/nonstrange) |
| `parent` | decaying meson, disambiguated by state (e.g. `Kstar2` = 1³P₂ K, `rho2` = 1³D₂ ρ) |
| `decay` | decay label as printed |
| `daughter1`, `daughter2` | **derived** daughters for kinematics (quasi-two-body `(AB)_M` → the resonance M) |
| `formula` | HO amplitude formula as printed (`qt` = q̃ = q/β) |
| `coefficient` | **derived** numeric value of the flavor/spin coefficient |
| `amp_class` | reduced-amplitude class from the formula letter: `A`,`Aprime`,`Adoubleprime`,`A0`,`S`,`D`,`P`,`Ac`,`Sc`,`mixing_only`,`unlisted` |
| `qbar_power` | explicit q̃^L power |
| `amp_MeV` | **the tabulated (MeV^1/2) numeric value** (string: allows `below_threshold`, `~0`, two-part `-0.11-1.1`) |
| `realistic_bracket` | the class-ratio annotation in the "Realistic factor" column: `[D/S]`, `[P/S]`, `[A0/A]`, `[A'/A]`, `[A''/A]`, `[Ac/A]`, `[Sc/S]` |
| `realistic_factor` | the **numeric** realistic multiplier in that column: `(1.x)`/`(2.1)` |
| `experiment` | experiment column (`see Sec. VA`, `x +- y`, `not seen`, `seen`, `~14`) |
| `footnotes` | footnote letters (quoted when multiple: `"b,f,i"`) |
| `confidence` | `high` (image-verified) / `medium` / `low` (flagged ambiguity in `notes`) |
| `notes` | per-row provenance / caveats |

**Key point (`realistic_bracket` vs `realistic_factor`):** the paper's "Realistic
factor" column holds *two* different things — a class-ratio **bracket** like
`[D/S]` (a display annotation, per the Table-caption text "we calculate with A
and S but explicitly display factors of A'/A, D/S, etc.") and a numeric
**multiplier** like `(1.7)`. The old provisional files stored only one field and
mixed them; here they are separate.

**`amp_MeV` is NOT raw-HO for the structure-dependent (S/D/P) rows.** The paper
computes the S/D/P numeric column with the *leading* constant reduced amplitude
S₀ = 3hβ (dropping the −k·A·q̄² polynomial); this reproduces both the reported
S₀=3.27 fit and the tabulated D/P numbers to ~1% (see the audit). The `formula`
column still shows the full class letter (`D q̃`, `P q̃`).

### Errors in the old provisional files that this transcription fixes
- `h → η η` (1³F₄): old `+1.0`, page shows `−(1/241920)^½ → −0.6`.
- `h' → η η` (1³F₄): old `+1.4`, page shows `−(1/120960)^½ → −1.4`.
- `K* → K π` (1³F₄, K₄*): old had `+` coeff with a spurious "sign mismatch" note;
  page shows `−(1/40320)^½ → −2.7` (signs consistent, no mismatch).

### Genuine on-page sign mismatches (formula sign ≠ printed value sign), flagged `low`
- `Q2 → [K*η]_P` (1³D₂ strange): formula `+`, printed value `−1.6`.
- `omegaS → (ππ)_ρ π` (2³S₁): formula `+`, printed value `−5.8`.

## Table V footnotes (page 210)
- **a** From the Particle Data Group, Rev. Mod. Phys. 56, S1 (1984), unless otherwise noted.
- **b** We quote the pure Q_B(A) [pure singlet (triplet)] amplitude formulas under Q_1(2), where Q_1(2) is the lower (higher) state in mass.
- **c** We quote the pure ¹D₂ (³D₂) amplitude formulas under Q_1(2), lower (higher) state in mass.
- **d** A charmed P-wave meson is denoted by its I=½ u-d-s analog with subscript c; β_c (≡ β numerically) is the charmed HO parameter; the form factor is exp[−¼ (m_c/(m_c+m_d))² q²/β_c²].
- **e** δ₂ is the qq̄ state at ~1100 MeV decaying to ηπ and KK̄ (widths ~200, ~125 MeV); the observed KK̄π and ηππ are assumed to arise from this channel.
- **f** Final-state-particle width allowed for; notation M* → (AB)_M P for M* → MP → ABP.
- **g** Reference 14.
- **h** Reference 9.
- **i** Reference 15.
- **j** This result is extremely sensitive to the Q_A–Q_B mixing angle of Fig. 4.

## `table_vi_photon_decays.csv` — Table VI, photon decays (pages 24–26, 79 rows)

Two blocks, distinguished by the `multipole` column:
- **M1** (magnetic-dipole, 42 rows): `formula` is the magnetic moment μ in units of
  e/2; `predicted` and `experiment` are μ in units of μ_N. Parenthetical predicted
  values (e.g. `(+0.06)`) are mixing/forbidden-induced order-of-magnitude estimates
  (footnote d). Hindered quarkonium rows carry a `q^2/(24 m_Q) E2` recoil term.
- **E1 / M2** (37 rows, "Other electric and magnetic multipole decays"): `formula`
  is the amplitude in units of √(αq); `predicted`(="Theory") and `experiment` are in
  MeV^1/2. `multipole` is classified by the q-power in the formula (`q^L/m` → M2, `q`
  → E1); the `I_m`/`E_n^i` factors are defined in the Table VI caption (Appendix D).

Columns: `page,journal_pg,multipole,decay,parent,daughter,formula,predicted,experiment,footnotes,confidence,notes`.
`formula` fields are quoted (they contain commas from `I_x(a,b)`). `daughter` is the
non-photon final meson.

### Table VI footnotes (pages 26–27)
- **a** This moment includes a contribution of +0.01 from π⁰–η mixing.
- **b** Absolute widths unknown, but predicted branching ratios compare favorably to experiment via Table V.
- **c** Recoil term included in this decay since the direct term is so small (hindered M1).
- **d** Amplitudes from photon emission off the non-ss̄ (cc̄,bb̄) component of the initial φ(ψ,Υ) or to the ss̄(cc̄,bb̄) component of the final meson (Table III); Zweig-violating contributions are omitted, so the parenthetical "predicted" values are order-of-magnitude only.
- **e** We tentatively identify η_r with the ι(1440); see Sec. VA.
- **f** F_b denotes a bs̄ meson.
- **g** Recoil absorbed by a light-quark system, so a form factor exp(−q²/16β²) with β=0.40 GeV is included (as in Table V).
- **h** Reference 20.
- **i** Reference 23.

### GI-predicted masses for the 2S/1D vector states (Sec. VA, page 215), used as Table V/VII parents
ρ_S 1.45, ω_S 1.46, K*_S 1.58, φ_S 1.69 (2³S₁); ρ_D 1.66, ω_D 1.66, K*_D 1.78, φ_D 1.88 (1³D₁) GeV.

## `table_vii_annihilation_em.csv` — Table VII, leptonic/γγ/gluonic decays + charge radii (pages 28–30, 61 rows)

Four sub-tables via the `subtable` column:
- **leptonic** (27 rows): `quantity` is the decay constant (f_π, f_ρ, …); `formula`
  is the amplitude (e.g. `2 sqrt3 P_pi`, `sqrt6 V_rho`) in terms of the P_P/P'_A1/
  V_V/V'_V factors defined in the caption; `predicted`/`experiment` are the amplitude
  (units as in the paper).
- **gamma_gamma** (13 rows): P→γγ / ³P₂→γγ amplitudes; `predicted`/`experiment` carry
  their unit inline (`2.6 eV^1/2`, `0.50 keV^1/2`) because it varies by row.
- **gluonic** (18 rows): QQ̄→2g/3g amplitudes in MeV^1/2 (noted per row).
- **charge_radius** (3 rows): `predicted`/`experiment` are r_E² in fm², written as the
  paper does, `±(radius)^2` (sign = sign of r_E²).

Columns: `page,journal_pg,subtable,decay,quantity,formula,predicted,experiment,footnotes,confidence,notes`.
`ζ`/`η_t` rows are hypothetical t-t̄ mesons (kept for completeness, flagged in notes).

### Table VII footnotes (page 218)
- **a** Reference 23.
- **b** This amplitude is very sensitive to f–f' mixing.
- **c** Reference 24.
- **d** Reference 25.
- **e** Reference 26.
- **f** Reference 27.
