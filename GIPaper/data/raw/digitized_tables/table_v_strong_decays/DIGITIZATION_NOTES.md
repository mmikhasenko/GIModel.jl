# Table V (Strong Decays) — Digitization Notes

> **SUPERSEDED.** The `*.provisional.csv` files in this folder have been folded
> into the single canonical, page-image-verified
> [`../table_v_strong_decays.csv`](../table_v_strong_decays.csv) (220 rows; see
> [`../README_canonical_tables.md`](../README_canonical_tables.md)). The
> canonical transcription also **corrected three sign/value errors** in the old
> `1³F₄` rows here (`h→ηη`, `h'→ηη`, `K₄*→Kπ`) and split the conflated
> `[D/S]`/`(1.x)` realistic-factor columns. These files are retained only until
> `audit_table_v_decays.jl` and `audit_table_v_2s_1d_decays.jl` are rewired to
> the canonical CSV; do not use them for new work.

## Charmed / open-charm section (unit W2a-3)

### Source
- OCR: `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`, lines 965–993
  ("Charmed mesons" blocks under `## PDF Page 22`, `TABLE V. (Continued)`).
- Image verified against `paper/vision_ocr/page_images/page-022.png`
  (printed page **210**, journal page number `32` top-right).
- Output: `table_v_charmed.provisional.csv` (10 rows).

### Row count / structure
The charmed section is 10 decays across four spin blocks, all on page 210:
- `1³S₁`: D*⁺→D⁰π⁺, D*⁺→D⁺π⁰, D*⁰→D⁰π⁰ (3 rows)
- `1³P₂`: K*_c→Dπ, K*_c→D*π (2 rows)
- `1³P₁ and 1¹P₁`: Q_1c→[D*π]_S, Q_1c→[D*π]_D, Q_2c→[D*π]_S, Q_2c→[D*π]_D (4 rows)
- `1³P₀`: κ_c→Dπ (1 row)

### IMPORTANT: OCR omissions recovered from the image
The vision-OCR markdown for this page captured **only two columns**
(Decay, HO amplitude). The page image actually has the full 6-column
Table V layout. The following columns were read directly off `page-022.png`
and are NOT present in the OCR text:
- **(MeV¹ᐟ²)** numeric amplitude column → `paper_MeV12`
- **"Realistic" factor** column → `realistic_factor`
- **References, footnotes** column → `footnotes`

### A_c / S_c / β_c convention as PRINTED (feeds W2c)
- The charmed amplitudes use **A_c** (spatial/D-type class) and **S_c**
  (S-type class) in place of the light-meson **A** and **S** classes.
  They are the charm-sector analogues of the same amplitude classes; the
  algebraic Clebsch-Gordan prefactors (e.g. `(2/3)^(1/2)`, `(1/6)^(1/2)`)
  are identical in form to the light rows — only the class letter carries a
  `_c` subscript.
- The **"Realistic" factor column** for every charmed row is NOT a number.
  It is the bracketed *class ratio* `[A_c/A]` (for A_c rows) or `[S_c/S]`
  (for S_c rows). I.e. the realistic correction is "replace the light class
  value by its charm counterpart." This is structurally different from the
  light section, where the realistic factor is a parenthesized number like
  `(1.5)`. Recorded verbatim in `realistic_factor` as `[A_c/A]` / `[S_c/S]`.
- **β_c** (footnote d, verbatim): "β_c denotes the harmonic-oscillator
  parameter appropriate to charmed mesons which in our numerical results we
  have taken equal to β." So numerically **β_c = β**, but it is kept as a
  distinct symbol in the amplitude formulas.
- Footnote d also gives the **modified form factor** for these charmed
  states:  exp[ −(1/4)(m_c/(m_c+m_d))² q²/β_c² ]
  (note it is written with β_c in the exponent denominator, taken = β).
- The three P-wave amplitudes that are A_c and carry a recoil/unequal-mass
  factor print an explicit multiplier:
  `A_c q̃² · m_c β / ((m_c + m_d) β_c)`  (q̃ ≡ q/β).
  Rows: K*_c→Dπ, K*_c→D*π, Q_1c→[D*π]_D, Q_2c→[D*π]_D. The S-wave S_c rows
  and the 1³S₁ A_c q̃ rows do NOT carry this factor.
- Footnote **d** on charmed P-wave rows: "We denote a charmed P-wave meson
  by the name of its I=1/2 u-d-s analog with a subscript c" — hence the
  names K*_c, Q_1c, Q_2c, κ_c mirror the strange K*, Q_1, Q_2, κ.
- Footnote **b** on the Q_ic rows: "We quote the pure Q_B(A) [pure
  singlet(triplet)] amplitude formulas under Q_1(2), where Q_1(2) is the
  lower(higher) state in mass." Same mixing caveat as the strange Q_1/Q_2.

### Low-confidence cells (flagged, not guessed)
1. **Q_1c → [D*π]_S** (`confidence=low`): coefficient expression is printed
   **+(1/6)^(1/2)** (positive) yet the paper's (MeV¹ᐟ²) value is **−1.5**
   (negative). Sign mismatch between the algebraic coefficient and the
   numeric column. Both taken verbatim from the image; the discrepancy is
   the paper's, flagged for W2c to resolve (likely an S_c sign convention or
   a print sign on the coefficient).
2. **Q_2c → [D*π]_D** (`confidence=low`): coefficient printed **−(1/6)^(1/2)**
   (negative) yet the (MeV¹ᐟ²) value is **+0.7** (positive). Same
   coefficient-vs-numeric sign mismatch as above.

Both sign mismatches are real on the page image (not OCR noise) — the
coefficient sign and the numeric sign genuinely disagree for these two rows.
All other rows have consistent signs between coefficient and numeric value.

### Skipped rows
None. All 10 printed charmed-section decays are present.

### Numeric values verified off page-022.png (MeV¹ᐟ²)
-0.34, +0.24, -0.27 (1³S₁); -7.3, -5.1 (1³P₂);
-1.5, +5.3, +15, +0.7 (1³P₁/1¹P₁); -15 (1³P₀).

## Light 2S and 1D-remainder section (unit W2a-1)

### Source
- OCR: `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`, lines ~721–930
  (`## PDF Page 17`–`## PDF Page 21`, `TABLE V. (Continued)`).
- Every row image-verified against `page-017.png` (printed p.205),
  `page-018.png` (p.206), `page-019.png` (p.207), `page-020.png` (p.208),
  `page-021.png` (p.209).
- Output: `table_v_light_2s_1d.provisional.csv` (**94 rows**).
- The light `1³D₃` rows (page 16) are NOT here — they are in
  `table_v_light_1d3.provisional.csv` (owned by another unit). Not duplicated.

### Row count / structure (per `section` column)
- `1^3D_2 nonstrange` — 18 rows (ρ, ω, φ isoscalars; page 17). ρ/ω/φ→[X]_P
  (D-wave) and →[X]_F (F-wave) pairs; φ→ρπ/ωη P&F are zero (ideal mixing).
- `1^3D_2 1^1D_2 strange` — 16 rows (Q₁, Q₂; pages 17–18). Footnote **c**:
  "we quote the pure ¹D₂(³D₂) amplitude formulas under Q_1(2), where Q_1(2)
  is the lower(higher) state in mass."
- `1^1D_2 nonstrange` — 8 rows (A₃, ω, φ; page 18). Header on the image
  reads **1¹D₂ (nonstrange)** (OCR mislabeled; corrected here).
- `1^3D_1` — 21 rows (ρ_D, ω_D, φ_D, K*_D; pages 18–19). All D-wave, all
  experiment column = "see Sec. VA".
- `2^1S_0` — 9 rows (π′, η_r, η′_r, K′; pages 19–20). All P-wave (class P).
- `2^3S_1` — 22 rows (ρ_S, ω_S, φ_S, K*_S; pages 20–21). All P-wave.

### Amplitude-class conventions as PRINTED
- **D-wave (P-wave partial wave, `_P` label)** rows in the 1D nonets carry
  amplitude class **D** with `q̄¹`, and the "Realistic" factor column prints
  the class ratio **[D/S]** (recorded verbatim in `realistic_factor`).
- **F-wave (`_F` label)** rows carry class **A** with `q̄³`; their realistic
  factor is a parenthesized number, e.g. `(1.7)`, `(1.6)`.
- **2S nonets** (2¹S₀, 2³S₁) use amplitude class **P** with `q̄¹`; the
  realistic factor column prints the class ratio **[P/S]**.
- These class letters D and P are the structure-dependent amplitudes defined
  in Table IV / the ledger; `qbar_power` records the orbital power of q̄.
- The two OCR column-order quirks (P-wave rows put [D/S] in the realistic
  column with no experiment value) were resolved by the image.

### New classes / footnotes encountered
- `[D/S]` and `[P/S]` class-ratio realistic factors (new vs the 1S/1P file,
  which used numeric factors and `mixing_only`).
- Footnote **c** (¹D₂/³D₂ pure-formula caveat for Q₁/Q₂) — new usage here.
- Footnote **f** appears on several 2S/1D rows (final-state-width `(AB)_M P`
  convention), same footnote as the 1S/1P file.

### Low-confidence cells (flagged, not guessed) — 4 total
1. **Q₂→[K*η]_P** (`1^3D_2 1^1D_2 strange`, p.206): coefficient printed
   **+[(√2−1)/√288]** (positive) but the (MeV¹ᐟ²) value is **−1.6**
   (negative). Coefficient-vs-numeric sign mismatch, real on the image.
2. **ω_S→(ππ)_ρπ** (`2^3S_1`, p.208): coefficient printed **+(1/27)^(1/2)**
   (positive) but the (MeV¹ᐟ²) value is **−5.8** (negative). Real sign
   mismatch on the image.
3. **η_r→(Kπ)_K*Kbar** (`2^1S_0`, p.207): the (MeV¹ᐟ²) cell is a two-part
   value **"−0.11 − 1.1"** (mixing-angle-dependent range, not a single
   number). Stored verbatim in `paper_MeV12`.
4. **η′_r→K*Kbar** (`2^1S_0`, p.207): same — two-part value
   **"−0.46 + 2.4"**. Stored verbatim.

### Other image corrections vs OCR (not low-confidence)
- **K*_S→ρK** (`2^3S_1`, p.209): OCR rendered the coefficient as
  `+(−1/108)^(1/2)` (an artifact); the image clearly shows `+(1/108)^(1/2)`,
  value +2.9. Recorded as `+(1/108)^(1/2)`.
- **φ→[K*Kbar]_F** (`1^1D_2 nonstrange`, p.206): coefficient `−(1/400)^(1/2)`
  read directly from the image (denominator 400), value −2.8.

### Skipped rows (with reasons)
- The **1³F₄** block on page 209 (δ, h, h′, K* with `A q̄⁴`, F-wave/G-wave)
  is present on the same page as the 2³S₁ tail but is a **1F** nonet, OUT of
  the 2S/1D scope of this unit. NOT digitized here — belongs to a 1F/2P unit.
- No other rows skipped; every printed 2S and 1D-remainder decay is present.

### Strange 2S/1D rows in Table V? YES
- Strange **1D**: Q₁/Q₂ (`1^3D_2 1^1D_2 strange`, 16 rows) and K*_D
  (`1^3D_1`, part of its 21 rows).
- Strange **2S**: K′ (`2^1S_0`, part of 9 rows) and K*_S (`2^3S_1`, part of
  22 rows). Included in this file.

### Structural surprises
- The `1^1D_2 (nonstrange)` section header (A₃/ω/φ) was mislabeled in the
  vision OCR (it attached "(nonstrange)" to the A₃ row and omitted the header
  between Q₂ and A₃). Corrected from the image.
- The 2¹S₀ / 2³S₁ boundary falls mid-page: the `2^3S_1` header appears on
  page 208 after `K′→φK`, so K′ (ρK/ωK/φK) rows on page 208 are still 2¹S₀.
- η_r / η′_r 2¹S₀ rows carry two-part mixing-dependent numeric values rather
  than single numbers.
