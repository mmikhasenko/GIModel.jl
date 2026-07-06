# Table V (Strong Decays) — Digitization Notes

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
