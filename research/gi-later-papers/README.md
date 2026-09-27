# What the GI model says about J=L mixing, 1985 → 2016

**Status: closed as an implementation investigation.** The historical cause of
GI85's caption values remains unknown. The [GIPaper benchmark](../../GIPaper/data/mixing_angle_targets.md)
now preserves the original targets alongside the six GK91 values; small
quantitative residuals remain visible. The separate ψ(3.82) tensor question
below remains open.

**Question.** GIModel, after the L>0 contact fix, does not reproduce the J=L
singlet–triplet mixing angles printed in the Godfrey–Isgur (1985) figure
captions. Is that a defect of our calculation, or of the 1985 numbers? What did
the authors, and later reimplementers of the same model, report over time?

**Short answer.** The first follow-up by the model's first author already gives
the small angles, in a draft dated July 1986, a year after GI85. That is
Godfrey & Kokoski (preprint TRI-PP-86-51, revised as GIPP 90-6A, published as
PRD 43, 1679 in 1991). It presents its table as the calculation of GI85,
"Some of these mesons where [sic] omitted in ref. 10 for brevity". Its angles agree
with GIModel in sign and magnitude for all six 1P systems, and its contact,
tensor and spin–orbit expectation values agree with ours to about 1 MeV. It
states that the model's spin–orbit mixing amplitude "is so small" that the
physical mixing may come from decay-channel couplings instead. Every later
GI-model calculation we found, by Godfrey and by an independent group, gives
the same small angles. No publication, including this one, remarks that the
GI85 captions differ. The contact-term fix to GIModel (commit 4c84a9e) is
confirmed by the same paper: its Table IV lists P-wave ⟨H_cont⟩ = 33, 20, 15,
8, 7, 4 MeV for K, D, Ds, B, Bs, Bc.

Numerical evidence on our side:
[GIPaper/docs/investigations/mixing_composition_investigation.md](../../GIPaper/docs/investigations/mixing_composition_investigation.md).

## Timeline

| Year | Work | Model / parameters | What it reports on J=L mixing |
|---|---|---|---|
| 1985 | Godfrey & Isgur, PRD 32, 189 (*GI85*) | Defines the model, Table II parameters | Figure captions: K θ₁P, θ₁D, θ₂P, θ₁F, θ₂D, θ₁G = 34°, 33°, 15°, 32°, 25°, 33°. D: −41°, −39°. Ds: −44°, −39°. B, Bs, Bc: −43°, −45°, −53°. D/Ds J=1 splittings 50/40 MeV. |
| 1986→1991 | Godfrey & Kokoski, PRD 43, 1679 (*GK91*). Earlier version TRI-PP-86-51, July 1986; revised preprint GIPP 90-6A, July 1990 ([scan](Godfrey–Kokoski.pdf)) | GI85 model and parameters (b = 0.18 GeV², m_u = 0.22, m_s = 0.419, m_c = 1.628, m_b = 4.977 GeV). Same convention as GI85 Eq. (10), Q_low = cos θ ¹P₁ + sin θ ³P₁, with heavy quark first (K = sū) | Table I: θ = −5° (K), −26° (D), −38° (Ds), −31° (B), −40° (Bs), +68° (Bc). J=1 pairs K 1.35/1.37, D 2.46/2.47, Ds 2.55/2.56, B, Bs degenerate, Bc 6.74/6.75 GeV. ³P₀, ³P₂ as in GI85. Table IV: ⟨H_cont⟩, ⟨H⁺_so⟩, ⟨H_ten⟩ = 33/47/56 (K), 20/26/27 (D), 15/27/28 (Ds) MeV. Sec. IV: the spin–orbit mixing amplitude is small, so mixing may be due to "couplings to common decay channels" (refs: Lipkin 1977; "S. Godfrey and N. Isgur, in progress"). |
| 1991/92 | Lu, Wise & Isgur, PRD 45, 1553, CEBAF-TH-91-16 ([scan](Lu-Wise-Isgur.pdf)) (*LWI92*) | Heavy-quark symmetry, not a GI-model calculation | No model mixing angle quoted. Notes that since quark-model S-wave widths are ~10× the D-wave ones, "even a weak mixing" of the j = 3/2, 1/2 D₁ states would matter. Consistent with GK91's small-mixing picture; says nothing about GI85's captions. |
| 1996 | Blundell, Godfrey & Phelps, PRD 53, 3712, [hep-ph/9510245](https://arxiv.org/abs/hep-ph/9510245) (*BGP96*) | b = 0.18 GeV², m_u = 0.22, m_s = 0.419 GeV | Unmixed K₁ masses 1.37 (³P₁) and 1.35 GeV (¹P₁). "In this model spin-orbit mixing results in θ_K = −5°" [GK91], with the K₁ masses unchanged at the quoted precision. Called inconsistent with experiment (θ_K ≈ 45° from decays), not with GI85. Also quotes model ⟨H_cont⟩ = 33 MeV for the K 1P. |
| 2004 | Godfrey, PRD 70, 054017, [hep-ph/0406228](https://arxiv.org/abs/hep-ph/0406228) | m_c = 1.628, m_b = 4.977 GeV | Bc 1P: 6706, 6741, 6750, 6768 MeV. θ₁P = 22.4° for the upper state, i.e. the lower state is ≈ 15 % singlet (GI85: 36 %). |
| 2005 | Godfrey, PRD 72, 054029, [hep-ph/0508078](https://arxiv.org/abs/hep-ph/0508078) | "GI [GK91, GI85]" | D 1P: 2399, 2460/2470, 2502 MeV. Radiative widths "calculated using θ = −26°". |
| 2006 | Strange axial-vector mixing review, [hep-ph/0606297](https://arxiv.org/abs/hep-ph/0606297) | — | Lists "θ_K ∼ 34° [GI85], θ_K ∼ 5° [GK91]", both as relativized-quark-model values, without comment. |
| 2014 | Lü & Li, [arXiv:1407.3092](https://arxiv.org/abs/1407.3092) — *independent reimplementation* | "parameters taken from" GI85; Gaussian expansion basis, same two-stage procedure | D 1P 2398, 2455, 2467, 2501 MeV. θ₁P = −25.5°, θ₁D = −38.2°, θ₁F = −39.5°, θ₁G = −40.2°. States that the GI85 1S–1F masses are "well reproduced"; no comment on the angles. |
| 2016 | Godfrey & Moats, PRD 93, 034035, [arXiv:1510.08305](https://arxiv.org/abs/1510.08305) | b = 0.18, m_q = 0.22, m_s = 0.419, m_c = 1.628 | D: θ₁P = −25.68°, θ₁D = −38.17°, θ₁F = −39.52°, θ₁G = −40.18°. Ds: −37.48°, −38.47°, … Full mass tables. |
| 2016 | Godfrey, Moats & Swanson, PRD 94, 054025, [arXiv:1607.02169](https://arxiv.org/abs/1607.02169) | m_b = 4.977 GeV; GI columns | B: \|θ₁P\| = 30.28°, Bs 39.12°, with J=1 pairs 5777/5784 and 5857/5861 MeV. |

Later GI-based work mostly cites GM16/GMS16/G04 for angles, e.g.
[arXiv:2411.02976](https://arxiv.org/abs/2411.02976), which compares with
θ(1P) = −25.68°, θ(1D) = −38.17°. The "modified GI" (screened) strange-meson
study [arXiv:1705.03144](https://arxiv.org/abs/1705.03144) fits θ₁P to K₁ decays
instead of predicting it.

## Comparison with GIModel

Lower-state singlet probability cos²θ, which does not depend on sign or
ordering conventions. Signs differ between papers: GM16 and GMS16 use opposite
signs for D and B with the same quark ordering, so only magnitudes are compared.

| State | GI85 caption | GK91 (1986/91) | later GI-model result | GIModel now |
|---|---:|---:|---:|---:|
| K 1P | 0.687 | 0.992 | ≈ 0.99 (BGP96) | 0.994 |
| D 1P | 0.570 | 0.808 | 0.812 (GM16), 0.814 (Lü–Li) | 0.811 |
| Ds 1P | 0.517 | 0.621 | 0.630 (GM16) | 0.594 |
| D 1D | 0.604 | – | 0.618 (GM16, Lü–Li) | 0.619 |
| Ds 1D | 0.604 | – | 0.614 (GM16) | 0.613 |
| B 1P | 0.535 | 0.735 | 0.746 (GMS16) | 0.781 |
| Bs 1P | 0.500 | 0.587 | 0.602 (GMS16) | 0.567 |
| Bc 1P | 0.362 | 0.140 | 0.146 (G04) | 0.127 |

Signed angles in GI85's convention and quark ordering are directly comparable
for GI85, GK91 and GIModel. GK91's kaon is sū, which flips its sign.

| | K | D | Ds | B | Bs | Bc |
|---|---:|---:|---:|---:|---:|---:|
| GI85 caption | +34° | −41° | −44° | −43° | −45° | −53° |
| GK91 Table I (K sign flipped to us̄) | +5° | −26° | −38° | −31° | −40° | +68° |
| GIModel | +4.3° | −25.8° | −39.6° | −27.9° | −41.1° | +69.2° |

GK91 and GIModel agree in sign for all six, including Bc, whose GI85 caption
angle has the opposite sign. The spin-dependent expectation values also match:
GK91 Table IV gives ⟨H_cont⟩, ⟨H⁺_so⟩, ⟨H_ten⟩ = 33, 47, 56 (K); 20, 26, 27 (D);
15, 27, 28 (Ds) MeV. Our inverted S, L, 2T are 33.9, 48.0, 56.0; 19.8, 27.5,
27.4; 15.4, 28.3, 28.4.

Masses: GIModel's 1P multiplets (³P₀, J=1 pair, ³P₂) for D, Ds, B, Bs and Bc
agree with GM16/GMS16/G04 within 0–4 MeV in all 20 states. Ours are
systematically 2–4 MeV low in ³P₀. The unmixed K₁ masses agree with BGP96 at
their 10 MeV precision.

## Why the later values are the ones to trust

1. **The earliest follow-up already disagrees with the captions.** GK91 is the
   same model, parameters, author and convention. It is dated a year after
   GI85 and presented as the same calculation. It reproduces GI85's ³P₀, ³P₂
   and fine-structure strengths but has the small V: the J=1 pairs are split by
   10–20 MeV, not 40–50. Our S, T, L, V pattern is the same.
2. **Internal consistency.** Inverting GI85's own printed multiplets
   (GI Eqs. 23–26) gives contact S, tensor T and spin–orbit L that GIModel
   reproduces within 1–4 MeV in all 11 P-wave multiplets. Only the antisymmetric
   element V disagrees: 4–18× larger in the caption for radial-ground-state
   K/D/Ds. No rescaling or sign change of the A15/A16 spin–orbit kernels yields
   GI85's V while keeping GI85's S, T and L. So the caption angles do not follow
   from the Hamiltonian that produced GI85's other numbers, as printed.
3. **Author's own practice.** From GK91 on, Godfrey uses the small angles with
   the same parameters and cites GI85 as the model.
4. **Independent reproduction.** Lü & Li (2014) independently obtain
   θ₁P(D) = −25.5° with GI85's parameters and a different numerical basis.

## What we could not establish

- *Why* GI85 printed larger angles. No publication we found comments on it:
  GK91 presents its Table I as GI85's calculation without noting the
  difference, and LWI92, BGP96, G04, G05, GM16 and GMS16 do not either.
  Possible causes (a pre-publication code version or treatment of the
  antisymmetric term, or caption transcription) cannot be distinguished from
  the literature. Only the authors' notes could settle it.
- Whether the "Godfrey and Isgur, in progress" work on decay-channel-induced
  mixing (GK91 ref. 20) was ever published.
- The ψ(3.82) S–D admixtures (GI85 Fig. 6): no later GI-model amplitude found.

## Compact paragraph (for the report)

> The J=L singlet–triplet mixing angles printed in the figure captions of
> Godfrey and Isgur (1985) are not reproduced by later calculations in the same
> model with the same parameters. The first draft of Godfrey and Kokoski's
> P-wave study (July 1986; published 1991) uses the 1985 convention and is
> presented as the 1985 calculation. It lists θ = −5° (K, as sū), −26° (D),
> −38° (Ds), −31° (B), −40° (Bs) and +68° (Bc), against the captions' 34° (us̄),
> −41°, −44°, −43°, −45° and −53°. Its contact, tensor and spin–orbit
> expectation values are unchanged, and it describes the model's spin–orbit
> mixing amplitude as small. All of Godfrey's later calculations and an
> independent reimplementation (Lü and Li, 2014) agree with these values, and
> none of these papers comments on the difference. GIModel includes the
> smeared contact term in every orbital sector, as Eqs. (A15) and (23)–(26)
> require. It reproduces the 1986–2016 angles, including the sign of the Bc
> angle, to within 0.1–3.5°, and the published 1P masses to within 4 MeV. It
> also reproduces the contact, tensor and spin–orbit strengths implied by the
> 1985 multiplets. Only the 1985 antisymmetric spin–orbit element is
> inconsistent with the rest of the 1985 spectrum. We therefore compare mixing
> with the GI-model values published from 1986 on, and record the 1985 caption
> angles as a historical discrepancy whose origin the publications do not
> explain.

## Files in this folder

- `Godfrey–Kokoski.pdf` and `Godfrey–Kokoski.ocr.txt`: preprint GIPP 90-6A (July
  1990) of PRD 43, 1679, with the text layer from `pdftotext -layout`. The OCR is
  two-up and noisy; check tables against the PDF (Table I on PDF page 9,
  Table IV on page 11, Sec. IV discussion on page 6).
- `Lu-Wise-Isgur.pdf` and `Lu-Wise-Isgur.ocr.txt`: preprint CEBAF-TH-91-16 of
  PRD 45, 1553, with its text layer.
- `references.bib`: BibTeX for every reference above, with the same keys as
  `report/references.bib`.

The report paragraph built from this note is in `report/main.tex`, Sec. on
spectrum reproduction, after the mixing-angle convergence paragraph.

## Reproduce

All arXiv texts were fetched with `curl https://arxiv.org/pdf/<id>` and
`pdftotext -layout`. They are not stored here. GIModel numbers come from
`GIPaper/scripts/investigate_mixing_composition.jl` and
`GIPaper/scripts/investigate_multiplet_parameters.jl`.
