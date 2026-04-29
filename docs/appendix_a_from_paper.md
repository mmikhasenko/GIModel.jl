# Appendix A: structure from the 1985 paper (navigation)

This note orients the repository toward the **original** Godfrey–Isgur
implementation of relativistic smearing. For full equation typography and
factors, use `paper/Godfrey-Isgur-1985.pdf` (PDF **pp. 36–38**). The
`paper/text/*` extractions are **OCR-fragile** in (A12)–(A15); do not use them
for coefficients without checking the scan.

## Smearing setup

- **(A7)** – Definition of a **smearing function** in three dimensions.
- **(A8)** – The smeared form of a potential:
  `f̃(r) = ∫ d³r' ρ(|r - r'|) f(r')` (schematic: convolution with a normalized
  kernel; see the PDF for the exact definition and notation for `G` and `S`).
- **(A9)** – Universal **σ** and **s** (Table II). In this repo they feed
  `contact_smearing_sigma` and the (experimental) 3D blur
  `smear_3d_radial` / `smeared_central_values`.

## Pointwise “text” potentials

- **(A10)** and **(A11)** – Relativized strengths for spin-dependent
  interactions, tied to the **ε** factors in Table II (already applied in
  `src/GIModel/` for contact and first-order fine structure where coded).
- The **Coulomb** and **confinement** building blocks in coordinate space
  (orientation around Eqs. (11)–(13) in the main text) match our `static_coulomb_G`
  and `static_confinement_S` (see `docs/formula_map.md`).

## What we do **not** yet implement (main gap)

- **(A12)** – **Expanded** effective **Coulomb** piece after the (A7)–(A8) smearing,
  with **derivative** structure (the HO diagonalization in the paper uses this
  form, not a naive 3D average of the pointwise $-4\alpha_s/(3r)$ on a line).
- **(A13)** – Analogous **expanded** **confinement** $S$; the paper also states
  a pragmatic simplification: **$S$ may be left unmodified** in their QED
  one-dimensional analog while **$G$** is expanded ((A12)–(A14) region in the PDF).
- **(A10)**-style Q-dependent **modifications** to spin–orbit and tensor
  *operators* are only partially represented via constant **ε** in the current
  code.

**Repository consequence:** the flag `appendix_a_smearing` enables only the
**experimental (A7)–(A8) style 3D blur** in `smeared_central_values` (see
`src/GIModel/GIModel.jl`), which `docs/formula_map.md` states is **not** the
same as (A12)–(A13). The flag `appendix_a_derivative_g` enables a modular
finite-difference **proxy** for the first derivative-smearing term,
`G + ∇²G/(4σ²)`, while leaving `S` pointwise. This is useful for comparison and
testing the central-potential interface, but it is not yet the fully audited
(A12)–(A13) expression. The next **serious** coding milestone is a paper-faithful
(A12) path (or a verified HO-basis port), then aligning spin-dependent radial
integrals with the same smeared $G$.

## Quick PDF map

| PDF page (approx.) | Content |
| ---: | --- |
| 36 | (A7)–(A11) |
| 37 | (A12) start, (A9) text |
| 37–38 | (A12)–(A16) expanded potentials |

## Related code

- `smeared_central_values` / `smear_3d_radial` — `src/GIModel/GIModel.jl`
- `appendix_a_derivative_central_values` —
  `src/GIModel/appendix_a_derivative_potential.jl`
- `src/GIModel/appendix_a_status.jl` — explicit status (which path is active)
- `scripts/compare_central_paths.jl` — pointwise, 3D blur, 1D G blur, and
  derivative proxy on one grid
- `scripts/compare_central_pointwise_vs_a7a8.jl` — pointwise vs (A7)–(A8) blur on
  a grid (diagnostic only)
