# Reproduction Audit

The single gate for the claim "the original Godfrey-Isgur model is fully
reproduced". Two inventories: every **table/figure** with its computed
local-vs-paper asset, and every **tagged equation** (main text 1a-29,
A1-A17, B1-B37, C1-C3, D1-D11 — 99 tags total) with its implementation
status. Detailed per-topic ledgers stay in `docs/observable_ledger.md`,
`docs/appendix_a_equation_audit.md`, `docs/formula_map.md`, and
`docs/paper_gap_ledger.md`; this file only aggregates their verdicts.

Status vocabulary:

- **reproduced** — computed locally from the model and compared against the
  paper's numbers in a checked-in asset.
- **partial** — a subset is reproduced; the missing part is itemized.
- **implemented** — the physics is in `src/` (or an audit script) but no
  systematic number-vs-paper asset exists yet.
- **folded** — not coded as a standalone formula; its content enters through
  a convention, an input table, or a more general implementation.
- **context** — derivational or discussion material with nothing executable
  to reproduce.
- **missing** — needed for the victory claim and not present.

## 1. Tables and figures

| Item | Content | Status | Computed asset | Gap to close |
| --- | --- | --- | --- | --- |
| Table I | Motivational: importance of confinement in heavy QQ̄ | context | — | Nothing to reproduce; candidate for a TIL demo. |
| Table II | Model parameters | reproduced (as audited input) | `data/table_ii_parameters.csv` + TOML, sync-checked by `scripts/data_checks.py validate` | — |
| Table III | Isoscalar mixed compositions, n = 1, 2 | reproduced | `GIPaper/docs/residual_reports/table_iii_mixing_audit.md`, `annihilation_model_scorecard.md` | General Eq. (16) blocks match the ω/φ and f₂/f₂′ eigenvectors; literal P1/P2 pseudoscalar modes reproduce the structure at amplitude RMS 0.180/0.114 while the calibrated P1 control carries the mass scoring. Literal-mode mass scoring is the residual item. |
| Table IV | Reduced partial-wave amplitude classes | reproduced | `src/strong_decays.jl` + tests | Light classes (A, A′, A″, A₀, S, D, P) coded and calibrated; charm classes `A_c`, `S_c` now coded (`f2dadf2`) with the footnote-d form factor and the P-wave recoil multiplier, β_c≡β. Key result: the light calibration (A=1.665, S₀=3.918) transfers to charm with **no refit**. |
| Table V | Strong decay amplitudes (multipage) | partial | `GIPaper/docs/residual_reports/table_v_light_decays.md` | Light 1S+1P audited (33 headline rows, 6% median). **Q1/Q2 K1-mixing rows now scored** three ways (`c5d7ccd`): paper column matches its own +34° quote (10% median) not our +19.5° model angle. Light 1³D₃ digitized (`8823643`). **Charmed section digitized + scored** (`44e1d50`, `f2dadf2`): 6 clean rows at 12% median; Q1c/Q2c reproduce at the paper's −41° charm 1P angle (model's own −24.8° misses the cancellation rows — same story as strange K1). **Light 2S + 1D-remainder digitized** (`824eaed`, 94 rows incl. strange 2S/1D, page-verified). Digitized-but-not-yet-scored: the 2S/1D/1³D₃ rows (audit currently scores 1S+1P + charmed). Missing: charmonium (ψ) section, the 1³F₄ nonet (δ/h/h′/K*, found on page 209), scoring the digitized D/2S rows, strange √3 recoil (W2c-2). |
| Table VI | Photon decays: M1 moments, E1 amplitudes | partial | `GIPaper/docs/residual_reports/table_vi_photon_decays.md` | 28 mixing-free rows reproduced (quarkonium M1 0.1-2%, open-flavor M1 1-4%, light ~6%, E1 3-7%). **Open-flavor M1 column now image-verified against printed page 24** (`07ad617`): the one-row OCR shift was correct, all 8 values confirmed, caveat retired. **2S → χ₀ is not a q-convention artifact**: model-mass q helps Υ′→χ_b0 (30%→16%) but hurts ψ′→χ_c0 (27%→43%) — the c-c̄ 2S-1P gap is genuinely compressed. **6 isoscalar M1 rows added** (`4a8645c`) via the audited Table III mixing (P1 pseudoscalar + S1 vector blocks, no new constants): all signs reproduced including the nontrivial φ→η(+)/φ→η′(−) relative sign; dominant-flavor magnitudes at the ~15-24% light-sector I-overlap residual; sub-dominant rows ~2-3× low because η↔η′ ordering is inverted vs the paper (needs the paper's per-row I mock-mass evaluation). **All printed blocks now audited** (`d678136`): light E1/M2 (A₂→πγ fit row +0.510 vs +0.55 ~7%; A₂→ργ/ωγ within 1%; A₁/B ~8-10%; f′→φγ sign+mag), strange K*(1420)→Kγ E1 +0.436 vs +0.48 (~9%; no separate charmed P-wave E1 rows exist), hindered bottomonium M1 (Υ′→η_bγ +0.009 vs +0.007; **Υ″→η_bγ sign-flips** −0.004 vs +0.007 — deepest 3S→1S/E2 cancellation, points at W6). Open items: η↔η′ ordering and the Υ″ sign both trace to the central-solve vs paper-HO-order radial residual (W6). ψ/Υ→(light)γ order-of-magnitude rows (footnote d) deliberately not modeled. |
| Table VII | Leptonic, γγ, gluonic decays; charge radii | missing | — | Not started. Needs D4-D8 (f_P, f_V definitions, leptonic widths), γγ and gluonic width formulas, and ⟨r²⟩ from the model wavefunctions. |
| Table VIII | Input masses/widths for the 0⁺⁺ coupled-channel discussion (Fig. 16) | context | — | Discussion-layer phenomenology (Sec. V D), not a model output. Optional. |
| Figs. 3-9 | Full meson spectra, all flavor sectors | reproduced | `GIPaper/docs/residual_reports/scorecard.md` + 7 sector reports | Mean abs residuals: bottomonium 4.3, b-flavored 5.2 (bū/bd̄ + bs̄ + bc̄), charmonium 6.0, strange 10.8, isoscalar 11.5, charmed 12.4, isovector 14.5 MeV. Regenerated byte-identically by `run_all_spectrum_checks.jl`. |
| Fig. 2 | Saturating α_s | reproduced | regression tests on `alpha_s_r` and derivatives | — |
| Mixing angles (text, Secs. V) | θ_nL quoted near Figs. 5-9: strange θ₁P≈34°, θ₁D≈33°, θ₂P≈15°, θ₁F≈32°, θ₂D≈25°, θ₁G≈33°; charm θ₁P^cū≈−41°, θ₁D^cū≈−39°, θ₁P^cs̄≈−44°, θ₁D^cs̄≈−39°; bottom θ₁P^bū≈−43°, θ₁P^bs̄≈−45°, θ₁P^bc̄≈−53° | partial | `GIPaper/docs/residual_reports/mixing_angles.md` (all 13 angles, convention derived once) | Large-offdiagonal blocks reproduced (charm/strange 1D at 0.3-2.5°, strange 1F 4.7°); b-sector 1P content agrees to 1-5° once the low/high labeling of the nearly degenerate eigenstates is exchanged. Open: small-offdiagonal 1P blocks (K1 +19.5° vs +34°, cs̄ −14.1° vs −44°) inherit few-MeV singlet-triplet splitting sensitivity, and the strange radial trend is inverted (paper decreases with n, ours increases) — resolution likely sits with W6 (HO-order spin blocks), not with the angle machinery. |

## 2. Equation inventory

### Main text (1a-29)

| Eq. | Content | Status | Where |
| --- | --- | --- | --- |
| 1a | Rest-frame Schrödinger equation H\|Ψ⟩ = E\|Ψ⟩ | implemented | `channel_solution` / `hamiltonian.jl` (FD), HO basis path |
| 1b | H₀ = √(p²+m₁²) + √(p²+m₂²) | implemented | `sqrt_kinetic_matrix_from_eigen` on the FD p² eigenbasis; regression-tested |
| 2a | Nonrelativistic kinetic limit | implemented (diagnostic) | `kinetic = :nonrelativistic` option |
| 2b | V = H^conf + H^hyp + H^so + H_A decomposition | implemented | the staged `Spectrum` pipeline + annihilation blocks mirror exactly this split |
| 3 | H^conf: linear + Coulomb with color factor | implemented | `central_potential`; ⟨F·F⟩ = −4/3 folded in |
| 4 | H^hyp: contact + tensor | implemented | contact: smeared, sandwiched, nonperturbative S-wave; tensor: diagonal shifts + same-J off-diagonal blocks |
| 5 | H^so = so(cm) + so(tp) split | implemented | `fine_structure_components` returns both pieces separately |
| 6 | Color-magnetic spin-orbit | implemented | pointwise α_s/r³ and smeared (1/r)dG̃/dr branches; symmetric L·S + antisymmetric same-J block |
| 7 | Thomas-precession spin-orbit | implemented | (1/2r) dH^conf/dr with smeared G̃+S̃ derivatives in the active branch |
| 8 | F_i color matrices (quark/antiquark) | folded | enters only through ⟨F·F⟩ |
| 9 | ⟨F_i·F_j⟩ = −4/3 (meson) | folded | sign/normalization audited in `formula_map.md` |
| 10 | Gaussian smearing function ρ_ij | implemented | `delta_sigma_3d` (≡ A7 kernel); normalization tested |
| 11 | Perturbative α_s(Q²) | context | only motivates the Fig. 2 fit; nothing to code |
| 12 | Saturating α_s(Q²) = Σ a_k exp(−Q²/4γ_k²) | folded | Fig. 2 caption coefficients enter through Eq. (13) |
| 13 | α_s(r) = Σ a_k erf(γ_k r) | implemented | `alpha_s_r` + analytic derivatives; regression-tested |
| 14 | Staged diagonalization of H̃₁ then perturbative mixing matrices | implemented | this is precisely the `central_spectrum → add_spin_corrections → add_intra_meson_mixing` pipeline; FD is the headline basis, HO central comparison complete (sub-MeV); HO-order validation of the spin/mixing blocks is the residual item |
| 15 | Annihilation order-of-magnitude α_sⁿ\|Ψ(0)\|²/M² | context | motivates Eq. (16) |
| 16 | General annihilation matrix element | implemented | `isoscalar_general_annihilation_solution`: 4π(2L+1), (α_iα_j/π²)^(n/2), S_L factors, 1/(m_im_j) all audited |
| 17 | S_L(Ψ) wavefunction factor | implemented | j_L momentum transform on phase-fixed HO waves |
| 18a | Pseudoscalar P1 bracket | implemented | `PaperP1Annihilation` |
| 18b | Pseudoscalar P2 bracket | implemented | `PaperP2Annihilation` |
| 19 | Pseudoscalar-emission amplitude (g σ·q ± h σ·p′) | partial | evaluated in the SU(6)/single-β SHO limit via the Table IV classes (`strong_decays.jl`); the operator itself is never applied to model wavefunctions |
| 20 | Realistic-factor ratio ⟨³S₁\|r^{L−1}\|M*⟩ / ⟨¹S₀\|…⟩ | missing | the parenthetical correction column of Table V; needs radial moments on model waves |
| 21 | Realistic-factor ratio for S-type (p matrix element) | missing | as Eq. (20) |
| 22 | Photon-emission helicity amplitude | implemented | Appendix-D mock-meson form in `audit_table_vi_photon_decays.jl` (to be promoted into `src/`) |
| 23-26 | Perturbative P-wave multiplet formulas E(J^PC) = E₀ + aS + bT + cL | folded | the full diagonalization supersedes them; the S/T/L decomposition is exposed by `fine_structure_components`, so a formula-level demo is cheap (TIL candidate) |
| 27-29 | 0⁺⁺ coupled-channel S-matrix / Breit-Wigner unitarity | context | Sec. V D data phenomenology around Table VIII / Fig. 16, outside the model proper |

### Appendix A (A1-A17) — detailed ledger in `docs/appendix_a_equation_audit.md`

| Eq. | Content | Status |
| --- | --- | --- |
| A1-A4 | Scattering-amplitude setup motivating the effective potential | context |
| A5-A6 | Vector/scalar effective kernels | context (deliberately not a coding source; OCR spin labels unsafe) |
| A7 | Gaussian smearing kernel | implemented (contact + diagnostic 3D smear; normalization tested) |
| A8 | Smeared potential definition | implemented (diagnostic convolution; active path uses closed forms) |
| A9 | Mass-dependent width σ(m₁,m₂) | implemented (`contact_smearing_sigma`, shared by all smeared kernels) |
| A10-A11 | Pointwise G(r), S(r) | implemented |
| A12-A14 | Closed-form G̃, S̃, τ_k | implemented — the **active** central path; tested against quadrature and derivatives |
| post-A14 | Momentum factors: central A(p)G̃A(p); spin-side (m/E)^(1/2+ε) sandwiches | implemented (`appendix_a_momentum_sandwich_matrix`; contact/fine-structure sandwiches) |
| A15 | Effective Coulomb-side spin operators | implemented (smeared-kernel branch); paper-order HO validation pending |
| A16 | Scalar/Thomas spin-orbit operator | partial (equal-mass exact; unequal-mass antisymmetric piece handled at the mixing stage) |
| A17 | HO matrix-element factorization | partial (central FD/HO comparison complete at sub-MeV; spin/mixing HO-order comparison pending) |

### Appendix B (B1-B37)

| Eq. | Content | Status |
| --- | --- | --- |
| B1-B9 | Flavor wavefunctions (π, K, η₈, η₁) with sign conventions | folded (flavor labels and the ns̄/ss̄ basis conventions in `meson.jl` / annihilation basis labels) |
| B10-B13 | Ideal-mixing basis M_ns, M_s and mixing rotation | implemented (ideal mixing is the default isoscalar scheme; the rotation is the Table III eigenvector convention) |
| B14-B15 | Perfect-mixing η, η′ | implemented (used by the Table V ¹S₀ formula column audit) |
| B16-B25 | Flavor operators X_q^i (and antiquark rule B25) | folded — their matrix elements are baked into the digitized Table V coefficient column, not derived independently; independent derivation would upgrade the Table V audit from "coefficients transcribed" to "coefficients derived" |
| B26-B30 | Spin wavefunctions χ | folded (standard angular algebra: `spin_dot`, `LdotS`, tensor factors — tested via sum rules) |
| B31-B36 | SHO wavefunctions Ψ_nLM | implemented (`harmonic_oscillator_basis.jl`; also implicit in the Table IV class formulas) |
| B37 | Example composed state vector | context (convention demo) |
| — | ⁴Note: B-appendix β powers appear inside Table IV classes | folded |

### Appendix C (C1-C3)

| Eq. | Content | Status |
| --- | --- | --- |
| C1 | Helicity amplitude definition H_m | folded (Table IV already lists partial waves; we work directly in the partial-wave basis) |
| C2 | Width from partial-wave amplitudes, (q/2π) factor | implemented (the `sqrt(q/2π)` amplitude convention in `strong_decay_amplitude`) |
| C3 | Helicity → partial-wave sum rule | context (not needed while amplitudes are compared, not widths) |

### Appendix D (D1-D11)

| Eq. | Content | Status |
| --- | --- | --- |
| D1 | μ_πω definition from the current matrix element | folded (the μ/(e/2) convention of the Table VI audit) |
| D2 | Mock-meson state definition | implemented (mock mass M̃ = ⟨E₁⟩+⟨E₂⟩ in the Table VI audit) |
| D3 | Relativized moment with (m+2E)/3E² → (m/E)^f prescription | implemented (I_i overlap with (m/E)^0.7, E_n moments with 0.5 — the paper's fitted exponents, no new constants) |
| D4-D6 | f_P, f_V, f_{A₁} decay-constant definitions | missing (needed for Table VII) |
| D7 | Γ(P → lν) | missing (Table VII) |
| D8 | Γ(V → l⁺l⁻) | missing (Table VII) |
| D9 | Γ(τ → A₁ν) | missing (Table VII) |
| D10-D11 | Γ(V→Pγ), Γ(P→Vγ) from μ | folded (trivial conversions; Table VI audit compares moments/amplitudes directly — apply when width columns are scored) |

## 3. Victory checklist

Ordered by what blocks the claim "the original model is fully reproduced":

1. **Mixing-angle audit asset** (small): one computed-vs-paper table for the
   13 quoted θ values (strange 6, charm 4, bottom 3), with the paper's
   rotation convention pinned once. All machinery exists.
2. **Table V completion** (largest): digitize and audit the remaining
   sections (light 2S/1D, strange with recoil factors, charmed, charmonium);
   wire the model K1 angle into the Q1/Q2 rows; code the charm classes
   `A_c`, `S_c`; resolve the strange-parent √3 normalization.
3. **Table VI completion** (medium): isoscalar mixing rows (mixing layer
   already provides amplitudes), remaining E1/M2 blocks, PDF crop audit of
   the OCR-shifted column, promote the Appendix-D kernels into `src/`.
4. **Table VII** (medium): D4-D8 decay constants and widths, γγ/gluonic
   widths, charge radii — all on existing wavefunctions.
5. **Eqs. (20)-(21)** (small): realistic-factor ratios on model waves —
   also directly improves flagged Table V rows.
6. **HO-order validation of spin/mixing blocks** (background): the last
   Appendix-A fidelity item (A15-A17); central operator already validated.

Items that are **not** reproduction blockers: Table I, Table VIII,
Eqs. (11), (15), (23)-(29), A1-A6, B37 (context/superseded); they are
candidates for the direction-2 demos instead.
