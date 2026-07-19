# Reproduction Audit

> **Retired into the manifest.** The per-unit inventory that used to live here —
> every table, figure, and tagged equation with its status, code, tests, report,
> and page image — is now the machine-checked source of truth in
> `docs/paper_manifest/*.toml`, rendered as the drill-down dashboard
> (`docs/paper_dashboard.qmd`; `cd docs && make dashboard`) and validated by
> `GIPaper/scripts/check_manifest.jl`. Go there for "is unit X done, and where is
> the code?".
>
> What stays in this file is the **cross-cutting synthesis** the per-unit manifest
> can't hold: the provenance-tier honesty ledger for the fitted constants (§0),
> where the residuals concentrate and why (§1), and the remaining-work costing
> (§2). The granular per-row numbers behind every "reproduced" verdict live in the
> residual reports under `GIPaper/docs/residual_reports/`, which the manifest units
> link directly.

## 0. Global scales and fitted constants (honesty ledger)

"Reproduced" means different things depending on where a number comes from.
Every non-derived constant in the pipeline falls into one of three tiers; the
victory claim must be read against this table, not around it.

| Tier | Constants | Source | Bearing on the claim |
| --- | --- | --- | --- |
| **A. Paper's own inputs** (reproduced as inputs, not refit) | Table II: `b=0.18 GeV²`, `c=−0.253 GeV`, `σ₀=1.80 GeV`, `s=1.55`, quark masses, the four `ε` factors; Fig. 2 `α_s` coefficients `(0.25,0.15,0.20)`/denominators `(1,10,1000)`; Table III annihilation amplitudes `A(³S₁)=2.5`, `A(³P₂)=−0.8`, `P1/P2` `A_np=0.50/0.55`; photon smearing exponents `0.7`/`0.5`; decay oscillator `β=0.40 GeV` | Godfrey-Isgur 1985, audited into `data/` and sync-checked | These are the paper's; using them is faithful reproduction, not fitting. |
| **B. Paper's own fit, re-performed identically** | Strong-decay strengths `A=1.665`, `S₀=3.287` (leading convention) | Two-point fit to the *same* rows the paper fits (`ρ→ππ`, `B→[ωπ]_S`) | Methodologically identical to the paper (which gets `A=1.67`, `S=3.27`); no extra freedom introduced. (The earlier `table_iv`-convention value `S₀=3.918` is superseded by the leading-S₀ convention that reproduces Table V — see the `fit-A-S0` manifest unit.) |
| **C. Ours, not in the paper — diagnostic bridging scales** | `k_spin_orbit = 0.48`, `k_tensor = 0.42` (`[fine_structure]` in the TOML) | Two global O(1) multipliers on the **first-order fine-structure expectations evaluated on the FD radial mesh**, chosen to align the FD splittings with the paper's HO-basis result | **This is the one genuine caveat on the spectrum reproduction.** The paper diagonalizes in a harmonic-oscillator basis (Eq. 14) and needs no such factor; our FD first-order fine structure does. The mass *centres*, radial/orbital spacings, contact hyperfine, and mixing *structure* carry no such scale — only the spin-orbit and tensor splitting magnitudes do. **W6** (HO-order validation, A15-A17) is exactly the work that would remove tier-C scales; until then, "spectra reproduced" means centres+spacings to few-MeV and fine-structure splittings up to these two global scales. |

Everything else in the model is derived (kinetic operator, potentials,
Clebsch-Gordan/angular factors, overlap kernels) with no adjustable constant.
The photon (Table VI) and general-annihilation (Table III/16) audits add **no
new fitted constants** on our side — they reuse tier-A/B values only.

## 1. Where the residuals live (the logic)

Across every audited table the deviations concentrate in exactly **two**
places, and neither is a coding error:

1. **Small-offdiagonal same-J mixing angles.** The model's own K1 (strange 1P
   +19.5° vs paper +34°) and charm 1P (−24.8° vs −41°) angles deviate when the
   singlet-triplet splitting is only a few MeV. This *alone* explains the
   Q1/Q2 and Q1c/Q2c decay misses: those rows reproduce the paper's numbers
   exactly when evaluated at the *paper's* angle. Large-offdiagonal angles
   (1D, 1F blocks) agree to a few degrees.
2. **Light-sector wavefunction overlaps (~15-24%)** and the deepest
   multi-node cancellations (η↔η′ ordering, Υ″→η_bγ sign, 2S→χ₀ magnitude).
   All trace to the **FD central-solve vs paper HO-order radial residual** —
   the same root as the tier-C scales in §0.

Both buckets point at the same fix: **W6**, paper-order perturbation theory in
the HO basis. The model reproduces the paper's *algebra, signs, centres, and
spacings*; the open gaps are the physics limits the paper itself discusses,
modulated by our FD-vs-HO basis choice.

## 2. Remaining work and costs

Per-unit status is in the dashboard; this is the workstream-level costing of
what still stands between the current state and "the whole paper, reproduced".

| Remaining | Cost | Blocks victory? |
| --- | --- | --- |
| **Table VII** (D4-D8 decay constants + leptonic/γγ/gluonic widths + charge radii; digitize + audit) | **large** (2-3 units: 1 data, 1-2 code) | yes — last untouched table |
| **W6 HO-order validation** (A15-A17; removes the tier-C scales and the §1 residuals) | **large / research** (background) | it is the deepest fidelity item; "reproduced up to global scales" holds without it, "fully reproduced" is stronger with it |
| **Strange √3 recoil** (W2c-2: derive the Appendix-B isospin factor for `K*₂→Kπ`; now precisely characterized as non-kinematic) | small (1 unit) | soft — one characterized outlier |
| **Eq. (21) type-S momentum ratio** (companion to the reproduced type-A factors) | small (1 unit) | soft — quality, not coverage |
| Promote Table VI/Appendix-D kernels from audit scripts into `src/` + tests | med | no — engineering hygiene |

**Not reproduction blockers** (context/superseded — candidates for the
direction-2 demos instead): Table I, Table VIII, Eqs. (11), (15), (23)-(29),
A1-A6, B37.

**Honest one-line status:** spectra and decay/EM tables are reproduced to the
few-MeV / ~10-20% level the paper works at, with two global fine-structure
scales (§0 tier C) standing in for the paper's HO-order treatment; **Table V is
digitized across every section and reproduced end-to-end (160/178 scoreable rows)
under the leading-S₀ convention**, so the only untouched table is **Table VII**,
and the remaining *fidelity* gap is W6.
