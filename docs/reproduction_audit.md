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
| **C. Ours, not in the paper — diagnostic bridging scales** | `k_spin_orbit = 0.48`, `k_tensor = 0.42` (`[fine_structure]` in the TOML) | Two global O(1) multipliers inherited from the earlier FD calibration and currently applied in both backend-native fixed-sector matrices | **This is the remaining caveat on literal spectrum reproduction.** Native HO now performs the paper-order fixed-sector diagonalization and needs no mesh bridge, but these two non-paper strengths still scale its spin-orbit and tensor blocks. PA-18 must remove them or isolate them in an explicitly named comparator mode, then certify the formula normalization with paper inputs. |

Everything else in the model is derived (kinetic operator, potentials,
Clebsch-Gordan/angular factors, overlap kernels) with no adjustable constant.
The photon (Table VI) and general-annihilation (Table III/16) audits add **no
new fitted constants** on our side — they reuse tier-A/B values only.

## 1. Where the residuals live (the logic)

Across the current audits the largest remaining deviations concentrate in two
places:

1. **Small-offdiagonal same-J mixing angles.** The model's own K1 (strange 1P
   +19.5° vs paper +34°) and charm 1P (−24.8° vs −41°) angles deviate when the
   singlet-triplet splitting is only a few MeV. This *alone* explains the
   Q1/Q2 and Q1c/Q2c decay misses: those rows reproduce the paper's numbers
   exactly when evaluated at the *paper's* angle. Large-offdiagonal angles
   (1D, 1F blocks) agree to a few degrees.
2. **Flavor-mixed and cancellation-sensitive observables**, especially radial
   η/η′ signs, Υ″→η_bγ, and 2S→χ₀ magnitudes. These still use partially bespoke
   consumer composition and are the PA-17 migration target.

W6 is complete: converged native HO and FD now agree on the central smeared
origin functional to 0.36% and on the light pion full solve to 0.1 MeV. The
remaining work is therefore final-state consumer integration (PA-17) and
removal/isolation of the two bridge strengths (PA-18), not another radial-wave
or mesh abstraction.

## 2. Remaining work and costs

Per-unit status is in the dashboard; this is the workstream-level costing of
what still stands between the current state and "the whole paper, reproduced".

| Remaining | Cost | Blocks victory? |
| --- | --- | --- |
| **PA-17 final-state consumers** | medium | yes for flavor-mixed observable certification |
| **PA-18 literal paper strengths and end-to-end HO certification** | research/calibration | yes for an unqualified original-algorithm claim |
| **FD-COMP independent convergence audit** | medium | no; comparator evidence only |

**Not reproduction blockers** (context/superseded — candidates for the
direction-2 demos instead): Table I, Table VIII, Eqs. (11), (15), (23)-(29),
A1-A6, B37.

**Honest one-line status:** the three-stage native-HO spectrum algorithm,
adaptive basis convergence, and all paper table audit machinery are present;
the remaining qualification is that some flavor-sensitive observables have not
yet migrated to the final physical composition and two non-paper fine-structure
bridge scales remain active.
