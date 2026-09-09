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
| **C. Ours, not in the paper** | No active spectrum parameters | The former `k_spin_orbit`/`k_tensor` bridge factors were deleted in PA-18 | The native HO spectrum uses only the paper's Table-II inputs and literal A15-A16 coefficients. Calibrated annihilation and strong-decay controls remain explicitly named, report-side comparisons rather than hidden spectrum strengths. |

Everything else in the model is derived (kinetic operator, potentials,
Clebsch-Gordan/angular factors, overlap kernels) with no adjustable constant.
The photon (Table VI) and general-annihilation (Table III/16) audits add **no
new fitted constants** on our side — they reuse tier-A/B values only.

## 1. Where the residuals live (the logic)

Across the current audits the largest remaining deviations concentrate in two
places:

1. **Same-J mixing angles.** The literal A15-A16 calculation leaves substantial
   convention-aware residuals in several strange, charm, and bottom 1P blocks;
   the current values and complementary-angle comparison are generated in
   `mixing_angles.md`. These propagate directly into the Q1/Q2 and Q1c/Q2c
   strong-decay rows, which remain sensitive to the paper's quoted angles.
2. **Flavor-mixed and cancellation-sensitive observables**, especially radial
   excited η/η′ channels and Υ″→χ_b0γ. These now consume the
   final physical composition. Their visible residuals are model/convention
   discrepancies, not parallel consumer state.

W6 is complete: converged native HO and FD agree on the central smeared origin
functional to 0.44% and on the light-pion full solve to 0.1 MeV. PA-17 and
PA-18 are complete; no additional radial-wave or mesh abstraction is needed.
The optional FD-COMP follow-up is also complete: independent spacing/domain
sweeps certify a precision comparator profile, mixed eigenspaces, and wave-
sensitive observables without making FD part of the paper route.

## 2. Coverage completion

Table VI now computes all 79 canonical rows (42 M1, 35 E1, 2 M2), including
excited bottomonium and heavy-to-light mixing-induced amplitudes. Shared
four-flavor compositions and native fixed-channel waves replace the old
partial audit. The report records signed and magnitude residuals, explicit
kinematic inputs, and the paper-supplied +0.01 μN correction to phi -> pi gamma.

No model-output table remains partially implemented. This is a coverage claim:
the excited eta amplitudes and cancellation-sensitive transitions do not all
agree with the paper. Their residuals remain available in the generated reports.

Context/input tables I and VIII, superseded equations, and applying Eq. (19)
directly to calculated physical waves beyond the paper's single-beta SHO
approximation remain outside the reproduction boundary.
