# Charmed-meson counting scratch pad

Generated and checked by `python3 GIPaper/scripts/audit_charmed_census.py`.

The follow-up [phase-space selection](charmed_thresholds.md) applies calculated charm masses and physical emitted-meson masses; its JSON ledger retains open and closed decisions for every entry.

## Verdict

**504 = 18 partial-wave templates × 28 charge/flavor routes.**
**448 = 16 spectroscopic parent/daughter templates × 28 charge/flavor routes.**
The difference is **56 = 2 × 28**: the two axial-basis -> D* templates each have both an S wave and a D wave. Those are two matrix elements of one final channel, not two decays.
An independently written list of 28 routes and 16 templates reproduces every one of the original 504 keys exactly; there are no duplicated keys. Charge conservation and preservation of the heavy spectator are checked separately.

**Important correction to the interpretation:** these are counts in the unmixed spectroscopic basis, not a completed inventory of physical D1/D1′ channels. The prior description of 448 as physical final channels was too strong. Mixing expands the isolated A-double-prime daughter into both physical axial states (see below). No numerical widths were computed in this census.

## Names and finite scope

| Spectroscopic state | Scratch-pad name |
|---|---|
| 1¹S₀ | D |
| 1³S₁ | D* |
| 1³P₀ | D₀* |
| 1³P₂ | D₂* |
| 1³P₁ | D₁[t], a triplet **basis** state |
| 1¹P₁ | D₁′[s], a singlet **basis** state |
| 1³D₁, 1³D₃ | D(1³D₁), D(1³D₃) |
| 2¹S₀, 2³S₁ | D(2¹S₀), D(2³S₁) |

The primes above are bookkeeping labels, not assignments of measured resonances to pure spins. Physical D₁ and D₁′ are mixtures of [t] and [s]. Replace D by D_s for c sbar, and use the indicated antiparticle and charge for cbar q. No resonance masses are assigned.
The “charmed” sector includes D and D_s, particles and antiparticles. It includes the three 2S templates printed in Table IV. It does not include all possible D-wave parents: 1¹D₂ and 1³D₂ are absent from the printed template list. Thresholds are deliberately not imposed.

## The 28 flavor routes for any one template X -> Y + P

For positive charm there are 5 + 5 + 4 = 14; their distinct charge conjugates give another 14. X and Y stand for whichever spectroscopic parent and daughter are selected below.

| Route | Parent flavor slot | Daughter flavor slot | Emitted meson |
|---|---|---|---|
| R01 | D⁰ | D⁰ | π⁰ |
| R02 | D⁰ | D⁰ | η |
| R03 | D⁰ | D⁰ | η′ |
| R04 | D⁰ | D⁺ | π⁻ |
| R05 | D⁰ | D_s⁺ | K⁻ |
| R06 | D⁺ | D⁰ | π⁺ |
| R07 | D⁺ | D⁺ | π⁰ |
| R08 | D⁺ | D⁺ | η |
| R09 | D⁺ | D⁺ | η′ |
| R10 | D⁺ | D_s⁺ | K̄⁰ |
| R11 | D_s⁺ | D⁰ | K⁺ |
| R12 | D_s⁺ | D⁺ | K⁰ |
| R13 | D_s⁺ | D_s⁺ | η |
| R14 | D_s⁺ | D_s⁺ | η′ |
| R15 | D̄⁰ | D̄⁰ | π⁰ |
| R16 | D̄⁰ | D̄⁰ | η |
| R17 | D̄⁰ | D̄⁰ | η′ |
| R18 | D̄⁰ | D⁻ | π⁺ |
| R19 | D̄⁰ | D_s⁻ | K⁺ |
| R20 | D⁻ | D̄⁰ | π⁻ |
| R21 | D⁻ | D⁻ | π⁰ |
| R22 | D⁻ | D⁻ | η |
| R23 | D⁻ | D⁻ | η′ |
| R24 | D⁻ | D_s⁻ | K⁰ |
| R25 | D_s⁻ | D̄⁰ | K⁻ |
| R26 | D_s⁻ | D⁻ | K̄⁰ |
| R27 | D_s⁻ | D_s⁻ | η |
| R28 | D_s⁻ | D_s⁻ | η′ |

For example, the first five routes mean X⁰ -> Y⁰π⁰, Y⁰η, Y⁰η′, Y⁺π⁻, Y_s⁺K⁻. For X_s⁺ there are Y⁰K⁺, Y⁺K⁰, Y_s⁺η, Y_s⁺η′. The channel X_s⁺ -> Y_s⁺π⁰ vanishes in exact isospin and is not counted.
Only the light antiquark emits for c qbar; only the light quark emits for q cbar. There is no additional factor of two for emitter topology. Neutral π⁰, η and η′ are each counted once as a meson, not once per flavor component. A heavy daughter and a light emitted pseudoscalar are distinct, so their interchange contributes no extra charm channel.

## All 16 templates (18 partial waves)

| ID | Representative neutral channel | Relative waves | Table-IV entries | Final channels | Partial-wave channels |
|---|---|---|---|---:|---:|
| T01 | D*⁰ -> D⁰ + P | P | A1 | 28 | 28 |
| T02 | D₂*⁰ -> D⁰ + P | D | A2 | 28 | 28 |
| T03 | D₂*⁰ -> D*⁰ + P | D | A3 | 28 | 28 |
| T04 | D₁[t]⁰ -> D*⁰ + P | S, D | A4, S1 | 28 | 56 |
| T05 | D₁′[s]⁰ -> D*⁰ + P | S, D | A5, S2 | 28 | 56 |
| T06 | D(1³D₃)⁰ -> D⁰ + P | F | A6 | 28 | 28 |
| T07 | D(1³D₃)⁰ -> D*⁰ + P | F | A7 | 28 | 28 |
| T08 | D₁′[s]⁰ -> D₀*⁰ + P | P | Aprime1 | 28 | 28 |
| T09 | D(1³D₃)⁰ -> D₁′[s]⁰ + P | D | Adoubleprime1 | 28 | 28 |
| T10 | D₁[t]⁰ -> D₀*⁰ + P | P | A01 | 28 | 28 |
| T11 | D₀*⁰ -> D⁰ + P | S | S3 | 28 | 28 |
| T12 | D(1³D₁)⁰ -> D⁰ + P | P | D1 | 28 | 28 |
| T13 | D(1³D₁)⁰ -> D*⁰ + P | P | D2 | 28 | 28 |
| T14 | D(2¹S₀)⁰ -> D*⁰ + P | P | P1 | 28 | 28 |
| T15 | D(2³S₁)⁰ -> D⁰ + P | P | P2 | 28 | 28 |
| T16 | D(2³S₁)⁰ -> D*⁰ + P | P | P3 | 28 | 28 |
| Total | | | | **448** | **504** |

Rows T04 and T05 are the only double-wave rows. T08 and T10 have different axial basis parents; T06 and T12 have J=3 and J=1 parents, respectively. None are duplicates.
The original table also prints two 1³D₃ -> 1S entries with relative L=1; both fail the J triangle and neither is counted here. The retained 1³D₃ -> D and D* entries have L=3.

## Useful reductions of the large numbers

| Scope | Routes per template | Partial-wave entries | Distinct basis final channels |
|---|---:|---:|---:|
| Full nonet, charm and anticharm | 28 | 504 | 448 |
| Full nonet, positive charm only | 14 | 252 | 224 |
| π, K, η emission, both conjugates | 22 | 396 | 352 |
| π, η emission, both conjugates | 14 | 252 | 224 |
| Full nonet, both conjugates, omit the three 2S templates | 28 | 420 | 364 |

Charge-conjugate channels are separately named channels, not independent reduced amplitudes. Exact isospin further relates charge channels. Neither 504 nor 448 counts independent model parameters or independent experimental constraints.

## What changes when D₁ and D₁′ mean physical mesons?

Write D₁ = cos(θ) D₁[t] + sin(θ) D₁′[s], and D₁′ = -sin(θ) D₁[t] + cos(θ) D₁′[s], separately in each light-flavor sector. The phase convention is illustrative; only the two-dimensional span matters here.
T04/T05 already cover both axial basis parents -> D* with S and D waves; rotation leaves two physical parents with these wave possibilities. T08/T10 likewise already cover both axial parents -> D₀*.
**T09 is different:** A-double-prime includes only D(1³D₃) -> D₁′[s] + P with relative L=2. A generic nonzero mixing angle projects this daughter onto both physical D₁ and D₁′. Thus the represented subspace supplies two physical daughter channels per route, not one.
Consequently a generic mixing completion of this restricted template set gives **19 × 28 = 532 partial-wave channels** and **17 × 28 = 476 distinct final channels**. This is a projection/counting statement, not a numerical prediction: a complete amplitude can also require omitted basis matrix elements. Special mixing angles or dynamical cancellations can remove channels. It is still not an exhaustive through-D survey.
For a paper, label 504/448 explicitly as an **unmixed-basis census**. Do not attach those totals unqualified to named physical D₁/D₁′ resonances.

## Fully expanded ledger: all 448 basis final channels

Each line below is one distinct final channel. The bracket lists every retained relative partial wave; two-wave lines account for the extra 56 matrix elements.

### T01: D*⁰ -> D⁰ + P

- T01/R01: D*⁰ -> D⁰ + π⁰ [P]
- T01/R02: D*⁰ -> D⁰ + η [P]
- T01/R03: D*⁰ -> D⁰ + η′ [P]
- T01/R04: D*⁰ -> D⁺ + π⁻ [P]
- T01/R05: D*⁰ -> D_s⁺ + K⁻ [P]
- T01/R06: D*⁺ -> D⁰ + π⁺ [P]
- T01/R07: D*⁺ -> D⁺ + π⁰ [P]
- T01/R08: D*⁺ -> D⁺ + η [P]
- T01/R09: D*⁺ -> D⁺ + η′ [P]
- T01/R10: D*⁺ -> D_s⁺ + K̄⁰ [P]
- T01/R11: D_s*⁺ -> D⁰ + K⁺ [P]
- T01/R12: D_s*⁺ -> D⁺ + K⁰ [P]
- T01/R13: D_s*⁺ -> D_s⁺ + η [P]
- T01/R14: D_s*⁺ -> D_s⁺ + η′ [P]
- T01/R15: D̄*⁰ -> D̄⁰ + π⁰ [P]
- T01/R16: D̄*⁰ -> D̄⁰ + η [P]
- T01/R17: D̄*⁰ -> D̄⁰ + η′ [P]
- T01/R18: D̄*⁰ -> D⁻ + π⁺ [P]
- T01/R19: D̄*⁰ -> D_s⁻ + K⁺ [P]
- T01/R20: D*⁻ -> D̄⁰ + π⁻ [P]
- T01/R21: D*⁻ -> D⁻ + π⁰ [P]
- T01/R22: D*⁻ -> D⁻ + η [P]
- T01/R23: D*⁻ -> D⁻ + η′ [P]
- T01/R24: D*⁻ -> D_s⁻ + K⁰ [P]
- T01/R25: D_s*⁻ -> D̄⁰ + K⁻ [P]
- T01/R26: D_s*⁻ -> D⁻ + K̄⁰ [P]
- T01/R27: D_s*⁻ -> D_s⁻ + η [P]
- T01/R28: D_s*⁻ -> D_s⁻ + η′ [P]

### T02: D₂*⁰ -> D⁰ + P

- T02/R01: D₂*⁰ -> D⁰ + π⁰ [D]
- T02/R02: D₂*⁰ -> D⁰ + η [D]
- T02/R03: D₂*⁰ -> D⁰ + η′ [D]
- T02/R04: D₂*⁰ -> D⁺ + π⁻ [D]
- T02/R05: D₂*⁰ -> D_s⁺ + K⁻ [D]
- T02/R06: D₂*⁺ -> D⁰ + π⁺ [D]
- T02/R07: D₂*⁺ -> D⁺ + π⁰ [D]
- T02/R08: D₂*⁺ -> D⁺ + η [D]
- T02/R09: D₂*⁺ -> D⁺ + η′ [D]
- T02/R10: D₂*⁺ -> D_s⁺ + K̄⁰ [D]
- T02/R11: D_s₂*⁺ -> D⁰ + K⁺ [D]
- T02/R12: D_s₂*⁺ -> D⁺ + K⁰ [D]
- T02/R13: D_s₂*⁺ -> D_s⁺ + η [D]
- T02/R14: D_s₂*⁺ -> D_s⁺ + η′ [D]
- T02/R15: D̄₂*⁰ -> D̄⁰ + π⁰ [D]
- T02/R16: D̄₂*⁰ -> D̄⁰ + η [D]
- T02/R17: D̄₂*⁰ -> D̄⁰ + η′ [D]
- T02/R18: D̄₂*⁰ -> D⁻ + π⁺ [D]
- T02/R19: D̄₂*⁰ -> D_s⁻ + K⁺ [D]
- T02/R20: D₂*⁻ -> D̄⁰ + π⁻ [D]
- T02/R21: D₂*⁻ -> D⁻ + π⁰ [D]
- T02/R22: D₂*⁻ -> D⁻ + η [D]
- T02/R23: D₂*⁻ -> D⁻ + η′ [D]
- T02/R24: D₂*⁻ -> D_s⁻ + K⁰ [D]
- T02/R25: D_s₂*⁻ -> D̄⁰ + K⁻ [D]
- T02/R26: D_s₂*⁻ -> D⁻ + K̄⁰ [D]
- T02/R27: D_s₂*⁻ -> D_s⁻ + η [D]
- T02/R28: D_s₂*⁻ -> D_s⁻ + η′ [D]

### T03: D₂*⁰ -> D*⁰ + P

- T03/R01: D₂*⁰ -> D*⁰ + π⁰ [D]
- T03/R02: D₂*⁰ -> D*⁰ + η [D]
- T03/R03: D₂*⁰ -> D*⁰ + η′ [D]
- T03/R04: D₂*⁰ -> D*⁺ + π⁻ [D]
- T03/R05: D₂*⁰ -> D_s*⁺ + K⁻ [D]
- T03/R06: D₂*⁺ -> D*⁰ + π⁺ [D]
- T03/R07: D₂*⁺ -> D*⁺ + π⁰ [D]
- T03/R08: D₂*⁺ -> D*⁺ + η [D]
- T03/R09: D₂*⁺ -> D*⁺ + η′ [D]
- T03/R10: D₂*⁺ -> D_s*⁺ + K̄⁰ [D]
- T03/R11: D_s₂*⁺ -> D*⁰ + K⁺ [D]
- T03/R12: D_s₂*⁺ -> D*⁺ + K⁰ [D]
- T03/R13: D_s₂*⁺ -> D_s*⁺ + η [D]
- T03/R14: D_s₂*⁺ -> D_s*⁺ + η′ [D]
- T03/R15: D̄₂*⁰ -> D̄*⁰ + π⁰ [D]
- T03/R16: D̄₂*⁰ -> D̄*⁰ + η [D]
- T03/R17: D̄₂*⁰ -> D̄*⁰ + η′ [D]
- T03/R18: D̄₂*⁰ -> D*⁻ + π⁺ [D]
- T03/R19: D̄₂*⁰ -> D_s*⁻ + K⁺ [D]
- T03/R20: D₂*⁻ -> D̄*⁰ + π⁻ [D]
- T03/R21: D₂*⁻ -> D*⁻ + π⁰ [D]
- T03/R22: D₂*⁻ -> D*⁻ + η [D]
- T03/R23: D₂*⁻ -> D*⁻ + η′ [D]
- T03/R24: D₂*⁻ -> D_s*⁻ + K⁰ [D]
- T03/R25: D_s₂*⁻ -> D̄*⁰ + K⁻ [D]
- T03/R26: D_s₂*⁻ -> D*⁻ + K̄⁰ [D]
- T03/R27: D_s₂*⁻ -> D_s*⁻ + η [D]
- T03/R28: D_s₂*⁻ -> D_s*⁻ + η′ [D]

### T04: D₁[t]⁰ -> D*⁰ + P

- T04/R01: D₁[t]⁰ -> D*⁰ + π⁰ [S, D]
- T04/R02: D₁[t]⁰ -> D*⁰ + η [S, D]
- T04/R03: D₁[t]⁰ -> D*⁰ + η′ [S, D]
- T04/R04: D₁[t]⁰ -> D*⁺ + π⁻ [S, D]
- T04/R05: D₁[t]⁰ -> D_s*⁺ + K⁻ [S, D]
- T04/R06: D₁[t]⁺ -> D*⁰ + π⁺ [S, D]
- T04/R07: D₁[t]⁺ -> D*⁺ + π⁰ [S, D]
- T04/R08: D₁[t]⁺ -> D*⁺ + η [S, D]
- T04/R09: D₁[t]⁺ -> D*⁺ + η′ [S, D]
- T04/R10: D₁[t]⁺ -> D_s*⁺ + K̄⁰ [S, D]
- T04/R11: D_s₁[t]⁺ -> D*⁰ + K⁺ [S, D]
- T04/R12: D_s₁[t]⁺ -> D*⁺ + K⁰ [S, D]
- T04/R13: D_s₁[t]⁺ -> D_s*⁺ + η [S, D]
- T04/R14: D_s₁[t]⁺ -> D_s*⁺ + η′ [S, D]
- T04/R15: D̄₁[t]⁰ -> D̄*⁰ + π⁰ [S, D]
- T04/R16: D̄₁[t]⁰ -> D̄*⁰ + η [S, D]
- T04/R17: D̄₁[t]⁰ -> D̄*⁰ + η′ [S, D]
- T04/R18: D̄₁[t]⁰ -> D*⁻ + π⁺ [S, D]
- T04/R19: D̄₁[t]⁰ -> D_s*⁻ + K⁺ [S, D]
- T04/R20: D₁[t]⁻ -> D̄*⁰ + π⁻ [S, D]
- T04/R21: D₁[t]⁻ -> D*⁻ + π⁰ [S, D]
- T04/R22: D₁[t]⁻ -> D*⁻ + η [S, D]
- T04/R23: D₁[t]⁻ -> D*⁻ + η′ [S, D]
- T04/R24: D₁[t]⁻ -> D_s*⁻ + K⁰ [S, D]
- T04/R25: D_s₁[t]⁻ -> D̄*⁰ + K⁻ [S, D]
- T04/R26: D_s₁[t]⁻ -> D*⁻ + K̄⁰ [S, D]
- T04/R27: D_s₁[t]⁻ -> D_s*⁻ + η [S, D]
- T04/R28: D_s₁[t]⁻ -> D_s*⁻ + η′ [S, D]

### T05: D₁′[s]⁰ -> D*⁰ + P

- T05/R01: D₁′[s]⁰ -> D*⁰ + π⁰ [S, D]
- T05/R02: D₁′[s]⁰ -> D*⁰ + η [S, D]
- T05/R03: D₁′[s]⁰ -> D*⁰ + η′ [S, D]
- T05/R04: D₁′[s]⁰ -> D*⁺ + π⁻ [S, D]
- T05/R05: D₁′[s]⁰ -> D_s*⁺ + K⁻ [S, D]
- T05/R06: D₁′[s]⁺ -> D*⁰ + π⁺ [S, D]
- T05/R07: D₁′[s]⁺ -> D*⁺ + π⁰ [S, D]
- T05/R08: D₁′[s]⁺ -> D*⁺ + η [S, D]
- T05/R09: D₁′[s]⁺ -> D*⁺ + η′ [S, D]
- T05/R10: D₁′[s]⁺ -> D_s*⁺ + K̄⁰ [S, D]
- T05/R11: D_s₁′[s]⁺ -> D*⁰ + K⁺ [S, D]
- T05/R12: D_s₁′[s]⁺ -> D*⁺ + K⁰ [S, D]
- T05/R13: D_s₁′[s]⁺ -> D_s*⁺ + η [S, D]
- T05/R14: D_s₁′[s]⁺ -> D_s*⁺ + η′ [S, D]
- T05/R15: D̄₁′[s]⁰ -> D̄*⁰ + π⁰ [S, D]
- T05/R16: D̄₁′[s]⁰ -> D̄*⁰ + η [S, D]
- T05/R17: D̄₁′[s]⁰ -> D̄*⁰ + η′ [S, D]
- T05/R18: D̄₁′[s]⁰ -> D*⁻ + π⁺ [S, D]
- T05/R19: D̄₁′[s]⁰ -> D_s*⁻ + K⁺ [S, D]
- T05/R20: D₁′[s]⁻ -> D̄*⁰ + π⁻ [S, D]
- T05/R21: D₁′[s]⁻ -> D*⁻ + π⁰ [S, D]
- T05/R22: D₁′[s]⁻ -> D*⁻ + η [S, D]
- T05/R23: D₁′[s]⁻ -> D*⁻ + η′ [S, D]
- T05/R24: D₁′[s]⁻ -> D_s*⁻ + K⁰ [S, D]
- T05/R25: D_s₁′[s]⁻ -> D̄*⁰ + K⁻ [S, D]
- T05/R26: D_s₁′[s]⁻ -> D*⁻ + K̄⁰ [S, D]
- T05/R27: D_s₁′[s]⁻ -> D_s*⁻ + η [S, D]
- T05/R28: D_s₁′[s]⁻ -> D_s*⁻ + η′ [S, D]

### T06: D(1³D₃)⁰ -> D⁰ + P

- T06/R01: D(1³D₃)⁰ -> D⁰ + π⁰ [F]
- T06/R02: D(1³D₃)⁰ -> D⁰ + η [F]
- T06/R03: D(1³D₃)⁰ -> D⁰ + η′ [F]
- T06/R04: D(1³D₃)⁰ -> D⁺ + π⁻ [F]
- T06/R05: D(1³D₃)⁰ -> D_s⁺ + K⁻ [F]
- T06/R06: D(1³D₃)⁺ -> D⁰ + π⁺ [F]
- T06/R07: D(1³D₃)⁺ -> D⁺ + π⁰ [F]
- T06/R08: D(1³D₃)⁺ -> D⁺ + η [F]
- T06/R09: D(1³D₃)⁺ -> D⁺ + η′ [F]
- T06/R10: D(1³D₃)⁺ -> D_s⁺ + K̄⁰ [F]
- T06/R11: D_s(1³D₃)⁺ -> D⁰ + K⁺ [F]
- T06/R12: D_s(1³D₃)⁺ -> D⁺ + K⁰ [F]
- T06/R13: D_s(1³D₃)⁺ -> D_s⁺ + η [F]
- T06/R14: D_s(1³D₃)⁺ -> D_s⁺ + η′ [F]
- T06/R15: D̄(1³D₃)⁰ -> D̄⁰ + π⁰ [F]
- T06/R16: D̄(1³D₃)⁰ -> D̄⁰ + η [F]
- T06/R17: D̄(1³D₃)⁰ -> D̄⁰ + η′ [F]
- T06/R18: D̄(1³D₃)⁰ -> D⁻ + π⁺ [F]
- T06/R19: D̄(1³D₃)⁰ -> D_s⁻ + K⁺ [F]
- T06/R20: D(1³D₃)⁻ -> D̄⁰ + π⁻ [F]
- T06/R21: D(1³D₃)⁻ -> D⁻ + π⁰ [F]
- T06/R22: D(1³D₃)⁻ -> D⁻ + η [F]
- T06/R23: D(1³D₃)⁻ -> D⁻ + η′ [F]
- T06/R24: D(1³D₃)⁻ -> D_s⁻ + K⁰ [F]
- T06/R25: D_s(1³D₃)⁻ -> D̄⁰ + K⁻ [F]
- T06/R26: D_s(1³D₃)⁻ -> D⁻ + K̄⁰ [F]
- T06/R27: D_s(1³D₃)⁻ -> D_s⁻ + η [F]
- T06/R28: D_s(1³D₃)⁻ -> D_s⁻ + η′ [F]

### T07: D(1³D₃)⁰ -> D*⁰ + P

- T07/R01: D(1³D₃)⁰ -> D*⁰ + π⁰ [F]
- T07/R02: D(1³D₃)⁰ -> D*⁰ + η [F]
- T07/R03: D(1³D₃)⁰ -> D*⁰ + η′ [F]
- T07/R04: D(1³D₃)⁰ -> D*⁺ + π⁻ [F]
- T07/R05: D(1³D₃)⁰ -> D_s*⁺ + K⁻ [F]
- T07/R06: D(1³D₃)⁺ -> D*⁰ + π⁺ [F]
- T07/R07: D(1³D₃)⁺ -> D*⁺ + π⁰ [F]
- T07/R08: D(1³D₃)⁺ -> D*⁺ + η [F]
- T07/R09: D(1³D₃)⁺ -> D*⁺ + η′ [F]
- T07/R10: D(1³D₃)⁺ -> D_s*⁺ + K̄⁰ [F]
- T07/R11: D_s(1³D₃)⁺ -> D*⁰ + K⁺ [F]
- T07/R12: D_s(1³D₃)⁺ -> D*⁺ + K⁰ [F]
- T07/R13: D_s(1³D₃)⁺ -> D_s*⁺ + η [F]
- T07/R14: D_s(1³D₃)⁺ -> D_s*⁺ + η′ [F]
- T07/R15: D̄(1³D₃)⁰ -> D̄*⁰ + π⁰ [F]
- T07/R16: D̄(1³D₃)⁰ -> D̄*⁰ + η [F]
- T07/R17: D̄(1³D₃)⁰ -> D̄*⁰ + η′ [F]
- T07/R18: D̄(1³D₃)⁰ -> D*⁻ + π⁺ [F]
- T07/R19: D̄(1³D₃)⁰ -> D_s*⁻ + K⁺ [F]
- T07/R20: D(1³D₃)⁻ -> D̄*⁰ + π⁻ [F]
- T07/R21: D(1³D₃)⁻ -> D*⁻ + π⁰ [F]
- T07/R22: D(1³D₃)⁻ -> D*⁻ + η [F]
- T07/R23: D(1³D₃)⁻ -> D*⁻ + η′ [F]
- T07/R24: D(1³D₃)⁻ -> D_s*⁻ + K⁰ [F]
- T07/R25: D_s(1³D₃)⁻ -> D̄*⁰ + K⁻ [F]
- T07/R26: D_s(1³D₃)⁻ -> D*⁻ + K̄⁰ [F]
- T07/R27: D_s(1³D₃)⁻ -> D_s*⁻ + η [F]
- T07/R28: D_s(1³D₃)⁻ -> D_s*⁻ + η′ [F]

### T08: D₁′[s]⁰ -> D₀*⁰ + P

- T08/R01: D₁′[s]⁰ -> D₀*⁰ + π⁰ [P]
- T08/R02: D₁′[s]⁰ -> D₀*⁰ + η [P]
- T08/R03: D₁′[s]⁰ -> D₀*⁰ + η′ [P]
- T08/R04: D₁′[s]⁰ -> D₀*⁺ + π⁻ [P]
- T08/R05: D₁′[s]⁰ -> D_s₀*⁺ + K⁻ [P]
- T08/R06: D₁′[s]⁺ -> D₀*⁰ + π⁺ [P]
- T08/R07: D₁′[s]⁺ -> D₀*⁺ + π⁰ [P]
- T08/R08: D₁′[s]⁺ -> D₀*⁺ + η [P]
- T08/R09: D₁′[s]⁺ -> D₀*⁺ + η′ [P]
- T08/R10: D₁′[s]⁺ -> D_s₀*⁺ + K̄⁰ [P]
- T08/R11: D_s₁′[s]⁺ -> D₀*⁰ + K⁺ [P]
- T08/R12: D_s₁′[s]⁺ -> D₀*⁺ + K⁰ [P]
- T08/R13: D_s₁′[s]⁺ -> D_s₀*⁺ + η [P]
- T08/R14: D_s₁′[s]⁺ -> D_s₀*⁺ + η′ [P]
- T08/R15: D̄₁′[s]⁰ -> D̄₀*⁰ + π⁰ [P]
- T08/R16: D̄₁′[s]⁰ -> D̄₀*⁰ + η [P]
- T08/R17: D̄₁′[s]⁰ -> D̄₀*⁰ + η′ [P]
- T08/R18: D̄₁′[s]⁰ -> D₀*⁻ + π⁺ [P]
- T08/R19: D̄₁′[s]⁰ -> D_s₀*⁻ + K⁺ [P]
- T08/R20: D₁′[s]⁻ -> D̄₀*⁰ + π⁻ [P]
- T08/R21: D₁′[s]⁻ -> D₀*⁻ + π⁰ [P]
- T08/R22: D₁′[s]⁻ -> D₀*⁻ + η [P]
- T08/R23: D₁′[s]⁻ -> D₀*⁻ + η′ [P]
- T08/R24: D₁′[s]⁻ -> D_s₀*⁻ + K⁰ [P]
- T08/R25: D_s₁′[s]⁻ -> D̄₀*⁰ + K⁻ [P]
- T08/R26: D_s₁′[s]⁻ -> D₀*⁻ + K̄⁰ [P]
- T08/R27: D_s₁′[s]⁻ -> D_s₀*⁻ + η [P]
- T08/R28: D_s₁′[s]⁻ -> D_s₀*⁻ + η′ [P]

### T09: D(1³D₃)⁰ -> D₁′[s]⁰ + P

- T09/R01: D(1³D₃)⁰ -> D₁′[s]⁰ + π⁰ [D]
- T09/R02: D(1³D₃)⁰ -> D₁′[s]⁰ + η [D]
- T09/R03: D(1³D₃)⁰ -> D₁′[s]⁰ + η′ [D]
- T09/R04: D(1³D₃)⁰ -> D₁′[s]⁺ + π⁻ [D]
- T09/R05: D(1³D₃)⁰ -> D_s₁′[s]⁺ + K⁻ [D]
- T09/R06: D(1³D₃)⁺ -> D₁′[s]⁰ + π⁺ [D]
- T09/R07: D(1³D₃)⁺ -> D₁′[s]⁺ + π⁰ [D]
- T09/R08: D(1³D₃)⁺ -> D₁′[s]⁺ + η [D]
- T09/R09: D(1³D₃)⁺ -> D₁′[s]⁺ + η′ [D]
- T09/R10: D(1³D₃)⁺ -> D_s₁′[s]⁺ + K̄⁰ [D]
- T09/R11: D_s(1³D₃)⁺ -> D₁′[s]⁰ + K⁺ [D]
- T09/R12: D_s(1³D₃)⁺ -> D₁′[s]⁺ + K⁰ [D]
- T09/R13: D_s(1³D₃)⁺ -> D_s₁′[s]⁺ + η [D]
- T09/R14: D_s(1³D₃)⁺ -> D_s₁′[s]⁺ + η′ [D]
- T09/R15: D̄(1³D₃)⁰ -> D̄₁′[s]⁰ + π⁰ [D]
- T09/R16: D̄(1³D₃)⁰ -> D̄₁′[s]⁰ + η [D]
- T09/R17: D̄(1³D₃)⁰ -> D̄₁′[s]⁰ + η′ [D]
- T09/R18: D̄(1³D₃)⁰ -> D₁′[s]⁻ + π⁺ [D]
- T09/R19: D̄(1³D₃)⁰ -> D_s₁′[s]⁻ + K⁺ [D]
- T09/R20: D(1³D₃)⁻ -> D̄₁′[s]⁰ + π⁻ [D]
- T09/R21: D(1³D₃)⁻ -> D₁′[s]⁻ + π⁰ [D]
- T09/R22: D(1³D₃)⁻ -> D₁′[s]⁻ + η [D]
- T09/R23: D(1³D₃)⁻ -> D₁′[s]⁻ + η′ [D]
- T09/R24: D(1³D₃)⁻ -> D_s₁′[s]⁻ + K⁰ [D]
- T09/R25: D_s(1³D₃)⁻ -> D̄₁′[s]⁰ + K⁻ [D]
- T09/R26: D_s(1³D₃)⁻ -> D₁′[s]⁻ + K̄⁰ [D]
- T09/R27: D_s(1³D₃)⁻ -> D_s₁′[s]⁻ + η [D]
- T09/R28: D_s(1³D₃)⁻ -> D_s₁′[s]⁻ + η′ [D]

### T10: D₁[t]⁰ -> D₀*⁰ + P

- T10/R01: D₁[t]⁰ -> D₀*⁰ + π⁰ [P]
- T10/R02: D₁[t]⁰ -> D₀*⁰ + η [P]
- T10/R03: D₁[t]⁰ -> D₀*⁰ + η′ [P]
- T10/R04: D₁[t]⁰ -> D₀*⁺ + π⁻ [P]
- T10/R05: D₁[t]⁰ -> D_s₀*⁺ + K⁻ [P]
- T10/R06: D₁[t]⁺ -> D₀*⁰ + π⁺ [P]
- T10/R07: D₁[t]⁺ -> D₀*⁺ + π⁰ [P]
- T10/R08: D₁[t]⁺ -> D₀*⁺ + η [P]
- T10/R09: D₁[t]⁺ -> D₀*⁺ + η′ [P]
- T10/R10: D₁[t]⁺ -> D_s₀*⁺ + K̄⁰ [P]
- T10/R11: D_s₁[t]⁺ -> D₀*⁰ + K⁺ [P]
- T10/R12: D_s₁[t]⁺ -> D₀*⁺ + K⁰ [P]
- T10/R13: D_s₁[t]⁺ -> D_s₀*⁺ + η [P]
- T10/R14: D_s₁[t]⁺ -> D_s₀*⁺ + η′ [P]
- T10/R15: D̄₁[t]⁰ -> D̄₀*⁰ + π⁰ [P]
- T10/R16: D̄₁[t]⁰ -> D̄₀*⁰ + η [P]
- T10/R17: D̄₁[t]⁰ -> D̄₀*⁰ + η′ [P]
- T10/R18: D̄₁[t]⁰ -> D₀*⁻ + π⁺ [P]
- T10/R19: D̄₁[t]⁰ -> D_s₀*⁻ + K⁺ [P]
- T10/R20: D₁[t]⁻ -> D̄₀*⁰ + π⁻ [P]
- T10/R21: D₁[t]⁻ -> D₀*⁻ + π⁰ [P]
- T10/R22: D₁[t]⁻ -> D₀*⁻ + η [P]
- T10/R23: D₁[t]⁻ -> D₀*⁻ + η′ [P]
- T10/R24: D₁[t]⁻ -> D_s₀*⁻ + K⁰ [P]
- T10/R25: D_s₁[t]⁻ -> D̄₀*⁰ + K⁻ [P]
- T10/R26: D_s₁[t]⁻ -> D₀*⁻ + K̄⁰ [P]
- T10/R27: D_s₁[t]⁻ -> D_s₀*⁻ + η [P]
- T10/R28: D_s₁[t]⁻ -> D_s₀*⁻ + η′ [P]

### T11: D₀*⁰ -> D⁰ + P

- T11/R01: D₀*⁰ -> D⁰ + π⁰ [S]
- T11/R02: D₀*⁰ -> D⁰ + η [S]
- T11/R03: D₀*⁰ -> D⁰ + η′ [S]
- T11/R04: D₀*⁰ -> D⁺ + π⁻ [S]
- T11/R05: D₀*⁰ -> D_s⁺ + K⁻ [S]
- T11/R06: D₀*⁺ -> D⁰ + π⁺ [S]
- T11/R07: D₀*⁺ -> D⁺ + π⁰ [S]
- T11/R08: D₀*⁺ -> D⁺ + η [S]
- T11/R09: D₀*⁺ -> D⁺ + η′ [S]
- T11/R10: D₀*⁺ -> D_s⁺ + K̄⁰ [S]
- T11/R11: D_s₀*⁺ -> D⁰ + K⁺ [S]
- T11/R12: D_s₀*⁺ -> D⁺ + K⁰ [S]
- T11/R13: D_s₀*⁺ -> D_s⁺ + η [S]
- T11/R14: D_s₀*⁺ -> D_s⁺ + η′ [S]
- T11/R15: D̄₀*⁰ -> D̄⁰ + π⁰ [S]
- T11/R16: D̄₀*⁰ -> D̄⁰ + η [S]
- T11/R17: D̄₀*⁰ -> D̄⁰ + η′ [S]
- T11/R18: D̄₀*⁰ -> D⁻ + π⁺ [S]
- T11/R19: D̄₀*⁰ -> D_s⁻ + K⁺ [S]
- T11/R20: D₀*⁻ -> D̄⁰ + π⁻ [S]
- T11/R21: D₀*⁻ -> D⁻ + π⁰ [S]
- T11/R22: D₀*⁻ -> D⁻ + η [S]
- T11/R23: D₀*⁻ -> D⁻ + η′ [S]
- T11/R24: D₀*⁻ -> D_s⁻ + K⁰ [S]
- T11/R25: D_s₀*⁻ -> D̄⁰ + K⁻ [S]
- T11/R26: D_s₀*⁻ -> D⁻ + K̄⁰ [S]
- T11/R27: D_s₀*⁻ -> D_s⁻ + η [S]
- T11/R28: D_s₀*⁻ -> D_s⁻ + η′ [S]

### T12: D(1³D₁)⁰ -> D⁰ + P

- T12/R01: D(1³D₁)⁰ -> D⁰ + π⁰ [P]
- T12/R02: D(1³D₁)⁰ -> D⁰ + η [P]
- T12/R03: D(1³D₁)⁰ -> D⁰ + η′ [P]
- T12/R04: D(1³D₁)⁰ -> D⁺ + π⁻ [P]
- T12/R05: D(1³D₁)⁰ -> D_s⁺ + K⁻ [P]
- T12/R06: D(1³D₁)⁺ -> D⁰ + π⁺ [P]
- T12/R07: D(1³D₁)⁺ -> D⁺ + π⁰ [P]
- T12/R08: D(1³D₁)⁺ -> D⁺ + η [P]
- T12/R09: D(1³D₁)⁺ -> D⁺ + η′ [P]
- T12/R10: D(1³D₁)⁺ -> D_s⁺ + K̄⁰ [P]
- T12/R11: D_s(1³D₁)⁺ -> D⁰ + K⁺ [P]
- T12/R12: D_s(1³D₁)⁺ -> D⁺ + K⁰ [P]
- T12/R13: D_s(1³D₁)⁺ -> D_s⁺ + η [P]
- T12/R14: D_s(1³D₁)⁺ -> D_s⁺ + η′ [P]
- T12/R15: D̄(1³D₁)⁰ -> D̄⁰ + π⁰ [P]
- T12/R16: D̄(1³D₁)⁰ -> D̄⁰ + η [P]
- T12/R17: D̄(1³D₁)⁰ -> D̄⁰ + η′ [P]
- T12/R18: D̄(1³D₁)⁰ -> D⁻ + π⁺ [P]
- T12/R19: D̄(1³D₁)⁰ -> D_s⁻ + K⁺ [P]
- T12/R20: D(1³D₁)⁻ -> D̄⁰ + π⁻ [P]
- T12/R21: D(1³D₁)⁻ -> D⁻ + π⁰ [P]
- T12/R22: D(1³D₁)⁻ -> D⁻ + η [P]
- T12/R23: D(1³D₁)⁻ -> D⁻ + η′ [P]
- T12/R24: D(1³D₁)⁻ -> D_s⁻ + K⁰ [P]
- T12/R25: D_s(1³D₁)⁻ -> D̄⁰ + K⁻ [P]
- T12/R26: D_s(1³D₁)⁻ -> D⁻ + K̄⁰ [P]
- T12/R27: D_s(1³D₁)⁻ -> D_s⁻ + η [P]
- T12/R28: D_s(1³D₁)⁻ -> D_s⁻ + η′ [P]

### T13: D(1³D₁)⁰ -> D*⁰ + P

- T13/R01: D(1³D₁)⁰ -> D*⁰ + π⁰ [P]
- T13/R02: D(1³D₁)⁰ -> D*⁰ + η [P]
- T13/R03: D(1³D₁)⁰ -> D*⁰ + η′ [P]
- T13/R04: D(1³D₁)⁰ -> D*⁺ + π⁻ [P]
- T13/R05: D(1³D₁)⁰ -> D_s*⁺ + K⁻ [P]
- T13/R06: D(1³D₁)⁺ -> D*⁰ + π⁺ [P]
- T13/R07: D(1³D₁)⁺ -> D*⁺ + π⁰ [P]
- T13/R08: D(1³D₁)⁺ -> D*⁺ + η [P]
- T13/R09: D(1³D₁)⁺ -> D*⁺ + η′ [P]
- T13/R10: D(1³D₁)⁺ -> D_s*⁺ + K̄⁰ [P]
- T13/R11: D_s(1³D₁)⁺ -> D*⁰ + K⁺ [P]
- T13/R12: D_s(1³D₁)⁺ -> D*⁺ + K⁰ [P]
- T13/R13: D_s(1³D₁)⁺ -> D_s*⁺ + η [P]
- T13/R14: D_s(1³D₁)⁺ -> D_s*⁺ + η′ [P]
- T13/R15: D̄(1³D₁)⁰ -> D̄*⁰ + π⁰ [P]
- T13/R16: D̄(1³D₁)⁰ -> D̄*⁰ + η [P]
- T13/R17: D̄(1³D₁)⁰ -> D̄*⁰ + η′ [P]
- T13/R18: D̄(1³D₁)⁰ -> D*⁻ + π⁺ [P]
- T13/R19: D̄(1³D₁)⁰ -> D_s*⁻ + K⁺ [P]
- T13/R20: D(1³D₁)⁻ -> D̄*⁰ + π⁻ [P]
- T13/R21: D(1³D₁)⁻ -> D*⁻ + π⁰ [P]
- T13/R22: D(1³D₁)⁻ -> D*⁻ + η [P]
- T13/R23: D(1³D₁)⁻ -> D*⁻ + η′ [P]
- T13/R24: D(1³D₁)⁻ -> D_s*⁻ + K⁰ [P]
- T13/R25: D_s(1³D₁)⁻ -> D̄*⁰ + K⁻ [P]
- T13/R26: D_s(1³D₁)⁻ -> D*⁻ + K̄⁰ [P]
- T13/R27: D_s(1³D₁)⁻ -> D_s*⁻ + η [P]
- T13/R28: D_s(1³D₁)⁻ -> D_s*⁻ + η′ [P]

### T14: D(2¹S₀)⁰ -> D*⁰ + P

- T14/R01: D(2¹S₀)⁰ -> D*⁰ + π⁰ [P]
- T14/R02: D(2¹S₀)⁰ -> D*⁰ + η [P]
- T14/R03: D(2¹S₀)⁰ -> D*⁰ + η′ [P]
- T14/R04: D(2¹S₀)⁰ -> D*⁺ + π⁻ [P]
- T14/R05: D(2¹S₀)⁰ -> D_s*⁺ + K⁻ [P]
- T14/R06: D(2¹S₀)⁺ -> D*⁰ + π⁺ [P]
- T14/R07: D(2¹S₀)⁺ -> D*⁺ + π⁰ [P]
- T14/R08: D(2¹S₀)⁺ -> D*⁺ + η [P]
- T14/R09: D(2¹S₀)⁺ -> D*⁺ + η′ [P]
- T14/R10: D(2¹S₀)⁺ -> D_s*⁺ + K̄⁰ [P]
- T14/R11: D_s(2¹S₀)⁺ -> D*⁰ + K⁺ [P]
- T14/R12: D_s(2¹S₀)⁺ -> D*⁺ + K⁰ [P]
- T14/R13: D_s(2¹S₀)⁺ -> D_s*⁺ + η [P]
- T14/R14: D_s(2¹S₀)⁺ -> D_s*⁺ + η′ [P]
- T14/R15: D̄(2¹S₀)⁰ -> D̄*⁰ + π⁰ [P]
- T14/R16: D̄(2¹S₀)⁰ -> D̄*⁰ + η [P]
- T14/R17: D̄(2¹S₀)⁰ -> D̄*⁰ + η′ [P]
- T14/R18: D̄(2¹S₀)⁰ -> D*⁻ + π⁺ [P]
- T14/R19: D̄(2¹S₀)⁰ -> D_s*⁻ + K⁺ [P]
- T14/R20: D(2¹S₀)⁻ -> D̄*⁰ + π⁻ [P]
- T14/R21: D(2¹S₀)⁻ -> D*⁻ + π⁰ [P]
- T14/R22: D(2¹S₀)⁻ -> D*⁻ + η [P]
- T14/R23: D(2¹S₀)⁻ -> D*⁻ + η′ [P]
- T14/R24: D(2¹S₀)⁻ -> D_s*⁻ + K⁰ [P]
- T14/R25: D_s(2¹S₀)⁻ -> D̄*⁰ + K⁻ [P]
- T14/R26: D_s(2¹S₀)⁻ -> D*⁻ + K̄⁰ [P]
- T14/R27: D_s(2¹S₀)⁻ -> D_s*⁻ + η [P]
- T14/R28: D_s(2¹S₀)⁻ -> D_s*⁻ + η′ [P]

### T15: D(2³S₁)⁰ -> D⁰ + P

- T15/R01: D(2³S₁)⁰ -> D⁰ + π⁰ [P]
- T15/R02: D(2³S₁)⁰ -> D⁰ + η [P]
- T15/R03: D(2³S₁)⁰ -> D⁰ + η′ [P]
- T15/R04: D(2³S₁)⁰ -> D⁺ + π⁻ [P]
- T15/R05: D(2³S₁)⁰ -> D_s⁺ + K⁻ [P]
- T15/R06: D(2³S₁)⁺ -> D⁰ + π⁺ [P]
- T15/R07: D(2³S₁)⁺ -> D⁺ + π⁰ [P]
- T15/R08: D(2³S₁)⁺ -> D⁺ + η [P]
- T15/R09: D(2³S₁)⁺ -> D⁺ + η′ [P]
- T15/R10: D(2³S₁)⁺ -> D_s⁺ + K̄⁰ [P]
- T15/R11: D_s(2³S₁)⁺ -> D⁰ + K⁺ [P]
- T15/R12: D_s(2³S₁)⁺ -> D⁺ + K⁰ [P]
- T15/R13: D_s(2³S₁)⁺ -> D_s⁺ + η [P]
- T15/R14: D_s(2³S₁)⁺ -> D_s⁺ + η′ [P]
- T15/R15: D̄(2³S₁)⁰ -> D̄⁰ + π⁰ [P]
- T15/R16: D̄(2³S₁)⁰ -> D̄⁰ + η [P]
- T15/R17: D̄(2³S₁)⁰ -> D̄⁰ + η′ [P]
- T15/R18: D̄(2³S₁)⁰ -> D⁻ + π⁺ [P]
- T15/R19: D̄(2³S₁)⁰ -> D_s⁻ + K⁺ [P]
- T15/R20: D(2³S₁)⁻ -> D̄⁰ + π⁻ [P]
- T15/R21: D(2³S₁)⁻ -> D⁻ + π⁰ [P]
- T15/R22: D(2³S₁)⁻ -> D⁻ + η [P]
- T15/R23: D(2³S₁)⁻ -> D⁻ + η′ [P]
- T15/R24: D(2³S₁)⁻ -> D_s⁻ + K⁰ [P]
- T15/R25: D_s(2³S₁)⁻ -> D̄⁰ + K⁻ [P]
- T15/R26: D_s(2³S₁)⁻ -> D⁻ + K̄⁰ [P]
- T15/R27: D_s(2³S₁)⁻ -> D_s⁻ + η [P]
- T15/R28: D_s(2³S₁)⁻ -> D_s⁻ + η′ [P]

### T16: D(2³S₁)⁰ -> D*⁰ + P

- T16/R01: D(2³S₁)⁰ -> D*⁰ + π⁰ [P]
- T16/R02: D(2³S₁)⁰ -> D*⁰ + η [P]
- T16/R03: D(2³S₁)⁰ -> D*⁰ + η′ [P]
- T16/R04: D(2³S₁)⁰ -> D*⁺ + π⁻ [P]
- T16/R05: D(2³S₁)⁰ -> D_s*⁺ + K⁻ [P]
- T16/R06: D(2³S₁)⁺ -> D*⁰ + π⁺ [P]
- T16/R07: D(2³S₁)⁺ -> D*⁺ + π⁰ [P]
- T16/R08: D(2³S₁)⁺ -> D*⁺ + η [P]
- T16/R09: D(2³S₁)⁺ -> D*⁺ + η′ [P]
- T16/R10: D(2³S₁)⁺ -> D_s*⁺ + K̄⁰ [P]
- T16/R11: D_s(2³S₁)⁺ -> D*⁰ + K⁺ [P]
- T16/R12: D_s(2³S₁)⁺ -> D*⁺ + K⁰ [P]
- T16/R13: D_s(2³S₁)⁺ -> D_s*⁺ + η [P]
- T16/R14: D_s(2³S₁)⁺ -> D_s*⁺ + η′ [P]
- T16/R15: D̄(2³S₁)⁰ -> D̄*⁰ + π⁰ [P]
- T16/R16: D̄(2³S₁)⁰ -> D̄*⁰ + η [P]
- T16/R17: D̄(2³S₁)⁰ -> D̄*⁰ + η′ [P]
- T16/R18: D̄(2³S₁)⁰ -> D*⁻ + π⁺ [P]
- T16/R19: D̄(2³S₁)⁰ -> D_s*⁻ + K⁺ [P]
- T16/R20: D(2³S₁)⁻ -> D̄*⁰ + π⁻ [P]
- T16/R21: D(2³S₁)⁻ -> D*⁻ + π⁰ [P]
- T16/R22: D(2³S₁)⁻ -> D*⁻ + η [P]
- T16/R23: D(2³S₁)⁻ -> D*⁻ + η′ [P]
- T16/R24: D(2³S₁)⁻ -> D_s*⁻ + K⁰ [P]
- T16/R25: D_s(2³S₁)⁻ -> D̄*⁰ + K⁻ [P]
- T16/R26: D_s(2³S₁)⁻ -> D*⁻ + K̄⁰ [P]
- T16/R27: D_s(2³S₁)⁻ -> D_s*⁻ + η [P]
- T16/R28: D_s(2³S₁)⁻ -> D_s*⁻ + η′ [P]
