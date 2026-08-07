# W6 — paper-order (finite HO-basis) validation of the spin-distorted waves

Scores the Table VII gluonic subtable under three treatments of the
spin-dependent operators (smeared contact for S-waves, calibrated
spin-orbit + tensor for ³P_J), against the spin-independent central-wave
baseline. Zero new parameters: the operators are the spectrum-calibrated
blocks assembled natively by `fixed_channel_solution`. The **paper-order**
treatment is full diagonalization of the fixed-(L,S,J) Hamiltonian in
the finite paper-beta HO basis; no mesh operator is projected into HO
for every subtable.

## 1. Basis fidelity control: central-wave S_L, HO vs FD

The paper-style HO diagonalization (24 states, paper β convention)
reproduces the FD smeared wavefunction-at-origin to ≤3%. The 15-20% row
residuals of the central-wave audit are NOT a basis-fidelity artifact —
and truncation lowers S_L, the wrong direction to explain them.

| flavor | L | n | S_FD | S_HO | S_HO/S_FD |
|:-:|:-:|:-:|---:|---:|:-:|
| c | 0 | 1 | +0.2820 | +0.2807 | 0.9954 |
| c | 0 | 2 | -0.2020 | -0.2006 | 0.9935 |
| c | 1 | 1 | +0.1438 | +0.1432 | 0.9956 |
| b | 0 | 1 | +0.7267 | +0.7093 | 0.9761 |
| b | 0 | 2 | -0.5171 | -0.5029 | 0.9726 |
| b | 0 | 3 | +0.4486 | +0.4353 | 0.9702 |
| b | 0 | 4 | -0.4117 | -0.3993 | 0.9699 |
| b | 1 | 1 | +0.1997 | +0.1958 | 0.9804 |
| b | 1 | 2 | -0.2114 | -0.2066 | 0.9775 |

## 1b. Light ¹S₀ full-solve cross-check

The residuals are spin-dependent wavefunction distortion (§2). But *how*
the paper carries the distortion is settled by the light `¹S₀` nonstrange
mass, where the contact term is enormous. Both entries below are complete
fixed-channel diagonalizations with independent representations:

| treatment | pion ¹S₀ mass (GeV) |
|---|---:|
| **native finite-HO full diagonalization** | **0.0969** |
| independent fine-grid FD full diagonalization | 0.0957 |

The two methods reproduce the light pion with no HO-to-mesh fallback.

## 2. Gluonic subtable under the native full treatments

`ratio` = |model|/|paper| per treatment. `central` = spin-independent wave;
**`paper` = native full diagonalization in the finite HO basis**;
`FD` = the same fixed-channel Hamiltonian on the fine grid.

| decay | paper amp | central | **paper HO** | full FD | M_paper |
|---|---:|:-:|:-:|:-:|---:|
| `eta_c -> 2g` | +4.700 | 0.87 | **1.09** | 1.11 | 2.960 |
| `psi -> 3g` | +0.420 | 1.11 | **1.02** | 1.02 | 3.091 |
| `eta'_c -> 2g` | -2.700 | 0.99 | **1.06** | 1.08 | 3.619 |
| `psi' -> 3g` | -0.280 | 1.06 | **1.01** | 1.01 | 3.680 |
| `eta_b -> 2g` | +2.500 | 0.98 | **1.14** | 1.21 | 9.394 |
| `Upsilon -> 3g` | +0.210 | 1.13 | **1.03** | 1.05 | 9.459 |
| `eta'_b -> 2g` | -1.700 | 1.01 | **1.10** | 1.18 | 9.974 |
| `Upsilon' -> 3g` | -0.150 | 1.11 | **1.03** | 1.05 | 10.005 |
| `Upsilon'' -> 3g` | +0.130 | 1.09 | **1.02** | 1.05 | 10.354 |
| `Upsilon''' -> 3g` | -0.110 | 1.18 | **1.10** | 1.13 | 10.633 |
| `chi_2c -> 2g` | +0.880 | 1.14 | **1.08** | 1.09 | 3.533 |
| `chi_0c -> 2g` | +2.500 | 0.78 | **0.95** | 0.95 | 3.458 |
| `chi_2b -> 2g` | +0.350 | 0.98 | **0.93** | 0.94 | 9.888 |
| `chi_0b -> 2g` | +0.820 | 0.81 | **0.98** | 1.01 | 9.854 |
| `chi'_2b -> 2g` | -0.370 | 0.98 | **0.92** | 0.94 | 10.255 |
| `chi'_0b -> 2g` | -0.820 | 0.85 | **0.98** | 1.01 | 10.231 |

- central: median 0.994, spread [0.78, 1.18]
- paper HO (full diag): median 1.021, spread [0.92, 1.14]
- full FD: median 1.046, spread [0.94, 1.21]

## 3. Splitting-pattern collapse

Ratio-of-ratios inside a multiplet sharing one central wave; 1.00 means
the treatment carries the paper's full spin-splitting of the origin.

| pattern | central | paper (full diag) |
|---|:-:|:-:|
| `eta_c/psi` | 0.78 | 1.07 |
| `eta'_c/psi'` | 0.94 | 1.05 |
| `eta_b/Upsilon` | 0.87 | 1.10 |
| `eta'_b/Upsilon'` | 0.92 | 1.07 |
| `chi_0c/chi_2c` | 0.68 | 0.87 |
| `chi_0b/chi_2b` | 0.83 | 1.06 |
| `chi'_0b/chi'_2b` | 0.87 | 1.06 |

## Conclusion

The W4 residual structure (singlet low / triplet high, ³P₀ low / ³P₂ high)
is the paper's spin-dependent wavefunction distortion, carried by **full
diagonalization of `H_central + V_spin` in the finite paper-β HO basis** —
the paper's literal method, now used across the whole harmonized Table VII
audit. The light `¹S₀` mass (§1b) confirms that native HO and FD solve the
same resummed contact problem without an intermediate grid projection. The former
attribution of the charm rows to FD-vs-HO wavefunction-at-origin infidelity
is refuted by §1.
