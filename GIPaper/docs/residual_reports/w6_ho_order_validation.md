# W6 — paper-order (finite HO-basis) validation of the spin-distorted waves

Scores the Table VII gluonic subtable under four treatments of the
spin-dependent operators (smeared contact for S-waves, calibrated
spin-orbit + tensor for ³P_J), against the spin-independent central-wave
baseline. Zero new parameters: the operators are the spectrum-calibrated
blocks (`contact_hyperfine_operator`, `fine_structure_grid_operator` with
k_spin_orbit/k_tensor). The **paper-order** treatment — full diagonalization
of `H_central + V_spin` in the finite paper-β HO basis
(`ho_full_distorted_states`) — is what the harmonized Table VII audit uses
for every subtable.

## 1. Basis fidelity control: central-wave S_L, HO vs FD

The paper-style HO diagonalization (24 states, paper β convention)
reproduces the FD smeared wavefunction-at-origin to ≤3%. The 15-20% row
residuals of the central-wave audit are NOT a basis-fidelity artifact —
and truncation lowers S_L, the wrong direction to explain them.

| flavor | L | n | S_FD | S_HO | S_HO/S_FD |
|:-:|:-:|:-:|---:|---:|:-:|
| c | 0 | 1 | +0.2820 | +0.2807 | 0.9955 |
| c | 0 | 2 | -0.2020 | -0.2007 | 0.9936 |
| c | 1 | 1 | +0.1438 | +0.1432 | 0.9956 |
| b | 0 | 1 | +0.7267 | +0.7097 | 0.9765 |
| b | 0 | 2 | -0.5171 | -0.5032 | 0.9732 |
| b | 0 | 3 | +0.4486 | +0.4356 | 0.9708 |
| b | 0 | 4 | -0.4117 | -0.3996 | 0.9706 |
| b | 1 | 1 | +0.1997 | +0.1958 | 0.9806 |
| b | 1 | 2 | -0.2114 | -0.2067 | 0.9777 |

## 1b. What order is "paper order"? The light ¹S₀ mass discriminator

The residuals are spin-dependent wavefunction distortion (§2). But *how*
the paper carries the distortion is settled by the light `¹S₀` nonstrange
mass, where the contact term is enormous. First-order PT and full
diagonalization agree for heavy spin splittings but diverge sharply here:

| treatment | pion ¹S₀ mass (GeV) |
|---|---:|
| first-order PT (`ho_first_order_distorted_states`) | 0.2844 |
| **finite-HO full diag (`ho_full_distorted_states`)** | **0.0968** |
| nonperturbative FD (fine grid) | 0.0957 |

First-order PT over-raises the pion to ≈0.28 GeV; the paper keeps it light
(≈0.10 GeV), which BOTH the finite-HO full diagonalization and the fine-grid
FD resummation reproduce. So the paper's mechanism is **full diagonalization**
(the contact is resummed, not truncated at first order); first-order PT is
only a heavy-quark proxy. The finite HO basis — not a perturbation order —
is why the heavy gluonic ratios sit *below* the fine-grid FD overshoot (§2):
a finite basis resums the contact less aggressively than the fine FD grid.

## 2. Gluonic subtable under the four treatments

`ratio` = |model|/|paper| per treatment. `central` = spin-independent wave;
`1st-PT` = first-order PT in the HO central eigenbasis; **`paper` = full
diagonalization in the finite HO basis** (the harmonized-audit treatment);
`nonpert` = the operator resummed on the fine FD grid.

| decay | paper amp | central | 1st-PT | **paper** | nonpert FD | M_paper |
|---|---:|:-:|:-:|:-:|:-:|---:|
| `eta_c -> 2g` | +4.700 | 0.87 | 1.06 | **1.09** | 1.11 | 2.960 |
| `psi -> 3g` | +0.420 | 1.11 | 1.02 | **1.02** | 1.02 | 3.091 |
| `eta'_c -> 2g` | -2.700 | 0.99 | 1.08 | **1.06** | 1.08 | 3.619 |
| `psi' -> 3g` | -0.280 | 1.06 | 1.01 | **1.01** | 1.01 | 3.680 |
| `eta_b -> 2g` | +2.500 | 0.98 | 1.11 | **1.14** | 1.21 | 9.394 |
| `Upsilon -> 3g` | +0.210 | 1.13 | 1.04 | **1.04** | 1.05 | 9.458 |
| `eta'_b -> 2g` | -1.700 | 1.01 | 1.10 | **1.10** | 1.18 | 9.974 |
| `Upsilon' -> 3g` | -0.150 | 1.11 | 1.03 | **1.03** | 1.05 | 10.005 |
| `Upsilon'' -> 3g` | +0.130 | 1.09 | 1.03 | **1.02** | 1.05 | 10.354 |
| `Upsilon''' -> 3g` | -0.110 | 1.18 | 1.11 | **1.10** | 1.13 | 10.633 |
| `chi_2c -> 2g` | +0.880 | 1.14 | 1.08 | **1.08** | 1.09 | 3.533 |
| `chi_0c -> 2g` | +2.500 | 0.78 | 0.93 | **0.95** | 0.95 | 3.458 |
| `chi_2b -> 2g` | +0.350 | 0.98 | 0.92 | **0.93** | 0.94 | 9.888 |
| `chi_0b -> 2g` | +0.820 | 0.81 | 0.95 | **0.98** | 1.01 | 9.854 |
| `chi'_2b -> 2g` | -0.370 | 0.98 | 0.92 | **0.92** | 0.94 | 10.255 |
| `chi'_0b -> 2g` | -0.820 | 0.85 | 0.96 | **0.98** | 1.01 | 10.231 |

- central: median 0.994, spread [0.78, 1.18]
- 1st-order PT: median 1.028, spread [0.92, 1.11]
- paper (full diag): median 1.021, spread [0.92, 1.14]
- nonpert FD: median 1.046, spread [0.94, 1.21]

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
audit. First-order PT reproduces the heavy gluonic ratios (its perturbative
limit) but is refuted as the mechanism by the light `¹S₀` mass (§1b): it
over-raises the pion, whereas full diagonalization keeps it light like the
FD resummation. The fine-grid FD *over*-resums (η_b 1.21 vs the finite-basis
1.14) — a grid-resolution effect, not perturbation order. The former
attribution of the charm rows to FD-vs-HO wavefunction-at-origin infidelity
is refuted by §1; the former attribution to a first-order treatment is
refuted by §1b.
