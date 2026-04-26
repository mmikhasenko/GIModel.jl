# Residuals: charmed (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.843 |   -37.2 | high |
| `2^1S_0` | 2.580 | 2.501 |   -79.2 | high |
| `1^3S_1` | 2.040 | 2.165 |  +124.9 | high |
| `2^3S_1` | 2.640 | 2.789 |  +149.1 | high |
| `1^3D_1` | 2.820 | 2.680 |  -139.6 | high |
| `1^3P_0` | 2.400 | 2.445 |   +45.1 | high |
| `1^1P_1` | 2.440 | 2.561 |  +121.0 | high |
| `1^3P_1` | 2.490 | 2.525 |   +35.5 | high |
| `1^3P_2` | 2.500 | 2.602 |  +102.4 | high |
| `1^3D_3` | 2.830 | 3.051 |  +220.8 | high |
| `1^3F_4` | 3.110 | 3.404 |  +294.4 | medium |
| `1^1S_0` | 1.980 | 2.039 |   +58.9 | high |
| `2^1S_0` | 2.670 | 2.688 |   +17.5 | high |
| `1^3S_1` | 2.130 | 2.241 |  +111.1 | high |
| `2^3S_1` | 2.730 | 2.864 |  +134.4 | high |
| `1^3D_1` | 2.900 | 2.923 |   +22.7 | high |
| `1^3P_0` | 2.480 | 2.619 |  +138.5 | high |
| `1^1P_1` | 2.530 | 2.654 |  +124.3 | high |
| `1^3P_1` | 2.570 | 2.647 |   +76.7 | high |
| `1^3P_2` | 2.590 | 2.665 |   +74.7 | high |
| `1^3D_3` | 2.920 | 3.031 |  +110.7 | high |
| `1^3F_4` | 3.190 | 3.329 |  +139.2 | medium |

Mean absolute residual: 107.2 MeV.
Max absolute residual: 294.4 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.969 |   +98.5 |
| `1F` | 2 | 3.150 | 3.367 |  +216.8 |
| `1P` | 8 | 2.518 | 2.607 |   +89.2 |
| `1S` | 4 | 2.046 | 2.137 |   +91.2 |
| `2S` | 4 | 2.670 | 2.769 |   +98.6 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
