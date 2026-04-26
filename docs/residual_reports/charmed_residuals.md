# Residuals: charmed (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.843 |   -37.2 | high |
| `2^1S_0` | 2.580 | 2.501 |   -79.2 | high |
| `1^3S_1` | 2.040 | 2.165 |  +124.9 | high |
| `2^3S_1` | 2.640 | 2.789 |  +149.1 | high |
| `1^3D_1` | 2.820 | 2.873 |   +52.6 | high |
| `1^3P_0` | 2.400 | 2.130 |  -269.6 | high |
| `1^1P_1` | 2.440 | 2.561 |  +121.0 | high |
| `1^3P_1` | 2.490 | 2.415 |   -75.2 | high |
| `1^3P_2` | 2.500 | 2.735 |  +234.8 | high |
| `1^3D_3` | 2.830 | 2.915 |   +85.0 | high |
| `1^3F_4` | 3.110 | 3.098 |   -12.0 | medium |
| `1^1S_0` | 1.980 | 2.039 |   +58.9 | high |
| `2^1S_0` | 2.670 | 2.688 |   +17.5 | high |
| `1^3S_1` | 2.130 | 2.241 |  +111.1 | high |
| `2^3S_1` | 2.730 | 2.864 |  +134.4 | high |
| `1^3D_1` | 2.900 | 2.971 |   +70.7 | high |
| `1^3P_0` | 2.480 | 2.495 |   +14.9 | high |
| `1^1P_1` | 2.530 | 2.654 |  +124.3 | high |
| `1^3P_1` | 2.570 | 2.615 |   +44.6 | high |
| `1^3P_2` | 2.590 | 2.710 |  +120.0 | high |
| `1^3D_3` | 2.920 | 2.994 |   +74.2 | high |
| `1^3F_4` | 3.190 | 3.241 |   +51.2 | medium |

Mean absolute residual: 93.8 MeV.
Max absolute residual: 269.6 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.945 |   +74.2 |
| `1F` | 2 | 3.150 | 3.170 |   +19.6 |
| `1P` | 8 | 2.518 | 2.608 |   +90.2 |
| `1S` | 4 | 2.046 | 2.137 |   +91.2 |
| `2S` | 4 | 2.670 | 2.769 |   +98.6 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
