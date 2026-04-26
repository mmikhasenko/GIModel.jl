# Residuals: isoscalar (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.520 | -0.219 |  -738.5 | medium |
| `1^1S_0` | 0.960 | -0.219 | -1178.5 | medium |
| `2^1S_0` | 1.440 | 0.221 | -1218.8 | medium |
| `2^1S_0` | 1.630 | 0.221 | -1408.8 | medium |
| `1^3S_1` | 0.780 | 1.297 |  +517.3 | medium |
| `1^3S_1` | 1.020 | 1.297 |  +277.3 | medium |
| `2^3S_1` | 1.460 | 2.083 |  +623.2 | medium |
| `1^3D_1` | 1.660 | 1.404 |  -255.6 | medium |
| `2^3S_1` | 1.690 | 2.083 |  +393.2 | medium |
| `1^3D_1` | 1.880 | 1.404 |  -475.6 | medium |
| `1^1P_1` | 1.220 | 1.445 |  +225.1 | medium |
| `1^1P_1` | 1.470 | 1.445 |   -24.9 | medium |
| `2^1P_1` | 1.780 | 1.979 |  +198.6 | medium |
| `2^1P_1` | 2.010 | 1.979 |   -31.4 | medium |
| `1^3P_0` | 1.090 | 1.091 |    +0.6 | medium |
| `1^3P_0` | 1.360 | 1.091 |  -269.4 | medium |
| `2^3P_0` | 1.780 | 1.753 |   -27.3 | medium |
| `2^3P_0` | 1.990 | 1.753 |  -237.3 | medium |
| `1^3P_1` | 1.240 | 1.358 |  +118.0 | medium |
| `1^3P_1` | 1.480 | 1.358 |  -122.0 | medium |
| `2^3P_1` | 1.820 | 1.952 |  +131.9 | medium |
| `2^3P_1` | 2.030 | 1.952 |   -78.1 | medium |
| `1^3P_2` | 1.280 | 1.556 |  +276.3 | medium |
| `1^3P_2` | 1.530 | 1.556 |   +26.3 | medium |
| `2^3P_2` | 1.820 | 2.028 |  +208.3 | medium |
| `2^3P_2` | 2.040 | 2.028 |   -11.7 | medium |
| `1^3F_2` | 2.050 | 1.648 |  -402.1 | medium |
| `1^3F_2` | 2.240 | 1.648 |  -592.1 | medium |
| `1^1D_2` | 1.680 | 1.834 |  +154.1 | medium |
| `1^1D_2` | 1.890 | 1.834 |   -55.9 | medium |
| `1^3D_2` | 1.700 | 1.704 |    +4.1 | medium |
| `1^3D_2` | 1.910 | 1.704 |  -205.9 | medium |
| `2^3D_2` | 2.260 | 2.178 |   -82.3 | low |
| `1^3D_3` | 1.680 | 2.109 |  +428.5 | medium |
| `1^3D_3` | 1.900 | 2.109 |  +208.5 | medium |
| `1^3G_3` | 2.370 | 1.853 |  -516.5 | medium |
| `1^3G_3` | 2.540 | 1.853 |  -686.5 | medium |
| `1^1F_3` | 2.030 | 2.158 |  +128.4 | medium |
| `1^1F_3` | 2.220 | 2.158 |   -61.6 | medium |
| `1^3F_3` | 2.050 | 2.028 |   -22.1 | medium |
| `1^3F_3` | 2.230 | 2.028 |  -202.1 | medium |
| `1^3F_4` | 2.010 | 2.545 |  +535.1 | medium |
| `1^3F_4` | 2.200 | 2.545 |  +345.1 | medium |
| `1^1G_4` | 2.330 | 2.443 |  +112.7 | medium |
| `1^1G_4` | 2.510 | 2.443 |   -67.3 | medium |
| `1^3G_4` | 2.340 | 2.324 |   -16.1 | medium |
| `1^3G_4` | 2.520 | 2.324 |  -196.1 | medium |
| `1^3G_5` | 2.300 | 2.915 |  +615.3 | medium |
| `1^3G_5` | 2.470 | 2.915 |  +445.3 | medium |

Mean absolute residual: 309.3 MeV.
Max absolute residual: 1408.8 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 8 | 1.789 | 1.833 |   +43.7 |
| `2D` | 1 | 2.260 | 2.178 |   -82.3 |
| `1F` | 8 | 2.126 | 2.159 |   +33.0 |
| `1G` | 8 | 2.419 | 2.443 |   +24.2 |
| `1P` | 8 | 1.364 | 1.440 |   +76.3 |
| `2P` | 8 | 1.916 | 1.974 |   +57.5 |
| `1S` | 4 | 0.860 | 0.918 |   +58.4 |
| `2S` | 4 | 1.565 | 1.618 |   +52.7 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
