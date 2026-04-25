# Bottomonium Baseline Residuals

Baseline model: radial finite-difference solver with `relativistic` kinetic energy, with smeared S-wave contact hyperfine, GI Table II masses, `b`, `c`, and Fig. 2 running Coulomb ansatz.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 9.400 | 9.471 |   +70.7 | high |
| `2^1S_0` | 9.980 | 10.028 |   +48.1 | high |
| `3^1S_0` | 10.340 | 10.380 |   +40.0 | high |
| `1^3S_1` | 9.460 | 9.522 |   +62.4 | high |
| `2^3S_1` | 10.000 | 10.058 |   +58.0 | high |
| `1^3D_1` | 10.140 | 10.176 |   +36.4 | high |
| `3^3S_1` | 10.350 | 10.404 |   +54.3 | high |
| `2^3D_1` | 10.440 | 10.481 |   +41.5 | high |
| `4^3S_1` | 10.630 | 10.682 |   +52.0 | medium |
| `3^3D_1` | 10.700 | 10.740 |   +39.6 | medium |
| `5^3S_1` | 10.880 | 10.923 |   +42.9 | medium |
| `6^3S_1` | 11.100 | 11.140 |   +40.0 | medium |
| `1^1P_1` | 9.880 | 9.920 |   +40.1 | high |
| `2^1P_1` | 10.250 | 10.291 |   +40.6 | high |
| `1^3P_0` | 9.850 | 9.920 |   +70.1 | high |
| `2^3P_0` | 10.230 | 10.291 |   +60.6 | high |
| `1^3P_1` | 9.880 | 9.920 |   +40.1 | high |
| `2^3P_1` | 10.250 | 10.291 |   +40.6 | high |
| `1^3P_2` | 9.900 | 9.920 |   +20.1 | high |
| `2^3P_2` | 10.260 | 10.291 |   +30.6 | high |
| `1^3F_2` | 10.350 | 10.379 |   +28.6 | medium |
| `1^1D_2` | 10.150 | 10.176 |   +26.4 | high |
| `2^1D_2` | 10.450 | 10.481 |   +31.5 | high |
| `1^3D_2` | 10.150 | 10.176 |   +26.4 | high |
| `2^3D_2` | 10.450 | 10.481 |   +31.5 | high |
| `1^3D_3` | 10.160 | 10.176 |   +16.4 | high |
| `2^3D_3` | 10.450 | 10.481 |   +31.5 | high |
| `1^1F_3` | 10.350 | 10.379 |   +28.6 | medium |
| `1^3F_3` | 10.350 | 10.379 |   +28.6 | medium |
| `1^3F_4` | 10.360 | 10.379 |   +18.6 | medium |

Mean absolute residual: 39.9 MeV.
Max absolute residual: 70.7 MeV.

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
