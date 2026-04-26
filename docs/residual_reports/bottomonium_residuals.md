# Residuals: bottomonium (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 9.400 | 9.477 |   +77.2 | high |
| `2^1S_0` | 9.980 | 10.032 |   +51.9 | high |
| `3^1S_0` | 10.340 | 10.383 |   +43.1 | high |
| `1^3S_1` | 9.460 | 9.520 |   +60.2 | high |
| `2^3S_1` | 10.000 | 10.057 |   +56.8 | high |
| `1^3D_1` | 10.140 | 10.171 |   +31.4 | high |
| `3^3S_1` | 10.350 | 10.403 |   +53.3 | high |
| `2^3D_1` | 10.440 | 10.477 |   +37.0 | high |
| `4^3S_1` | 10.630 | 10.681 |   +51.1 | medium |
| `3^3D_1` | 10.700 | 10.735 |   +35.4 | medium |
| `5^3S_1` | 10.880 | 10.922 |   +42.0 | medium |
| `6^3S_1` | 11.100 | 11.139 |   +39.2 | medium |
| `1^1P_1` | 9.880 | 9.920 |   +40.1 | high |
| `2^1P_1` | 10.250 | 10.291 |   +40.6 | high |
| `1^3P_0` | 9.850 | 9.903 |   +52.7 | high |
| `2^3P_0` | 10.230 | 10.277 |   +46.7 | high |
| `1^3P_1` | 9.880 | 9.917 |   +37.3 | high |
| `2^3P_1` | 10.250 | 10.288 |   +38.5 | high |
| `1^3P_2` | 9.900 | 9.925 |   +25.2 | high |
| `2^3P_2` | 10.260 | 10.295 |   +34.7 | high |
| `1^3F_2` | 10.350 | 10.376 |   +26.5 | medium |
| `1^1D_2` | 10.150 | 10.176 |   +26.4 | high |
| `2^1D_2` | 10.450 | 10.481 |   +31.5 | high |
| `1^3D_2` | 10.150 | 10.176 |   +25.9 | high |
| `2^3D_2` | 10.450 | 10.481 |   +31.0 | high |
| `1^3D_3` | 10.160 | 10.179 |   +19.0 | high |
| `2^3D_3` | 10.450 | 10.484 |   +33.7 | high |
| `1^1F_3` | 10.350 | 10.379 |   +28.6 | medium |
| `1^3F_3` | 10.350 | 10.379 |   +28.5 | medium |
| `1^3F_4` | 10.360 | 10.380 |   +19.8 | medium |

Mean absolute residual: 38.8 MeV.
Max absolute residual: 77.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 10.152 | 10.176 |   +24.4 |
| `2D` | 4 | 10.448 | 10.481 |   +33.0 |
| `3D` | 1 | 10.700 | 10.735 |   +35.4 |
| `1F` | 4 | 10.353 | 10.379 |   +25.4 |
| `1P` | 4 | 9.886 | 9.920 |   +34.2 |
| `2P` | 4 | 10.252 | 10.291 |   +38.1 |
| `1S` | 2 | 9.445 | 9.509 |   +64.4 |
| `2S` | 2 | 9.995 | 10.051 |   +55.5 |
| `3S` | 2 | 10.348 | 10.398 |   +50.7 |
| `4S` | 1 | 10.630 | 10.681 |   +51.1 |
| `5S` | 1 | 10.880 | 10.922 |   +42.0 |
| `6S` | 1 | 11.100 | 11.139 |   +39.2 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
