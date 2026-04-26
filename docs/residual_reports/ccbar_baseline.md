# Charmonium Baseline Residuals

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 2.970 | 3.099 |  +129.0 | high |
| `2^1S_0` | 3.620 | 3.712 |   +92.4 | high |
| `3^1S_0` | 4.060 | 4.137 |   +77.3 | high |
| `1^3S_1` | 3.100 | 3.197 |   +96.6 | high |
| `2^3S_1` | 3.680 | 3.786 |  +106.1 | high |
| `1^3D_1` | 3.820 | 3.892 |   +72.2 | high |
| `3^3S_1` | 4.100 | 4.203 |  +103.1 | high |
| `2^3D_1` | 4.190 | 4.268 |   +78.1 | high |
| `4^3S_1` | 4.450 | 4.549 |   +98.8 | medium |
| `3^3D_1` | 4.520 | 4.593 |   +72.6 | medium |
| `1^1P_1` | 3.520 | 3.601 |   +81.4 | high |
| `2^1P_1` | 3.960 | 4.047 |   +86.6 | high |
| `1^3P_0` | 3.440 | 3.557 |  +116.8 | high |
| `2^3P_0` | 3.920 | 4.003 |   +82.8 | high |
| `1^3P_1` | 3.510 | 3.596 |   +85.9 | high |
| `2^3P_1` | 3.950 | 4.041 |   +91.0 | high |
| `1^3P_2` | 3.550 | 3.614 |   +63.7 | high |
| `2^3P_2` | 3.980 | 4.059 |   +78.8 | high |
| `1^3F_2` | 4.090 | 4.150 |   +59.7 | medium |
| `1^1D_2` | 3.840 | 3.901 |   +61.5 | high |
| `2^1D_2` | 4.210 | 4.279 |   +69.0 | high |
| `1^3D_2` | 3.840 | 3.902 |   +61.8 | high |
| `2^3D_2` | 4.210 | 4.279 |   +68.9 | high |
| `1^3D_3` | 3.850 | 3.905 |   +55.2 | high |
| `2^3D_3` | 4.220 | 4.284 |   +63.7 | high |
| `1^1F_3` | 4.090 | 4.149 |   +59.3 | medium |
| `1^3F_3` | 4.100 | 4.151 |   +50.9 | medium |
| `1^3F_4` | 4.090 | 4.148 |   +57.9 | medium |

Mean absolute residual: 79.3 MeV.
Max absolute residual: 129.0 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 3.841 | 3.901 |   +61.0 |
| `2D` | 4 | 4.211 | 4.279 |   +68.5 |
| `3D` | 1 | 4.520 | 4.593 |   +72.6 |
| `1F` | 4 | 4.093 | 4.149 |   +56.8 |
| `1P` | 4 | 3.523 | 3.601 |   +78.1 |
| `2P` | 4 | 3.962 | 4.047 |   +84.1 |
| `1S` | 2 | 3.068 | 3.172 |  +104.7 |
| `2S` | 2 | 3.665 | 3.768 |  +102.7 |
| `3S` | 2 | 4.090 | 4.187 |   +96.7 |
| `4S` | 1 | 4.450 | 4.549 |   +98.8 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
