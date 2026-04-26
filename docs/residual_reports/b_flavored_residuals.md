# Residuals: b_flavored (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.236 |   -73.8 | high |
| `2^1S_0` | 5.900 | 5.838 |   -61.5 | high |
| `1^3S_1` | 5.370 | 5.386 |   +15.7 | high |
| `2^3S_1` | 5.930 | 5.967 |   +36.6 | high |
| `1^3P_2` | 5.800 | 6.104 |  +303.5 | high |
| `1^3D_3` | 6.110 | 6.213 |  +102.8 | medium |
| `1^3F_4` | 6.360 | 6.351 |    -8.6 | medium |
| `1^1S_0` | 5.390 | 5.373 |   -17.3 | high |
| `2^1S_0` | 5.980 | 5.970 |   -10.2 | high |
| `1^3S_1` | 5.450 | 5.470 |   +20.3 | high |
| `2^3S_1` | 6.010 | 6.052 |   +41.9 | high |
| `1^3P_2` | 5.880 | 5.988 |  +108.0 | high |
| `1^3D_3` | 6.180 | 6.238 |   +58.0 | medium |
| `1^3F_4` | 6.430 | 6.459 |   +28.7 | medium |
| `1^1S_0` | 6.270 | 6.335 |   +64.8 | high |
| `2^1S_0` | 6.850 | 6.910 |   +59.6 | high |
| `1^3S_1` | 6.340 | 6.392 |   +52.2 | high |
| `2^3S_1` | 6.890 | 6.950 |   +59.9 | high |
| `1^3P_2` | 6.770 | 6.808 |   +37.8 | high |
| `1^3D_3` | 7.040 | 7.080 |   +39.7 | medium |
| `1^3F_4` | 7.270 | 7.302 |   +32.4 | medium |

Mean absolute residual: 58.7 MeV.
Max absolute residual: 303.5 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.510 |   +66.8 |
| `1F` | 3 | 6.687 | 6.704 |   +17.5 |
| `1P` | 3 | 6.150 | 6.300 |  +149.8 |
| `1S` | 6 | 5.704 | 5.724 |   +19.9 |
| `2S` | 6 | 6.268 | 6.302 |   +33.6 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
