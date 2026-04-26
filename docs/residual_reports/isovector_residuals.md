# Residuals: isovector (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.150 | -0.219 |  -368.5 | high |
| `2^1S_0` | 1.300 | 0.221 | -1078.8 | high |
| `1^3S_1` | 0.770 | 1.297 |  +527.3 | high |
| `2^3S_1` | 1.450 | 2.083 |  +633.2 | high |
| `1^3D_1` | 1.660 | 2.216 |  +556.5 | high |
| `3^1S_0` | 1.880 | 0.723 | -1157.4 | high |
| `3^3S_1` | 2.000 | 2.609 |  +609.1 | high |
| `2^3D_1` | 2.150 | 2.594 |  +443.5 | high |
| `1^1P_1` | 1.220 | 1.445 |  +225.1 | high |
| `2^1P_1` | 1.780 | 1.979 |  +198.6 | high |
| `1^3P_0` | 1.090 | 1.559 |  +469.3 | high |
| `2^3P_0` | 1.780 | 1.974 |  +194.5 | high |
| `1^3P_1` | 1.240 | 1.622 |  +382.4 | high |
| `2^3P_1` | 1.820 | 2.092 |  +271.5 | high |
| `1^3P_2` | 1.310 | 1.316 |    +5.9 | high |
| `2^3P_2` | 1.820 | 1.912 |   +91.6 | high |
| `1^3F_2` | 2.050 | 2.657 |  +607.5 | high |
| `1^1D_2` | 1.680 | 1.834 |  +154.1 | high |
| `2^1D_2` | 2.130 | 2.284 |  +153.9 | high |
| `1^3D_2` | 1.700 | 2.004 |  +303.8 | high |
| `2^3D_2` | 2.150 | 2.434 |  +283.9 | high |
| `1^3D_3` | 1.680 | 1.549 |  -130.9 | high |
| `2^3D_3` | 2.130 | 2.044 |   -86.0 | high |
| `1^3G_3` | 2.370 | 3.024 |  +654.0 | high |
| `1^1F_3` | 2.030 | 2.158 |  +128.4 | high |
| `1^3F_3` | 2.050 | 2.306 |  +256.1 | high |
| `1^3F_4` | 2.010 | 1.766 |  -243.7 | high |
| `1^1G_4` | 2.330 | 2.443 |  +112.7 | high |
| `1^3G_4` | 2.340 | 2.574 |  +233.6 | high |
| `1^3G_5` | 2.300 | 1.966 |  -334.4 | high |

Mean absolute residual: 363.2 MeV.
Max absolute residual: 1157.4 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.682 | 1.834 |  +152.1 |
| `2D` | 4 | 2.138 | 2.284 |  +145.9 |
| `1F` | 4 | 2.032 | 2.158 |  +126.3 |
| `1G` | 4 | 2.331 | 2.443 |  +111.6 |
| `1P` | 4 | 1.252 | 1.445 |  +193.4 |
| `2P` | 4 | 1.807 | 1.979 |  +171.9 |
| `1S` | 2 | 0.615 | 0.918 |  +303.4 |
| `2S` | 2 | 1.412 | 1.618 |  +205.2 |
| `3S` | 2 | 1.970 | 2.137 |  +167.4 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
