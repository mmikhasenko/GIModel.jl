# Residuals: isovector (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.150 | 0.095 |   -55.0 | high |
| `2^1S_0` | 1.300 | 1.279 |   -21.3 | high |
| `1^3S_1` | 0.770 | 0.760 |    -9.9 | high |
| `2^3S_1` | 1.450 | 1.458 |    +8.2 | high |
| `1^3D_1` | 1.660 | 1.664 |    +4.4 | high |
| `3^1S_0` | 1.880 | 1.859 |   -21.0 | high |
| `3^3S_1` | 2.000 | 2.003 |    +3.2 | high |
| `2^3D_1` | 2.150 | 2.137 |   -13.3 | high |
| `1^1P_1` | 1.220 | 1.258 |   +37.6 | high |
| `2^1P_1` | 1.780 | 1.807 |   +27.1 | high |
| `1^3P_0` | 1.090 | 1.098 |    +8.1 | high |
| `2^3P_0` | 1.780 | 1.756 |   -24.5 | high |
| `1^3P_1` | 1.240 | 1.266 |   +26.2 | high |
| `2^3P_1` | 1.820 | 1.816 |    -3.7 | high |
| `1^3P_2` | 1.310 | 1.279 |   -31.3 | high |
| `2^3P_2` | 1.820 | 1.812 |    -8.2 | high |
| `1^3F_2` | 2.050 | 2.038 |   -12.3 | high |
| `1^1D_2` | 1.680 | 1.687 |    +7.0 | high |
| `2^1D_2` | 2.130 | 2.142 |   +12.1 | high |
| `1^3D_2` | 1.700 | 1.699 |    -1.4 | high |
| `2^3D_2` | 2.150 | 2.151 |    +1.4 | high |
| `1^3D_3` | 1.680 | 1.687 |    +6.8 | high |
| `2^3D_3` | 2.130 | 2.137 |    +7.5 | high |
| `1^3G_3` | 2.370 | 2.347 |   -22.6 | high |
| `1^1F_3` | 2.030 | 2.034 |    +4.4 | high |
| `1^3F_3` | 2.050 | 2.044 |    -6.4 | high |
| `1^3F_4` | 2.010 | 2.024 |   +13.5 | high |
| `1^1G_4` | 2.330 | 2.334 |    +3.6 | high |
| `1^3G_4` | 2.340 | 2.341 |    +0.8 | high |
| `1^3G_5` | 2.300 | 2.318 |   +17.6 | high |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.838 |  -742.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -742.5 | 0.095 |
| `2^1S_0` | 1.383 |  -104.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -104.4 | 1.279 |
| `1^3S_1` | 0.671 |   +88.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +88.9 | 0.760 |
| `2^3S_1` | 1.413 |   +46.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +46.4 | 1.458 |
| `1^3D_1` | 1.691 |    +0.0 |   -77.7 |   +61.7 |   -16.0 |   -15.3 |    +0.0 |   -31.3 | 1.664 |
| `3^1S_0` | 1.939 |   -80.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -80.1 | 1.859 |
| `3^3S_1` | 1.963 |   +40.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +40.3 | 2.003 |
| `2^3D_1` | 2.143 |    +0.0 |   -47.5 |   +48.6 |    +1.1 |    -8.5 |    +0.0 |    -7.4 | 2.137 |
| `1^1P_1` | 1.258 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.258 |
| `2^1P_1` | 1.807 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.807 |
| `1^3P_0` | 1.272 |    +0.0 |  -165.0 |   +79.9 |   -85.1 |   -88.7 |    +0.0 |  -173.8 | 1.098 |
| `2^3P_0` | 1.807 |    +0.0 |   -71.1 |   +54.5 |   -16.6 |   -35.3 |    +0.0 |   -51.9 | 1.756 |
| `1^3P_1` | 1.259 |    +0.0 |   -75.2 |   +40.7 |   -34.4 |   +41.9 |    +0.0 |    +7.5 | 1.266 |
| `2^3P_1` | 1.807 |    +0.0 |   -35.0 |   +26.8 |    -8.2 |   +17.4 |    +0.0 |    +9.2 | 1.816 |
| `1^3P_2` | 1.260 |    +0.0 |   +68.7 |   -41.8 |   +26.9 |    -7.9 |    +0.0 |   +19.0 | 1.279 |
| `2^3P_2` | 1.807 |    +0.0 |   +34.9 |   -26.6 |    +8.3 |    -3.5 |    +0.0 |    +4.8 | 1.812 |
| `1^3F_2` | 2.036 |    +0.0 |   -46.1 |   +52.3 |    +6.3 |    -5.6 |    +0.0 |    +0.6 | 2.038 |
| `1^1D_2` | 1.687 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.687 |
| `2^1D_2` | 2.142 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.142 |
| `1^3D_2` | 1.687 |    +0.0 |   -24.3 |   +21.1 |    -3.2 |   +14.6 |    +0.0 |   +11.4 | 1.699 |
| `2^3D_2` | 2.142 |    +0.0 |   -15.4 |   +16.2 |    +0.8 |    +8.4 |    +0.0 |    +9.2 | 2.151 |
| `1^3D_3` | 1.688 |    +0.0 |   +45.7 |   -43.3 |    +2.3 |    -4.0 |    +0.0 |    -1.7 | 1.687 |
| `2^3D_3` | 2.142 |    +0.0 |   +30.1 |   -32.7 |    -2.6 |    -2.4 |    +0.0 |    -5.0 | 2.137 |
| `1^3G_3` | 2.335 |    +0.0 |   -31.1 |   +46.3 |   +15.2 |    -2.7 |    +0.0 |   +12.4 | 2.347 |
| `1^1F_3` | 2.034 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.034 |
| `1^3F_3` | 2.034 |    +0.0 |   -11.0 |   +13.4 |    +2.4 |    +6.8 |    +0.0 |    +9.2 | 2.044 |
| `1^3F_4` | 2.035 |    +0.0 |   +31.6 |   -41.2 |    -9.6 |    -2.2 |    +0.0 |   -11.8 | 2.024 |
| `1^1G_4` | 2.334 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.334 |
| `1^3G_4` | 2.334 |    +0.0 |    -6.0 |    +9.4 |    +3.4 |    +3.8 |    +0.0 |    +7.2 | 2.341 |
| `1^3G_5` | 2.334 |    +0.0 |   +23.3 |   -38.7 |   -15.4 |    -1.3 |    +0.0 |   -16.7 | 2.318 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 1.460 | 1.458 |   +16.5 |  -0.997 |  +0.080 |
| `1^3D_1` | 1.660 | 1.664 |   +16.5 |  +0.080 |  +0.995 |
| `3^3S_1` | 2.003 | 2.003 |    +3.6 |  +1.000 |  -0.027 |
| `2^3D_1` | 2.136 | 2.137 |    +3.6 |  -0.027 |  -0.999 |
| `2^3P_2` | 1.812 | 1.812 |    +4.2 |  -1.000 |  +0.018 |
| `1^3F_2` | 2.037 | 2.038 |    +4.2 |  +0.018 |  +0.999 |
| `2^3D_3` | 2.137 | 2.137 |    +1.2 |  -1.000 |  +0.006 |
| `1^3G_3` | 2.347 | 2.347 |    +1.2 |  +0.006 |  +1.000 |

Mean absolute residual: 14.0 MeV.
Max absolute residual: 55.0 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.682 | 1.686 |    +4.5 |
| `2D` | 4 | 2.138 | 2.142 |    +4.0 |
| `1F` | 4 | 2.032 | 2.034 |    +1.7 |
| `1G` | 4 | 2.331 | 2.333 |    +2.1 |
| `1P` | 4 | 1.252 | 1.255 |    +3.6 |
| `2P` | 4 | 1.807 | 1.807 |    +0.4 |
| `1S` | 2 | 0.615 | 0.594 |   -21.1 |
| `2S` | 2 | 1.412 | 1.413 |    +0.8 |
| `3S` | 2 | 1.970 | 1.967 |    -2.9 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
