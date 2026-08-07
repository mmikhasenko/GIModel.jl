# Residuals: charmonium (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 2.970 | 2.958 |   -11.5 | high |
| `2^1S_0` | 3.620 | 3.619 |    -1.3 | high |
| `3^1S_0` | 4.060 | 4.058 |    -2.5 | high |
| `1^3S_1` | 3.100 | 3.091 |    -9.1 | high |
| `2^3S_1` | 3.680 | 3.680 |    +0.1 | high |
| `1^3D_1` | 3.820 | 3.824 |    +4.3 | high |
| `3^3S_1` | 4.100 | 4.101 |    +1.1 | high |
| `2^3D_1` | 4.190 | 4.198 |    +8.3 | high |
| `4^3S_1` | 4.450 | 4.450 |    -0.3 | medium |
| `3^3D_1` | 4.520 | 4.523 |    +2.9 | medium |
| `1^1P_1` | 3.520 | 3.523 |    +3.0 | high |
| `2^1P_1` | 3.960 | 3.962 |    +2.3 | high |
| `1^3P_0` | 3.440 | 3.458 |   +18.2 | high |
| `2^3P_0` | 3.920 | 3.924 |    +3.8 | high |
| `1^3P_1` | 3.510 | 3.526 |   +16.1 | high |
| `2^3P_1` | 3.950 | 3.964 |   +13.6 | high |
| `1^3P_2` | 3.550 | 3.533 |   -17.4 | high |
| `2^3P_2` | 3.980 | 3.969 |   -11.2 | high |
| `1^3F_2` | 4.090 | 4.091 |    +1.0 | medium |
| `1^1D_2` | 3.840 | 3.838 |    -1.6 | high |
| `2^1D_2` | 4.210 | 4.209 |    -0.9 | high |
| `1^3D_2` | 3.840 | 3.842 |    +1.7 | high |
| `2^3D_2` | 4.210 | 4.211 |    +1.3 | high |
| `1^3D_3` | 3.850 | 3.842 |    -8.1 | high |
| `2^3D_3` | 4.220 | 4.212 |    -8.0 | high |
| `1^1F_3` | 4.090 | 4.094 |    +4.4 | medium |
| `1^3F_3` | 4.100 | 4.097 |    -2.8 | medium |
| `1^3F_4` | 4.090 | 4.094 |    +4.1 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 3.087 |  -129.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -129.0 | 2.958 |
| `2^1S_0` | 3.668 |   -49.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -49.2 | 3.619 |
| `3^1S_0` | 4.091 |   -33.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -33.2 | 4.058 |
| `1^3S_1` | 3.066 |   +24.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +24.7 | 3.091 |
| `2^3S_1` | 3.666 |   +13.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.9 | 3.680 |
| `1^3D_1` | 3.839 |    +0.0 |   -26.1 |   +17.1 |    -9.0 |    -6.3 |    +0.0 |   -15.3 | 3.824 |
| `3^3S_1` | 4.091 |   +10.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +10.4 | 4.101 |
| `2^3D_1` | 4.209 |    +0.0 |   -19.5 |   +12.5 |    -6.9 |    -4.5 |    +0.0 |   -11.4 | 4.198 |
| `4^3S_1` | 4.441 |    +8.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.6 | 4.450 |
| `3^3D_1` | 4.532 |    +0.0 |   -16.0 |   +10.0 |    -6.0 |    -3.6 |    +0.0 |    -9.6 | 4.523 |
| `1^1P_1` | 3.523 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 3.523 |
| `2^1P_1` | 3.962 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 3.962 |
| `1^3P_0` | 3.530 |    +0.0 |   -57.2 |   +23.9 |   -33.4 |   -38.1 |    +0.0 |   -71.4 | 3.458 |
| `2^3P_0` | 3.963 |    +0.0 |   -33.1 |   +14.2 |   -18.9 |   -20.8 |    +0.0 |   -39.7 | 3.924 |
| `1^3P_1` | 3.523 |    +0.0 |   -24.1 |   +10.8 |   -13.3 |   +16.4 |    +0.0 |    +3.1 | 3.526 |
| `2^3P_1` | 3.962 |    +0.0 |   -15.7 |    +6.9 |    -8.8 |   +10.1 |    +0.0 |    +1.3 | 3.964 |
| `1^3P_2` | 3.523 |    +0.0 |   +22.9 |   -10.5 |   +12.4 |    -3.1 |    +0.0 |    +9.2 | 3.533 |
| `2^3P_2` | 3.962 |    +0.0 |   +15.2 |    -6.8 |    +8.4 |    -2.0 |    +0.0 |    +6.4 | 3.969 |
| `1^3F_2` | 4.095 |    +0.0 |   -16.3 |   +15.1 |    -1.3 |    -2.5 |    +0.0 |    -3.7 | 4.091 |
| `1^1D_2` | 3.838 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 3.838 |
| `2^1D_2` | 4.209 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 4.209 |
| `1^3D_2` | 3.838 |    +0.0 |    -8.2 |    +5.5 |    -2.6 |    +5.9 |    +0.0 |    +3.3 | 3.842 |
| `2^3D_2` | 4.209 |    +0.0 |    -6.2 |    +4.1 |    -2.1 |    +4.4 |    +0.0 |    +2.2 | 4.211 |
| `1^3D_3` | 3.839 |    +0.0 |   +15.9 |   -10.9 |    +5.0 |    -1.7 |    +0.0 |    +3.3 | 3.842 |
| `2^3D_3` | 4.209 |    +0.0 |   +12.2 |    -8.2 |    +4.0 |    -1.2 |    +0.0 |    +2.8 | 4.212 |
| `1^1F_3` | 4.094 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 4.094 |
| `1^3F_3` | 4.094 |    +0.0 |    -4.0 |    +3.7 |    -0.2 |    +3.0 |    +0.0 |    +2.8 | 4.097 |
| `1^3F_4` | 4.094 |    +0.0 |   +11.7 |   -11.1 |    +0.6 |    -1.0 |    +0.0 |    -0.4 | 4.094 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 3.680 | 3.680 |    +2.1 |  -1.000 |  +0.015 |
| `1^3D_1` | 3.824 | 3.824 |    +2.1 |  +0.015 |  +1.000 |
| `3^3S_1` | 4.101 | 4.101 |    +1.3 |  +1.000 |  -0.013 |
| `2^3D_1` | 4.198 | 4.198 |    +1.3 |  -0.013 |  -1.000 |
| `4^3S_1` | 4.450 | 4.450 |    +1.0 |  -1.000 |  +0.013 |
| `3^3D_1` | 4.523 | 4.523 |    +1.0 |  +0.013 |  +1.000 |
| `2^3P_2` | 3.969 | 3.969 |    +0.4 |  -1.000 |  +0.004 |
| `1^3F_2` | 4.091 | 4.091 |    +0.4 |  +0.004 |  +1.000 |

Mean absolute residual: 5.8 MeV.
Max absolute residual: 18.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 3.841 | 3.838 |    -2.2 |
| `2D` | 4 | 4.211 | 4.209 |    -1.5 |
| `3D` | 1 | 4.520 | 4.523 |    +2.9 |
| `1F` | 4 | 4.093 | 4.094 |    +1.9 |
| `1P` | 4 | 3.523 | 3.522 |    -1.0 |
| `2P` | 4 | 3.962 | 3.962 |    -0.3 |
| `1S` | 2 | 3.068 | 3.058 |    -9.7 |
| `2S` | 2 | 3.665 | 3.665 |    -0.2 |
| `3S` | 2 | 4.090 | 4.090 |    +0.2 |
| `4S` | 1 | 4.450 | 4.450 |    -0.3 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
