# Residuals: charmonium (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 2.970 | 2.967 |    -3.3 | high |
| `2^1S_0` | 3.620 | 3.625 |    +5.4 | high |
| `3^1S_0` | 4.060 | 4.063 |    +3.1 | high |
| `1^3S_1` | 3.100 | 3.091 |    -8.8 | high |
| `2^3S_1` | 3.680 | 3.679 |    -1.2 | high |
| `1^3D_1` | 3.820 | 3.818 |    -2.5 | high |
| `3^3S_1` | 4.100 | 4.100 |    -0.3 | high |
| `2^3D_1` | 4.190 | 4.193 |    +3.5 | high |
| `4^3S_1` | 4.450 | 4.448 |    -1.7 | medium |
| `3^3D_1` | 4.520 | 4.519 |    -1.2 | medium |
| `1^1P_1` | 3.520 | 3.523 |    +3.0 | high |
| `2^1P_1` | 3.960 | 3.962 |    +2.3 | high |
| `1^3P_0` | 3.440 | 3.439 |    -1.0 | high |
| `2^3P_0` | 3.920 | 3.914 |    -6.2 | high |
| `1^3P_1` | 3.510 | 3.505 |    -4.7 | high |
| `2^3P_1` | 3.950 | 3.951 |    +0.9 | high |
| `1^3P_2` | 3.550 | 3.546 |    -4.1 | high |
| `2^3P_2` | 3.980 | 3.978 |    -2.4 | high |
| `1^3F_2` | 4.090 | 4.091 |    +1.1 | medium |
| `1^1D_2` | 3.840 | 3.838 |    -1.6 | high |
| `2^1D_2` | 4.210 | 4.209 |    -0.9 | high |
| `1^3D_2` | 3.840 | 3.837 |    -3.3 | high |
| `2^3D_2` | 4.210 | 4.208 |    -2.4 | high |
| `1^3D_3` | 3.850 | 3.847 |    -2.9 | high |
| `2^3D_3` | 4.220 | 4.216 |    -4.0 | high |
| `1^1F_3` | 4.090 | 4.094 |    +4.4 | medium |
| `1^3F_3` | 4.100 | 4.096 |    -4.1 | medium |
| `1^3F_4` | 4.090 | 4.095 |    +4.5 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 3.078 |  -111.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -111.7 | 2.967 |
| `2^1S_0` | 3.667 |   -41.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -41.1 | 3.625 |
| `3^1S_0` | 4.090 |   -27.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -27.2 | 4.063 |
| `1^3S_1` | 3.066 |   +25.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +25.6 | 3.091 |
| `2^3S_1` | 3.666 |   +12.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +12.7 | 3.679 |
| `1^3D_1` | 3.841 |    +0.0 |   -42.5 |   +23.2 |   -19.4 |    -3.9 |    +0.0 |   -23.2 | 3.818 |
| `3^3S_1` | 4.091 |    +9.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +9.1 | 4.100 |
| `2^3D_1` | 4.210 |    +0.0 |   -31.3 |   +17.2 |   -14.1 |    -2.8 |    +0.0 |   -16.8 | 4.193 |
| `4^3S_1` | 4.441 |    +7.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +7.3 | 4.448 |
| `3^3D_1` | 4.533 |    +0.0 |   -25.6 |   +13.8 |   -11.8 |    -2.2 |    +0.0 |   -13.9 | 4.519 |
| `1^1P_1` | 3.523 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 3.523 |
| `2^1P_1` | 3.962 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 3.962 |
| `1^3P_0` | 3.536 |    +0.0 |   -95.9 |   +22.8 |   -73.0 |   -24.1 |    +0.0 |   -97.1 | 3.439 |
| `2^3P_0` | 3.965 |    +0.0 |   -52.9 |   +14.7 |   -38.2 |   -12.5 |    +0.0 |   -50.7 | 3.914 |
| `1^3P_1` | 3.524 |    +0.0 |   -40.3 |   +11.1 |   -29.2 |   +10.3 |    +0.0 |   -18.8 | 3.505 |
| `2^3P_1` | 3.963 |    +0.0 |   -25.2 |    +7.3 |   -17.9 |    +6.1 |    +0.0 |   -11.8 | 3.951 |
| `1^3P_2` | 3.524 |    +0.0 |   +34.0 |   -10.7 |   +23.2 |    -1.8 |    +0.0 |   +21.4 | 3.546 |
| `2^3P_2` | 3.963 |    +0.0 |   +23.1 |    -7.2 |   +15.9 |    -1.1 |    +0.0 |   +14.7 | 3.978 |
| `1^3F_2` | 4.095 |    +0.0 |   -26.1 |   +23.6 |    -2.5 |    -1.5 |    +0.0 |    -4.0 | 4.091 |
| `1^1D_2` | 3.838 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 3.838 |
| `2^1D_2` | 4.209 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 4.209 |
| `1^3D_2` | 3.839 |    +0.0 |   -13.0 |    +7.6 |    -5.4 |    +3.6 |    +0.0 |    -1.8 | 3.837 |
| `2^3D_2` | 4.209 |    +0.0 |    -9.9 |    +5.7 |    -4.2 |    +2.6 |    +0.0 |    -1.5 | 4.208 |
| `1^3D_3` | 3.839 |    +0.0 |   +24.1 |   -15.0 |    +9.1 |    -1.0 |    +0.0 |    +8.1 | 3.847 |
| `2^3D_3` | 4.209 |    +0.0 |   +18.6 |   -11.3 |    +7.3 |    -0.7 |    +0.0 |    +6.6 | 4.216 |
| `1^1F_3` | 4.094 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 4.094 |
| `1^3F_3` | 4.094 |    +0.0 |    -6.2 |    +5.9 |    -0.4 |    +1.8 |    +0.0 |    +1.4 | 4.096 |
| `1^3F_4` | 4.095 |    +0.0 |   +17.9 |   -17.5 |    +0.4 |    -0.6 |    +0.0 |    -0.2 | 4.095 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 3.679 | 3.679 |    +1.1 |  -1.000 |  +0.008 |
| `1^3D_1` | 3.817 | 3.818 |    +1.1 |  +0.008 |  +1.000 |
| `3^3S_1` | 4.100 | 4.100 |    +0.7 |  +1.000 |  -0.007 |
| `2^3D_1` | 4.193 | 4.193 |    +0.7 |  -0.007 |  -1.000 |
| `4^3S_1` | 4.448 | 4.448 |    +0.5 |  -1.000 |  +0.007 |
| `3^3D_1` | 4.519 | 4.519 |    +0.5 |  +0.007 |  +1.000 |
| `2^3P_2` | 3.978 | 3.978 |    +0.1 |  -1.000 |  +0.001 |
| `1^3F_2` | 4.091 | 4.091 |    +0.1 |  +0.001 |  +1.000 |

Mean absolute residual: 3.0 MeV.
Max absolute residual: 8.8 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 3.841 | 3.838 |    -2.6 |
| `2D` | 4 | 4.211 | 4.209 |    -1.7 |
| `3D` | 1 | 4.520 | 4.519 |    -1.2 |
| `1F` | 4 | 4.093 | 4.094 |    +1.7 |
| `1P` | 4 | 3.523 | 3.521 |    -2.2 |
| `2P` | 4 | 3.962 | 3.962 |    -0.7 |
| `1S` | 2 | 3.068 | 3.060 |    -7.4 |
| `2S` | 2 | 3.665 | 3.665 |    +0.5 |
| `3S` | 2 | 4.090 | 4.091 |    +0.5 |
| `4S` | 1 | 4.450 | 4.448 |    -1.7 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
