# Residuals: charmonium (GI-style)

Model: with GI momentum-sandwiched smeared contact hyperfine (every L), fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `OscillatorSolver` — nbasis = 24:8:80 adaptive, ΔE ≤ 0.1 MeV, beta in [0.25, 2.35] GeV (22 bracket points, refined to 0.002 GeV), 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 2.970 | 2.967 |    -2.9 | high |
| `2^1S_0` | 3.620 | 3.626 |    +5.8 | high |
| `3^1S_0` | 4.060 | 4.064 |    +3.6 | high |
| `1^3S_1` | 3.100 | 3.091 |    -8.6 | high |
| `2^3S_1` | 3.680 | 3.679 |    -0.9 | high |
| `1^3D_1` | 3.820 | 3.818 |    -1.7 | high |
| `3^3S_1` | 4.100 | 4.100 |    +0.0 | high |
| `2^3D_1` | 4.190 | 4.194 |    +4.1 | high |
| `4^3S_1` | 4.450 | 4.449 |    -1.2 | medium |
| `3^3D_1` | 4.520 | 4.520 |    -0.5 | medium |
| `1^1P_1` | 3.520 | 3.515 |    -4.9 | high |
| `2^1P_1` | 3.960 | 3.956 |    -3.8 | high |
| `1^3P_0` | 3.440 | 3.443 |    +2.9 | high |
| `2^3P_0` | 3.920 | 3.916 |    -3.5 | high |
| `1^3P_1` | 3.510 | 3.508 |    -1.8 | high |
| `2^3P_1` | 3.950 | 3.953 |    +3.2 | high |
| `1^3P_2` | 3.550 | 3.548 |    -1.8 | high |
| `2^3P_2` | 3.980 | 3.980 |    -0.5 | high |
| `1^3F_2` | 4.090 | 4.091 |    +1.4 | medium |
| `1^1D_2` | 3.840 | 3.837 |    -3.4 | high |
| `2^1D_2` | 4.210 | 4.208 |    -2.5 | high |
| `1^3D_2` | 3.840 | 3.837 |    -2.6 | high |
| `2^3D_2` | 4.210 | 4.208 |    -1.8 | high |
| `1^3D_3` | 3.850 | 3.848 |    -2.3 | high |
| `2^3D_3` | 4.220 | 4.217 |    -3.4 | high |
| `1^1F_3` | 4.090 | 4.094 |    +3.8 | medium |
| `1^3F_3` | 4.100 | 4.096 |    -3.9 | medium |
| `1^3F_4` | 4.090 | 4.095 |    +4.8 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 3.078 |  -111.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -111.4 | 2.967 |
| `2^1S_0` | 3.667 |   -41.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -41.0 | 3.626 |
| `3^1S_0` | 4.091 |   -27.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -27.1 | 4.064 |
| `1^3S_1` | 3.066 |   +25.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +25.6 | 3.091 |
| `2^3S_1` | 3.666 |   +12.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +12.7 | 3.679 |
| `1^3D_1` | 3.841 |    +0.7 |   -42.4 |   +23.2 |   -19.3 |    -3.9 |    +0.0 |   -22.4 | 3.818 |
| `3^3S_1` | 4.091 |    +9.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +9.0 | 4.100 |
| `2^3D_1` | 4.210 |    +0.6 |   -31.2 |   +17.2 |   -14.0 |    -2.7 |    +0.0 |   -16.2 | 4.194 |
| `4^3S_1` | 4.442 |    +7.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +7.2 | 4.449 |
| `3^3D_1` | 4.533 |    +0.5 |   -25.5 |   +13.8 |   -11.7 |    -2.2 |    +0.0 |   -13.3 | 4.520 |
| `1^1P_1` | 3.523 |    -8.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -8.2 | 3.515 |
| `2^1P_1` | 3.963 |    -6.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -6.4 | 3.956 |
| `1^3P_0` | 3.535 |    +3.7 |   -94.6 |   +22.8 |   -71.8 |   -23.8 |    +0.0 |   -91.9 | 3.443 |
| `2^3P_0` | 3.964 |    +2.5 |   -52.4 |   +14.7 |   -37.8 |   -12.4 |    +0.0 |   -47.8 | 3.916 |
| `1^3P_1` | 3.524 |    +2.8 |   -39.8 |   +11.1 |   -28.7 |   +10.2 |    +0.0 |   -15.6 | 3.508 |
| `2^3P_1` | 3.963 |    +2.2 |   -24.9 |    +7.3 |   -17.6 |    +6.0 |    +0.0 |    -9.4 | 3.953 |
| `1^3P_2` | 3.525 |    +2.2 |   +33.6 |   -10.7 |   +22.9 |    -1.8 |    +0.0 |   +23.3 | 3.548 |
| `2^3P_2` | 3.963 |    +1.8 |   +22.9 |    -7.2 |   +15.7 |    -1.1 |    +0.0 |   +16.3 | 3.980 |
| `1^3F_2` | 4.095 |    +0.2 |   -26.1 |   +23.6 |    -2.5 |    -1.5 |    +0.0 |    -3.7 | 4.091 |
| `1^1D_2` | 3.838 |    -1.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -1.9 | 3.837 |
| `2^1D_2` | 4.209 |    -1.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -1.6 | 4.208 |
| `1^3D_2` | 3.839 |    +0.6 |   -13.0 |    +7.6 |    -5.4 |    +3.6 |    +0.0 |    -1.2 | 3.837 |
| `2^3D_2` | 4.209 |    +0.5 |    -9.8 |    +5.7 |    -4.1 |    +2.6 |    +0.0 |    -1.0 | 4.208 |
| `1^3D_3` | 3.839 |    +0.6 |   +24.0 |   -15.0 |    +9.0 |    -1.0 |    +0.0 |    +8.6 | 3.848 |
| `2^3D_3` | 4.210 |    +0.5 |   +18.6 |   -11.3 |    +7.3 |    -0.7 |    +0.0 |    +7.0 | 4.217 |
| `1^1F_3` | 4.094 |    -0.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -0.7 | 4.094 |
| `1^3F_3` | 4.094 |    +0.2 |    -6.2 |    +5.9 |    -0.4 |    +1.8 |    +0.0 |    +1.6 | 4.096 |
| `1^3F_4` | 4.095 |    +0.2 |   +17.8 |   -17.4 |    +0.4 |    -0.6 |    +0.0 |    +0.0 | 4.095 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 3.679 | 3.679 |    +1.1 |  +1.000 |  -0.008 |
| `1^3D_1` | 3.818 | 3.818 |    +1.1 |  +0.008 |  +1.000 |
| `3^3S_1` | 4.100 | 4.100 |    +0.7 |  +1.000 |  -0.007 |
| `2^3D_1` | 4.194 | 4.194 |    +0.7 |  +0.007 |  +1.000 |
| `4^3S_1` | 4.449 | 4.449 |    +0.5 |  +1.000 |  -0.007 |
| `3^3D_1` | 4.519 | 4.520 |    +0.5 |  +0.007 |  +1.000 |
| `2^3P_2` | 3.980 | 3.980 |    +0.1 |  +1.000 |  -0.001 |
| `1^3F_2` | 4.091 | 4.091 |    +0.1 |  +0.001 |  +1.000 |

Mean absolute residual: 2.9 MeV.
Max absolute residual: 8.6 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 3.841 | 3.838 |    -2.6 |
| `2D` | 4 | 4.211 | 4.209 |    -1.6 |
| `3D` | 1 | 4.520 | 4.520 |    -0.5 |
| `1F` | 4 | 4.093 | 4.094 |    +1.7 |
| `1P` | 4 | 3.523 | 3.521 |    -2.2 |
| `2P` | 4 | 3.962 | 3.962 |    -0.6 |
| `1S` | 2 | 3.068 | 3.060 |    -7.2 |
| `2S` | 2 | 3.665 | 3.666 |    +0.8 |
| `3S` | 2 | 4.090 | 4.091 |    +0.9 |
| `4S` | 1 | 4.450 | 4.449 |    -1.2 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
