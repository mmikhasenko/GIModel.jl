# Residuals: bottomonium (GI-style)

Model: with GI momentum-sandwiched smeared contact hyperfine (every L), fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 9.400 | 9.392 |    -7.7 | high |
| `2^1S_0` | 9.980 | 9.974 |    -6.3 | high |
| `3^1S_0` | 10.340 | 10.331 |    -8.6 | high |
| `1^3S_1` | 9.460 | 9.458 |    -1.7 | high |
| `2^3S_1` | 10.000 | 10.004 |    +3.5 | high |
| `1^3D_1` | 10.140 | 10.137 |    -2.9 | high |
| `3^3S_1` | 10.350 | 10.353 |    +2.8 | high |
| `2^3D_1` | 10.440 | 10.441 |    +0.7 | high |
| `4^3S_1` | 10.630 | 10.632 |    +1.9 | medium |
| `3^3D_1` | 10.700 | 10.698 |    -2.0 | medium |
| `5^3S_1` | 10.880 | 10.874 |    -6.5 | medium |
| `6^3S_1` | 11.100 | 11.091 |    -8.8 | medium |
| `1^1P_1` | 9.880 | 9.880 |    +0.4 | high |
| `2^1P_1` | 10.250 | 10.249 |    -0.7 | high |
| `1^3P_0` | 9.850 | 9.845 |    -5.4 | high |
| `2^3P_0` | 10.230 | 10.225 |    -5.2 | high |
| `1^3P_1` | 9.880 | 9.875 |    -5.4 | high |
| `2^3P_1` | 10.250 | 10.245 |    -4.7 | high |
| `1^3P_2` | 9.900 | 9.896 |    -4.4 | high |
| `2^3P_2` | 10.260 | 10.260 |    +0.4 | high |
| `1^3F_2` | 10.350 | 10.350 |    -0.3 | medium |
| `1^1D_2` | 10.150 | 10.148 |    -2.3 | high |
| `2^1D_2` | 10.450 | 10.449 |    -0.6 | high |
| `1^3D_2` | 10.150 | 10.147 |    -3.3 | high |
| `2^3D_2` | 10.450 | 10.449 |    -1.4 | high |
| `1^3D_3` | 10.160 | 10.155 |    -5.5 | high |
| `2^3D_3` | 10.450 | 10.455 |    +5.1 | high |
| `1^1F_3` | 10.350 | 10.354 |    +4.3 | medium |
| `1^3F_3` | 10.350 | 10.354 |    +4.3 | medium |
| `1^3F_4` | 10.360 | 10.358 |    -2.3 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 9.452 |   -59.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -59.5 | 9.392 |
| `2^1S_0` | 9.999 |   -25.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -25.0 | 9.974 |
| `3^1S_0` | 10.349 |   -17.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -17.6 | 10.331 |
| `1^3S_1` | 9.445 |   +13.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.7 | 9.458 |
| `2^3S_1` | 9.997 |    +6.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +6.6 | 10.004 |
| `1^3D_1` | 10.149 |    +0.3 |   -15.4 |    +4.8 |   -10.6 |    -1.4 |    +0.0 |   -11.8 | 10.137 |
| `3^3S_1` | 10.348 |    +4.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +4.8 | 10.353 |
| `2^3D_1` | 10.450 |    +0.2 |   -12.5 |    +3.7 |    -8.7 |    -1.2 |    +0.0 |    -9.7 | 10.441 |
| `4^3S_1` | 10.628 |    +4.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +4.0 | 10.632 |
| `3^3D_1` | 10.706 |    +0.2 |   -10.8 |    +3.1 |    -7.7 |    -1.0 |    +0.0 |    -8.5 | 10.698 |
| `5^3S_1` | 10.870 |    +3.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +3.5 | 10.874 |
| `6^3S_1` | 11.088 |    +3.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +3.2 | 11.091 |
| `1^1P_1` | 9.883 |    -2.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -2.6 | 9.880 |
| `2^1P_1` | 10.251 |    -2.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -2.0 | 10.249 |
| `1^3P_0` | 9.888 |    +1.1 |   -38.7 |    +4.8 |   -33.9 |   -10.6 |    +0.0 |   -43.3 | 9.845 |
| `2^3P_0` | 10.254 |    +0.8 |   -26.1 |    +3.3 |   -22.9 |    -7.0 |    +0.0 |   -29.0 | 10.225 |
| `1^3P_1` | 9.883 |    +0.9 |   -16.3 |    +2.3 |   -14.0 |    +4.5 |    +0.0 |    -8.6 | 9.875 |
| `2^3P_1` | 10.251 |    +0.7 |   -11.5 |    +1.6 |    -9.9 |    +3.1 |    +0.0 |    -6.1 | 10.245 |
| `1^3P_2` | 9.884 |    +0.8 |   +14.2 |    -2.2 |   +12.0 |    -0.8 |    +0.0 |   +12.0 | 9.896 |
| `2^3P_2` | 10.252 |    +0.6 |   +10.3 |    -1.6 |    +8.7 |    -0.6 |    +0.0 |    +8.7 | 10.260 |
| `1^3F_2` | 10.355 |    +0.1 |    -9.5 |    +4.9 |    -4.6 |    -0.5 |    +0.0 |    -5.1 | 10.350 |
| `1^1D_2` | 10.148 |    -0.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -0.7 | 10.148 |
| `2^1D_2` | 10.450 |    -0.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -0.6 | 10.449 |
| `1^3D_2` | 10.148 |    +0.2 |    -4.8 |    +1.6 |    -3.3 |    +1.4 |    +0.0 |    -1.7 | 10.147 |
| `2^3D_2` | 10.450 |    +0.2 |    -3.9 |    +1.2 |    -2.7 |    +1.1 |    +0.0 |    -1.4 | 10.449 |
| `1^3D_3` | 10.149 |    +0.2 |    +9.2 |    -3.1 |    +6.1 |    -0.4 |    +0.0 |    +6.0 | 10.155 |
| `2^3D_3` | 10.450 |    +0.2 |    +7.5 |    -2.4 |    +5.1 |    -0.3 |    +0.0 |    +5.0 | 10.455 |
| `1^1F_3` | 10.355 |    -0.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -0.3 | 10.354 |
| `1^3F_3` | 10.355 |    +0.1 |    -2.3 |    +1.2 |    -1.1 |    +0.7 |    +0.0 |    -0.3 | 10.354 |
| `1^3F_4` | 10.355 |    +0.1 |    +6.7 |    -3.6 |    +3.2 |    -0.2 |    +0.0 |    +3.0 | 10.358 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 10.004 | 10.004 |    +0.2 |  +1.000 |  -0.001 |
| `1^3D_1` | 10.137 | 10.137 |    +0.2 |  +0.001 |  +1.000 |
| `3^3S_1` | 10.353 | 10.353 |    +0.0 |  -1.000 |  +0.000 |
| `2^3D_1` | 10.441 | 10.441 |    +0.0 |  -0.000 |  -1.000 |
| `4^3S_1` | 10.632 | 10.632 |    -0.0 |  +1.000 |  +0.001 |
| `3^3D_1` | 10.698 | 10.698 |    -0.0 |  -0.001 |  +1.000 |
| `2^3P_2` | 10.260 | 10.260 |    -0.0 |  +1.000 |  +0.000 |
| `1^3F_2` | 10.350 | 10.350 |    -0.0 |  -0.000 |  +1.000 |

Mean absolute residual: 3.6 MeV.
Max absolute residual: 8.8 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 10.152 | 10.148 |    -3.7 |
| `2D` | 4 | 10.448 | 10.450 |    +1.4 |
| `3D` | 1 | 10.700 | 10.698 |    -2.0 |
| `1F` | 4 | 10.353 | 10.355 |    +1.4 |
| `1P` | 4 | 9.886 | 9.882 |    -3.5 |
| `2P` | 4 | 10.252 | 10.251 |    -1.6 |
| `1S` | 2 | 9.445 | 9.442 |    -3.2 |
| `2S` | 2 | 9.995 | 9.996 |    +1.1 |
| `3S` | 2 | 10.348 | 10.347 |    -0.1 |
| `4S` | 1 | 10.630 | 10.632 |    +1.9 |
| `5S` | 1 | 10.880 | 10.874 |    -6.5 |
| `6S` | 1 | 11.100 | 11.091 |    -8.8 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
