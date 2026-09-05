# Residuals: bottomonium (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 9.400 | 9.390 |   -10.4 | high |
| `2^1S_0` | 9.980 | 9.971 |    -9.2 | high |
| `3^1S_0` | 10.340 | 10.329 |   -11.2 | high |
| `1^3S_1` | 9.460 | 9.458 |    -2.1 | high |
| `2^3S_1` | 10.000 | 10.004 |    +3.9 | high |
| `1^3D_1` | 10.140 | 10.141 |    +1.3 | high |
| `3^3S_1` | 10.350 | 10.353 |    +3.3 | high |
| `2^3D_1` | 10.440 | 10.444 |    +4.1 | high |
| `4^3S_1` | 10.630 | 10.632 |    +2.4 | medium |
| `3^3D_1` | 10.700 | 10.701 |    +1.1 | medium |
| `5^3S_1` | 10.880 | 10.874 |    -6.0 | medium |
| `6^3S_1` | 11.100 | 11.092 |    -8.4 | medium |
| `1^1P_1` | 9.880 | 9.883 |    +3.0 | high |
| `2^1P_1` | 10.250 | 10.251 |    +1.3 | high |
| `1^3P_0` | 9.850 | 9.853 |    +3.4 | high |
| `2^3P_0` | 10.230 | 10.231 |    +0.8 | high |
| `1^3P_1` | 9.880 | 9.884 |    +3.7 | high |
| `2^3P_1` | 10.250 | 10.252 |    +1.7 | high |
| `1^3P_2` | 9.900 | 9.888 |   -12.2 | high |
| `2^3P_2` | 10.260 | 10.255 |    -5.2 | high |
| `1^3F_2` | 10.350 | 10.352 |    +1.6 | medium |
| `1^1D_2` | 10.150 | 10.148 |    -1.6 | high |
| `2^1D_2` | 10.450 | 10.450 |    -0.0 | high |
| `1^3D_2` | 10.150 | 10.149 |    -0.9 | high |
| `2^3D_2` | 10.450 | 10.451 |    +0.5 | high |
| `1^3D_3` | 10.160 | 10.151 |    -9.1 | high |
| `2^3D_3` | 10.450 | 10.452 |    +2.0 | high |
| `1^1F_3` | 10.350 | 10.355 |    +4.6 | medium |
| `1^3F_3` | 10.350 | 10.355 |    +5.2 | medium |
| `1^3F_4` | 10.360 | 10.356 |    -4.2 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 9.456 |   -66.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -66.3 | 9.390 |
| `2^1S_0` | 10.000 |   -29.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -29.5 | 9.971 |
| `3^1S_0` | 10.350 |   -21.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -21.2 | 10.329 |
| `1^3S_1` | 9.445 |   +13.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.0 | 9.458 |
| `2^3S_1` | 9.997 |    +6.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +6.8 | 10.004 |
| `1^3D_1` | 10.149 |    +0.0 |    -9.6 |    +4.6 |    -5.0 |    -2.4 |    +0.0 |    -7.3 | 10.141 |
| `3^3S_1` | 10.348 |    +5.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +5.2 | 10.353 |
| `2^3D_1` | 10.450 |    +0.0 |    -7.8 |    +3.6 |    -4.1 |    -1.9 |    +0.0 |    -6.0 | 10.444 |
| `4^3S_1` | 10.628 |    +4.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +4.4 | 10.632 |
| `3^3D_1` | 10.706 |    +0.0 |    -6.7 |    +3.1 |    -3.7 |    -1.6 |    +0.0 |    -5.3 | 10.701 |
| `5^3S_1` | 10.870 |    +4.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +4.0 | 10.874 |
| `6^3S_1` | 11.088 |    +3.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +3.6 | 11.092 |
| `1^1P_1` | 9.883 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 9.883 |
| `2^1P_1` | 10.251 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 10.251 |
| `1^3P_0` | 9.886 |    +0.0 |   -23.4 |    +8.0 |   -15.5 |   -16.8 |    +0.0 |   -32.3 | 9.853 |
| `2^3P_0` | 10.253 |    +0.0 |   -16.0 |    +5.4 |   -10.7 |   -11.3 |    +0.0 |   -22.0 | 10.231 |
| `1^3P_1` | 9.883 |    +0.0 |    -9.9 |    +3.5 |    -6.4 |    +7.2 |    +0.0 |    +0.8 | 9.884 |
| `2^3P_1` | 10.251 |    +0.0 |    -7.1 |    +2.5 |    -4.6 |    +5.0 |    +0.0 |    +0.4 | 10.252 |
| `1^3P_2` | 9.883 |    +0.0 |    +9.6 |    -3.4 |    +6.2 |    -1.4 |    +0.0 |    +4.8 | 9.888 |
| `2^3P_2` | 10.251 |    +0.0 |    +6.9 |    -2.4 |    +4.5 |    -1.0 |    +0.0 |    +3.5 | 10.255 |
| `1^3F_2` | 10.355 |    +0.0 |    -6.0 |    +3.8 |    -2.2 |    -0.9 |    +0.0 |    -3.1 | 10.352 |
| `1^1D_2` | 10.148 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 10.148 |
| `2^1D_2` | 10.450 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 10.450 |
| `1^3D_2` | 10.148 |    +0.0 |    -3.1 |    +1.5 |    -1.6 |    +2.3 |    +0.0 |    +0.7 | 10.149 |
| `2^3D_2` | 10.450 |    +0.0 |    -2.5 |    +1.2 |    -1.3 |    +1.8 |    +0.0 |    +0.5 | 10.451 |
| `1^3D_3` | 10.148 |    +0.0 |    +6.0 |    -3.0 |    +3.1 |    -0.6 |    +0.0 |    +2.4 | 10.151 |
| `2^3D_3` | 10.450 |    +0.0 |    +4.9 |    -2.4 |    +2.6 |    -0.5 |    +0.0 |    +2.0 | 10.452 |
| `1^1F_3` | 10.355 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 10.355 |
| `1^3F_3` | 10.355 |    +0.0 |    -1.5 |    +0.9 |    -0.5 |    +1.1 |    +0.0 |    +0.6 | 10.355 |
| `1^3F_4` | 10.355 |    +0.0 |    +4.4 |    -2.8 |    +1.6 |    -0.4 |    +0.0 |    +1.2 | 10.356 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 10.004 | 10.004 |    +0.4 |  +1.000 |  -0.003 |
| `1^3D_1` | 10.141 | 10.141 |    +0.4 |  +0.003 |  +1.000 |
| `3^3S_1` | 10.353 | 10.353 |    +0.1 |  -1.000 |  +0.002 |
| `2^3D_1` | 10.444 | 10.444 |    +0.1 |  -0.002 |  -1.000 |
| `4^3S_1` | 10.632 | 10.632 |    -0.0 |  -1.000 |  -0.000 |
| `3^3D_1` | 10.701 | 10.701 |    -0.0 |  -0.000 |  +1.000 |
| `2^3P_2` | 10.255 | 10.255 |    +0.0 |  -1.000 |  +0.000 |
| `1^3F_2` | 10.352 | 10.352 |    +0.0 |  +0.000 |  +1.000 |

Mean absolute residual: 4.1 MeV.
Max absolute residual: 12.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 10.152 | 10.148 |    -3.6 |
| `2D` | 4 | 10.448 | 10.450 |    +1.4 |
| `3D` | 1 | 10.700 | 10.701 |    +1.1 |
| `1F` | 4 | 10.353 | 10.355 |    +1.4 |
| `1P` | 4 | 9.886 | 9.883 |    -3.1 |
| `2P` | 4 | 10.252 | 10.251 |    -1.3 |
| `1S` | 2 | 9.445 | 9.441 |    -4.2 |
| `2S` | 2 | 9.995 | 9.996 |    +0.6 |
| `3S` | 2 | 10.348 | 10.347 |    -0.4 |
| `4S` | 1 | 10.630 | 10.632 |    +2.4 |
| `5S` | 1 | 10.880 | 10.874 |    -6.0 |
| `6S` | 1 | 11.100 | 11.092 |    -8.4 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
