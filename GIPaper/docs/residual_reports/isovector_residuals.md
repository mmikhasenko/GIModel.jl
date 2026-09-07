# Residuals: isovector (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.150 | 0.149 |    -1.5 | high |
| `2^1S_0` | 1.300 | 1.292 |    -8.4 | high |
| `1^3S_1` | 0.770 | 0.771 |    +1.2 | high |
| `2^3S_1` | 1.450 | 1.455 |    +5.1 | high |
| `1^3D_1` | 1.660 | 1.663 |    +2.8 | high |
| `3^1S_0` | 1.880 | 1.872 |    -8.0 | high |
| `3^3S_1` | 2.000 | 1.998 |    -2.2 | high |
| `2^3D_1` | 2.150 | 2.151 |    +0.6 | high |
| `1^1P_1` | 1.220 | 1.258 |   +37.6 | high |
| `2^1P_1` | 1.780 | 1.807 |   +27.1 | high |
| `1^3P_0` | 1.090 | 1.068 |   -21.8 | high |
| `2^3P_0` | 1.780 | 1.767 |   -12.5 | high |
| `1^3P_1` | 1.240 | 1.223 |   -16.8 | high |
| `2^3P_1` | 1.820 | 1.808 |   -11.8 | high |
| `1^3P_2` | 1.310 | 1.298 |   -12.3 | high |
| `2^3P_2` | 1.820 | 1.815 |    -5.3 | high |
| `1^3F_2` | 2.050 | 2.054 |    +3.9 | high |
| `1^1D_2` | 1.680 | 1.687 |    +7.0 | high |
| `2^1D_2` | 2.130 | 2.142 |   +12.1 | high |
| `1^3D_2` | 1.700 | 1.693 |    -6.7 | high |
| `2^3D_2` | 2.150 | 2.152 |    +2.2 | high |
| `1^3D_3` | 1.680 | 1.682 |    +1.5 | high |
| `2^3D_3` | 2.130 | 2.129 |    -1.1 | high |
| `1^3G_3` | 2.370 | 2.370 |    +0.1 | high |
| `1^1F_3` | 2.030 | 2.034 |    +4.4 | high |
| `1^3F_3` | 2.050 | 2.045 |    -4.6 | high |
| `1^3F_4` | 2.010 | 2.007 |    -2.8 | high |
| `1^1G_4` | 2.330 | 2.334 |    +3.6 | high |
| `1^3G_4` | 2.340 | 2.344 |    +4.2 | high |
| `1^3G_5` | 2.300 | 2.296 |    -3.7 | high |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.764 |  -615.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -615.1 | 0.149 |
| `2^1S_0` | 1.397 |  -105.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -105.0 | 1.292 |
| `1^3S_1` | 0.667 |  +104.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +104.1 | 0.771 |
| `2^3S_1` | 1.411 |   +44.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +44.3 | 1.455 |
| `1^3D_1` | 1.701 |    +0.0 |  -128.3 |   +98.3 |   -30.0 |    -9.4 |    +0.0 |   -39.5 | 1.663 |
| `3^1S_0` | 1.948 |   -76.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -76.3 | 1.872 |
| `3^3S_1` | 1.962 |   +36.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +36.0 | 1.998 |
| `2^3D_1` | 2.145 |    +0.0 |   -75.7 |   +85.8 |   +10.1 |    -5.1 |    +0.0 |    +5.0 | 2.151 |
| `1^1P_1` | 1.258 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.258 |
| `2^1P_1` | 1.807 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.807 |
| `1^3P_0` | 1.290 |    +0.0 |  -273.0 |  +105.2 |  -167.8 |   -54.5 |    +0.0 |  -222.2 | 1.068 |
| `2^3P_0` | 1.808 |    +0.0 |  -111.3 |   +91.9 |   -19.4 |   -21.1 |    +0.0 |   -40.4 | 1.767 |
| `1^3P_1` | 1.264 |    +0.0 |  -123.4 |   +57.0 |   -66.4 |   +25.7 |    +0.0 |   -40.7 | 1.223 |
| `2^3P_1` | 1.807 |    +0.0 |   -54.7 |   +45.6 |    -9.2 |   +10.4 |    +0.0 |    +1.2 | 1.808 |
| `1^3P_2` | 1.265 |    +0.0 |  +101.8 |   -64.7 |   +37.1 |    -4.5 |    +0.0 |   +32.5 | 1.298 |
| `2^3P_2` | 1.807 |    +0.0 |   +54.5 |   -44.5 |    +9.9 |    -2.1 |    +0.0 |    +7.8 | 1.815 |
| `1^3F_2` | 2.041 |    +0.0 |   -75.1 |   +90.8 |   +15.7 |    -3.4 |    +0.0 |   +12.3 | 2.054 |
| `1^1D_2` | 1.687 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.687 |
| `2^1D_2` | 2.142 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.142 |
| `1^3D_2` | 1.688 |    +0.0 |   -38.9 |   +35.1 |    -3.8 |    +8.8 |    +0.0 |    +5.1 | 1.693 |
| `2^3D_2` | 2.142 |    +0.0 |   -24.3 |   +29.1 |    +4.8 |    +5.0 |    +0.0 |    +9.8 | 2.152 |
| `1^3D_3` | 1.692 |    +0.0 |   +68.3 |   -76.9 |    -8.6 |    -2.3 |    +0.0 |   -10.9 | 1.682 |
| `2^3D_3` | 2.143 |    +0.0 |   +46.4 |   -59.4 |   -13.0 |    -1.4 |    +0.0 |   -14.3 | 2.129 |
| `1^3G_3` | 2.338 |    +0.0 |   -50.2 |   +84.0 |   +33.8 |    -1.7 |    +0.0 |   +32.1 | 2.370 |
| `1^1F_3` | 2.034 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.034 |
| `1^3F_3` | 2.035 |    +0.0 |   -17.4 |   +24.0 |    +6.6 |    +4.1 |    +0.0 |   +10.7 | 2.045 |
| `1^3F_4` | 2.038 |    +0.0 |   +47.7 |   -77.3 |   -29.6 |    -1.3 |    +0.0 |   -30.9 | 2.007 |
| `1^1G_4` | 2.334 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.334 |
| `1^3G_4` | 2.334 |    +0.0 |    -9.5 |   +17.6 |    +8.2 |    +2.2 |    +0.0 |   +10.4 | 2.344 |
| `1^3G_5` | 2.336 |    +0.0 |   +35.4 |   -74.6 |   -39.3 |    -0.8 |    +0.0 |   -40.0 | 2.296 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 1.456 | 1.455 |    +9.5 |  -0.999 |  +0.046 |
| `1^3D_1` | 1.661 | 1.663 |    +9.5 |  +0.046 |  +0.998 |
| `3^3S_1` | 1.998 | 1.998 |    +2.0 |  +1.000 |  -0.013 |
| `2^3D_1` | 2.150 | 2.151 |    +2.0 |  -0.013 |  -1.000 |
| `2^3P_2` | 1.815 | 1.815 |    +1.5 |  -1.000 |  +0.006 |
| `1^3F_2` | 2.054 | 2.054 |    +1.5 |  +0.006 |  +1.000 |
| `2^3D_3` | 2.129 | 2.129 |    +0.3 |  -1.000 |  +0.001 |
| `1^3G_3` | 2.370 | 2.370 |    +0.3 |  +0.001 |  +1.000 |

Mean absolute residual: 7.8 MeV.
Max absolute residual: 37.6 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.682 | 1.683 |    +1.0 |
| `2D` | 4 | 2.138 | 2.141 |    +3.3 |
| `1F` | 4 | 2.032 | 2.032 |    -0.3 |
| `1G` | 4 | 2.331 | 2.332 |    +0.8 |
| `1P` | 4 | 1.252 | 1.250 |    -1.7 |
| `2P` | 4 | 1.807 | 1.807 |    +0.6 |
| `1S` | 2 | 0.615 | 0.616 |    +0.5 |
| `2S` | 2 | 1.412 | 1.414 |    +1.7 |
| `3S` | 2 | 1.970 | 1.966 |    -3.6 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
