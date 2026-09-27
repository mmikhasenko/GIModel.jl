# Residuals: isovector (GI-style)

Model: with GI momentum-sandwiched smeared contact hyperfine (every L), fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.150 | 0.149 |    -1.5 | high |
| `2^1S_0` | 1.300 | 1.292 |    -8.4 | high |
| `1^3S_1` | 0.770 | 0.771 |    +1.2 | high |
| `2^3S_1` | 1.450 | 1.455 |    +5.1 | high |
| `1^3D_1` | 1.660 | 1.666 |    +5.9 | high |
| `3^1S_0` | 1.880 | 1.872 |    -8.0 | high |
| `3^3S_1` | 2.000 | 1.998 |    -2.2 | high |
| `2^3D_1` | 2.150 | 2.153 |    +3.4 | high |
| `1^1P_1` | 1.220 | 1.218 |    -1.6 | high |
| `2^1P_1` | 1.780 | 1.777 |    -3.2 | high |
| `1^3P_0` | 1.090 | 1.087 |    -3.2 | high |
| `2^3P_0` | 1.780 | 1.778 |    -2.0 | high |
| `1^3P_1` | 1.240 | 1.237 |    -2.5 | high |
| `2^3P_1` | 1.820 | 1.818 |    -2.1 | high |
| `1^3P_2` | 1.310 | 1.307 |    -3.5 | high |
| `2^3P_2` | 1.820 | 1.823 |    +3.1 | high |
| `1^3F_2` | 2.050 | 2.055 |    +4.6 | high |
| `1^1D_2` | 1.680 | 1.680 |    +0.4 | high |
| `2^1D_2` | 2.130 | 2.135 |    +4.8 | high |
| `1^3D_2` | 1.700 | 1.696 |    -4.3 | high |
| `2^3D_2` | 2.150 | 2.155 |    +4.7 | high |
| `1^3D_3` | 1.680 | 1.683 |    +3.1 | high |
| `2^3D_3` | 2.130 | 2.131 |    +1.0 | high |
| `1^3G_3` | 2.370 | 2.370 |    +0.3 | high |
| `1^1F_3` | 2.030 | 2.033 |    +2.8 | high |
| `1^3F_3` | 2.050 | 2.046 |    -4.0 | high |
| `1^3F_4` | 2.010 | 2.008 |    -2.4 | high |
| `1^1G_4` | 2.330 | 2.333 |    +3.2 | high |
| `1^3G_4` | 2.340 | 2.344 |    +4.3 | high |
| `1^3G_5` | 2.300 | 2.296 |    -3.6 | high |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.764 |  -615.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -615.1 | 0.149 |
| `2^1S_0` | 1.397 |  -105.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -105.0 | 1.292 |
| `1^3S_1` | 0.667 |  +104.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +104.1 | 0.771 |
| `2^3S_1` | 1.411 |   +44.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +44.3 | 1.455 |
| `1^3D_1` | 1.700 |    +3.1 |  -127.5 |   +98.6 |   -29.0 |    -9.4 |    +0.0 |   -35.3 | 1.666 |
| `3^1S_0` | 1.948 |   -76.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -76.3 | 1.872 |
| `3^3S_1` | 1.962 |   +36.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +36.0 | 1.998 |
| `2^3D_1` | 2.145 |    +2.8 |   -75.2 |   +85.8 |   +10.5 |    -5.1 |    +0.0 |    +8.2 | 2.153 |
| `1^1P_1` | 1.262 |   -43.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -43.2 | 1.218 |
| `2^1P_1` | 1.810 |   -32.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -32.7 | 1.777 |
| `1^3P_0` | 1.284 |   +18.1 |  -267.7 |  +106.6 |  -161.1 |   -53.9 |    +0.0 |  -196.9 | 1.087 |
| `2^3P_0` | 1.807 |   +10.3 |  -109.8 |   +91.6 |   -18.1 |   -20.9 |    +0.0 |   -28.7 | 1.778 |
| `1^3P_1` | 1.261 |   +13.9 |  -120.8 |   +57.6 |   -63.2 |   +25.4 |    +0.0 |   -23.9 | 1.237 |
| `2^3P_1` | 1.807 |    +9.4 |   -53.9 |   +45.4 |    -8.5 |   +10.3 |    +0.0 |   +11.2 | 1.818 |
| `1^3P_2` | 1.268 |    +8.6 |   +99.7 |   -65.2 |   +34.6 |    -4.5 |    +0.0 |   +38.7 | 1.307 |
| `2^3P_2` | 1.808 |    +8.2 |   +53.5 |   -44.4 |    +9.1 |    -2.1 |    +0.0 |   +15.3 | 1.823 |
| `1^3F_2` | 2.041 |    +0.7 |   -74.9 |   +90.9 |   +16.0 |    -3.4 |    +0.0 |   +13.3 | 2.055 |
| `1^1D_2` | 1.687 |    -6.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -6.8 | 1.680 |
| `2^1D_2` | 2.142 |    -7.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -7.5 | 2.135 |
| `1^3D_2` | 1.688 |    +2.4 |   -38.7 |   +35.2 |    -3.5 |    +8.8 |    +0.0 |    +7.7 | 1.696 |
| `2^3D_2` | 2.142 |    +2.5 |   -24.1 |   +29.1 |    +4.9 |    +5.0 |    +0.0 |   +12.4 | 2.155 |
| `1^3D_3` | 1.693 |    +1.6 |   +68.0 |   -77.0 |    -9.1 |    -2.3 |    +0.0 |    -9.8 | 1.683 |
| `2^3D_3` | 2.143 |    +2.1 |   +46.2 |   -59.4 |   -13.2 |    -1.4 |    +0.0 |   -12.5 | 2.131 |
| `1^3G_3` | 2.338 |    +0.2 |   -50.2 |   +84.1 |   +33.9 |    -1.7 |    +0.0 |   +32.4 | 2.370 |
| `1^1F_3` | 2.034 |    -1.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -1.6 | 2.033 |
| `1^3F_3` | 2.035 |    +0.6 |   -17.4 |   +24.1 |    +6.7 |    +4.1 |    +0.0 |   +11.3 | 2.046 |
| `1^3F_4` | 2.038 |    +0.4 |   +47.6 |   -77.4 |   -29.7 |    -1.3 |    +0.0 |   -30.6 | 2.008 |
| `1^1G_4` | 2.334 |    -0.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -0.4 | 2.333 |
| `1^3G_4` | 2.334 |    +0.1 |    -9.5 |   +17.6 |    +8.2 |    +2.2 |    +0.0 |   +10.6 | 2.344 |
| `1^3G_5` | 2.336 |    +0.1 |   +35.4 |   -74.7 |   -39.3 |    -0.8 |    +0.0 |   -40.0 | 2.296 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 1.456 | 1.455 |    +9.6 |  -0.999 |  +0.045 |
| `1^3D_1` | 1.664 | 1.666 |    +9.6 |  +0.045 |  +0.998 |
| `3^3S_1` | 1.998 | 1.998 |    +2.1 |  +1.000 |  -0.013 |
| `2^3D_1` | 2.153 | 2.153 |    +2.1 |  -0.013 |  -1.000 |
| `2^3P_2` | 1.823 | 1.823 |    +1.4 |  -1.000 |  +0.006 |
| `1^3F_2` | 2.054 | 2.055 |    +1.4 |  +0.006 |  +1.000 |
| `2^3D_3` | 2.131 | 2.131 |    +0.3 |  -1.000 |  +0.001 |
| `1^3G_3` | 2.370 | 2.370 |    +0.3 |  +0.001 |  +1.000 |

Mean absolute residual: 3.3 MeV.
Max absolute residual: 8.4 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.682 | 1.683 |    +1.0 |
| `2D` | 4 | 2.138 | 2.141 |    +3.2 |
| `1F` | 4 | 2.032 | 2.032 |    -0.3 |
| `1G` | 4 | 2.331 | 2.332 |    +0.8 |
| `1P` | 4 | 1.252 | 1.249 |    -2.7 |
| `2P` | 4 | 1.807 | 1.806 |    -0.2 |
| `1S` | 2 | 0.615 | 0.616 |    +0.5 |
| `2S` | 2 | 1.412 | 1.414 |    +1.7 |
| `3S` | 2 | 1.970 | 1.966 |    -3.6 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
