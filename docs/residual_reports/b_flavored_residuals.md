# Residuals: b_flavored (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.281 |   -28.6 | high |
| `2^1S_0` | 5.900 | 5.873 |   -26.6 | high |
| `1^3S_1` | 5.370 | 5.380 |    +9.9 | high |
| `2^3S_1` | 5.930 | 5.944 |   +13.8 | high |
| `1^3P_2` | 5.800 | 5.782 |   -17.9 | high |
| `1^3D_3` | 6.110 | 6.054 |   -56.4 | medium |
| `1^3F_4` | 6.360 | 6.282 |   -78.0 | medium |
| `1^1S_0` | 5.390 | 5.374 |   -16.3 | high |
| `2^1S_0` | 5.980 | 5.963 |   -16.5 | high |
| `1^3S_1` | 5.450 | 5.454 |    +4.2 | high |
| `2^3S_1` | 6.010 | 6.020 |   +10.4 | high |
| `1^3P_2` | 5.880 | 5.869 |   -10.6 | high |
| `1^3D_3` | 6.180 | 6.165 |   -15.4 | medium |
| `1^3F_4` | 6.430 | 6.410 |   -20.4 | medium |
| `1^1S_0` | 6.270 | 6.261 |    -9.0 | high |
| `2^1S_0` | 6.850 | 6.846 |    -4.2 | high |
| `1^3S_1` | 6.340 | 6.337 |    -3.4 | high |
| `2^3S_1` | 6.890 | 6.892 |    +2.3 | high |
| `1^3P_2` | 6.770 | 6.758 |   -11.6 | high |
| `1^3D_3` | 7.040 | 7.041 |    +1.3 | medium |
| `1^3F_4` | 7.270 | 7.269 |    -0.6 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 5.355 |   -73.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -73.9 | 5.281 |
| `2^1S_0` | 5.926 |   -52.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -52.8 | 5.873 |
| `1^3S_1` | 5.355 |   +24.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +24.6 | 5.380 |
| `2^3S_1` | 5.926 |   +17.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +17.6 | 5.944 |
| `1^3P_2` | 5.786 |    +0.0 |  +151.6 |  -152.7 |    -1.1 |    -2.9 |    -4.0 | 5.782 |
| `1^3D_3` | 6.109 |    +0.0 |  +126.4 |  -180.2 |   -53.8 |    -2.0 |   -55.8 | 6.054 |
| `1^3F_4` | 6.375 |    +0.0 |   +99.0 |  -191.0 |   -92.0 |    -1.3 |   -93.3 | 6.282 |
| `1^1S_0` | 5.434 |   -60.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -60.4 | 5.374 |
| `2^1S_0` | 6.006 |   -42.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -42.7 | 5.963 |
| `1^3S_1` | 5.434 |   +20.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +20.1 | 5.454 |
| `2^3S_1` | 6.006 |   +14.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +14.2 | 6.020 |
| `1^3P_2` | 5.864 |    +0.0 |   +69.2 |   -61.5 |    +7.7 |    -2.5 |    +5.2 | 5.869 |
| `1^3D_3` | 6.181 |    +0.0 |   +54.4 |   -69.7 |   -15.3 |    -1.6 |   -16.9 | 6.165 |
| `1^3F_4` | 6.442 |    +0.0 |   +41.7 |   -73.0 |   -31.3 |    -1.0 |   -32.3 | 6.410 |
| `1^1S_0` | 6.318 |   -56.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -56.7 | 6.261 |
| `2^1S_0` | 6.881 |   -34.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -34.9 | 6.846 |
| `1^3S_1` | 6.318 |   +18.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +18.9 | 6.337 |
| `2^3S_1` | 6.881 |   +11.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +11.6 | 6.892 |
| `1^3P_2` | 6.751 |    +0.0 |   +20.1 |   -10.5 |    +9.6 |    -2.1 |    +7.5 | 6.758 |
| `1^3D_3` | 7.039 |    +0.0 |   +13.4 |   -10.1 |    +3.2 |    -1.0 |    +2.2 | 7.041 |
| `1^3F_4` | 7.270 |    +0.0 |    +9.8 |   -10.1 |    -0.2 |    -0.6 |    -0.9 | 7.269 |

## Fine-Structure Mass Convention (audit note)

Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 17.0 MeV.
Max absolute residual: 78.0 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.420 |   -23.5 |
| `1F` | 3 | 6.687 | 6.654 |   -33.0 |
| `1P` | 3 | 6.150 | 6.137 |   -13.3 |
| `1S` | 6 | 5.704 | 5.702 |    -1.8 |
| `2S` | 6 | 6.268 | 6.271 |    +2.7 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
