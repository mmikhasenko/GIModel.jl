# Residuals: b_flavored (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.305 |    -4.5 | high |
| `2^1S_0` | 5.900 | 5.900 |    +0.3 | high |
| `1^3S_1` | 5.370 | 5.370 |    +0.1 | high |
| `2^3S_1` | 5.930 | 5.935 |    +4.8 | high |
| `1^3P_2` | 5.800 | 5.784 |   -16.4 | high |
| `1^3D_3` | 6.110 | 6.082 |   -28.5 | medium |
| `1^3F_4` | 6.360 | 6.331 |   -29.1 | medium |
| `1^1S_0` | 5.390 | 5.387 |    -3.3 | high |
| `2^1S_0` | 5.980 | 5.980 |    +0.4 | high |
| `1^3S_1` | 5.450 | 5.448 |    -2.0 | high |
| `2^3S_1` | 6.010 | 6.015 |    +4.6 | high |
| `1^3P_2` | 5.880 | 5.867 |   -13.0 | high |
| `1^3D_3` | 6.180 | 6.170 |    -9.8 | medium |
| `1^3F_4` | 6.430 | 6.422 |    -8.4 | medium |
| `1^1S_0` | 6.270 | 6.261 |    -9.0 | high |
| `2^1S_0` | 6.850 | 6.851 |    +1.0 | high |
| `1^3S_1` | 6.340 | 6.333 |    -7.2 | high |
| `2^3S_1` | 6.890 | 6.889 |    -0.6 | high |
| `1^3P_2` | 6.770 | 6.757 |   -12.8 | high |
| `1^3D_3` | 7.040 | 7.041 |    +0.9 | medium |
| `1^3F_4` | 7.270 | 7.269 |    -0.5 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 5.355 |   -49.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -49.8 | 5.305 |
| `2^1S_0` | 5.926 |   -25.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -25.9 | 5.900 |
| `1^3S_1` | 5.355 |   +14.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +14.8 | 5.370 |
| `2^3S_1` | 5.926 |    +8.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.6 | 5.935 |
| `1^3P_2` | 5.786 |    +0.0 |   +77.2 |   -78.5 |    -1.2 |    -1.4 |    -2.6 | 5.784 |
| `1^3D_3` | 6.109 |    +0.0 |   +61.9 |   -88.8 |   -26.9 |    -0.9 |   -27.8 | 6.082 |
| `1^3F_4` | 6.375 |    +0.0 |   +47.1 |   -90.9 |   -43.8 |    -0.6 |   -44.4 | 6.331 |
| `1^1S_0` | 5.434 |   -47.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -47.4 | 5.387 |
| `2^1S_0` | 6.006 |   -25.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -25.7 | 5.980 |
| `1^3S_1` | 5.434 |   +13.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.9 | 5.448 |
| `2^3S_1` | 6.006 |    +8.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.4 | 6.015 |
| `1^3P_2` | 5.864 |    +0.0 |   +44.3 |   -39.9 |    +4.4 |    -1.5 |    +2.8 | 5.867 |
| `1^3D_3` | 6.181 |    +0.0 |   +34.1 |   -44.3 |   -10.3 |    -0.9 |   -11.2 | 6.170 |
| `1^3F_4` | 6.442 |    +0.0 |   +25.6 |   -45.3 |   -19.7 |    -0.6 |   -20.3 | 6.422 |
| `1^1S_0` | 6.318 |   -56.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -56.7 | 6.261 |
| `2^1S_0` | 6.881 |   -29.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -29.7 | 6.851 |
| `1^3S_1` | 6.318 |   +15.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +15.1 | 6.333 |
| `2^3S_1` | 6.881 |    +8.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.7 | 6.889 |
| `1^3P_2` | 6.751 |    +0.0 |   +17.2 |    -9.1 |    +8.1 |    -1.8 |    +6.3 | 6.757 |
| `1^3D_3` | 7.039 |    +0.0 |   +11.6 |    -8.8 |    +2.7 |    -0.9 |    +1.8 | 7.041 |
| `1^3F_4` | 7.270 |    +0.0 |    +8.5 |    -8.8 |    -0.3 |    -0.5 |    -0.8 | 7.269 |

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

Mean absolute residual: 7.5 MeV.
Max absolute residual: 29.1 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.431 |   -12.5 |
| `1F` | 3 | 6.687 | 6.674 |   -12.7 |
| `1P` | 3 | 6.150 | 6.136 |   -14.1 |
| `1S` | 6 | 5.704 | 5.700 |    -3.7 |
| `2S` | 6 | 6.268 | 6.271 |    +2.3 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
