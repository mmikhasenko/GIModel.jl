# Residuals: b_flavored (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.305 |    -4.5 | high |
| `2^1S_0` | 5.900 | 5.900 |    +0.3 | high |
| `1^3S_1` | 5.370 | 5.370 |    +0.1 | high |
| `2^3S_1` | 5.930 | 5.935 |    +4.8 | high |
| `1^3P_2` | 5.800 | 5.792 |    -8.4 | high |
| `1^3D_3` | 6.110 | 6.092 |   -17.8 | medium |
| `1^3F_4` | 6.360 | 6.343 |   -16.7 | medium |
| `1^1S_0` | 5.390 | 5.387 |    -3.3 | high |
| `2^1S_0` | 5.980 | 5.980 |    +0.4 | high |
| `1^3S_1` | 5.450 | 5.448 |    -2.0 | high |
| `2^3S_1` | 6.010 | 6.015 |    +4.6 | high |
| `1^3P_2` | 5.880 | 5.870 |   -10.2 | high |
| `1^3D_3` | 6.180 | 6.174 |    -6.2 | medium |
| `1^3F_4` | 6.430 | 6.426 |    -4.3 | medium |
| `1^1S_0` | 6.270 | 6.261 |    -9.0 | high |
| `2^1S_0` | 6.850 | 6.851 |    +1.0 | high |
| `1^3S_1` | 6.340 | 6.333 |    -7.2 | high |
| `2^3S_1` | 6.890 | 6.889 |    -0.6 | high |
| `1^3P_2` | 6.770 | 6.757 |   -12.7 | high |
| `1^3D_3` | 7.040 | 7.041 |    +1.1 | medium |
| `1^3F_4` | 7.270 | 7.270 |    -0.3 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 5.355 |   -49.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -49.8 | 5.305 |
| `2^1S_0` | 5.926 |   -25.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -25.9 | 5.900 |
| `1^3S_1` | 5.355 |   +14.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +14.8 | 5.370 |
| `2^3S_1` | 5.926 |    +8.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.6 | 5.935 |
| `1^3P_2` | 5.786 |    +0.0 |  +120.2 |  -112.5 |    +7.7 |    -2.2 |    +0.0 |    +5.4 | 5.792 |
| `1^3D_3` | 6.109 |    +0.0 |   +51.6 |   -68.0 |   -16.4 |    -0.8 |    +0.0 |   -17.1 | 6.092 |
| `1^3F_4` | 6.375 |    +0.0 |   +42.6 |   -74.1 |   -31.5 |    -0.5 |    +0.0 |   -32.0 | 6.343 |
| `1^1S_0` | 5.434 |   -47.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -47.4 | 5.387 |
| `2^1S_0` | 6.006 |   -25.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -25.7 | 5.980 |
| `1^3S_1` | 5.434 |   +13.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.9 | 5.448 |
| `2^3S_1` | 6.006 |    +8.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.4 | 6.015 |
| `1^3P_2` | 5.864 |    +0.0 |   +48.4 |   -41.1 |    +7.2 |    -1.7 |    +0.0 |    +5.5 | 5.870 |
| `1^3D_3` | 6.181 |    +0.0 |   +30.5 |   -37.3 |    -6.7 |    -0.8 |    +0.0 |    -7.6 | 6.174 |
| `1^3F_4` | 6.442 |    +0.0 |   +24.4 |   -39.9 |   -15.6 |    -0.5 |    +0.0 |   -16.1 | 6.426 |
| `1^1S_0` | 6.318 |   -56.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -56.7 | 6.261 |
| `2^1S_0` | 6.881 |   -29.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -29.7 | 6.851 |
| `1^3S_1` | 6.318 |   +15.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +15.1 | 6.333 |
| `2^3S_1` | 6.881 |    +8.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +8.7 | 6.889 |
| `1^3P_2` | 6.751 |    +0.0 |   +17.0 |    -8.8 |    +8.2 |    -1.8 |    +0.0 |    +6.4 | 6.757 |
| `1^3D_3` | 7.039 |    +0.0 |   +11.7 |    -8.8 |    +2.9 |    -0.9 |    +0.0 |    +2.0 | 7.041 |
| `1^3F_4` | 7.270 |    +0.0 |    +8.0 |    -8.1 |    -0.1 |    -0.5 |    +0.0 |    -0.6 | 7.270 |

## Fine-Structure Mass Convention (audit note)

Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 4.977000 | 0.220000 | `disabled` |
| `2^1S_0` | 4.977000 | 0.220000 | `disabled` |
| `1^3S_1` | 4.977000 | 0.220000 | `disabled` |
| `2^3S_1` | 4.977000 | 0.220000 | `disabled` |
| `1^3P_2` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 0.419000 | `disabled` |
| `2^1S_0` | 4.977000 | 0.419000 | `disabled` |
| `1^3S_1` | 4.977000 | 0.419000 | `disabled` |
| `2^3S_1` | 4.977000 | 0.419000 | `disabled` |
| `1^3P_2` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 1.628000 | `disabled` |
| `2^1S_0` | 4.977000 | 1.628000 | `disabled` |
| `1^3S_1` | 4.977000 | 1.628000 | `disabled` |
| `2^3S_1` | 4.977000 | 1.628000 | `disabled` |
| `1^3P_2` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 5.5 MeV.
Max absolute residual: 17.8 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.436 |    -7.6 |
| `1F` | 3 | 6.687 | 6.680 |    -7.1 |
| `1P` | 3 | 6.150 | 6.140 |   -10.4 |
| `1S` | 6 | 5.704 | 5.700 |    -3.7 |
| `2S` | 6 | 6.268 | 6.271 |    +2.3 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
