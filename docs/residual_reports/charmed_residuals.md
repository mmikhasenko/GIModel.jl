# Residuals: charmed (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.860 |   -20.2 | high |
| `2^1S_0` | 2.580 | 2.575 |    -5.2 | high |
| `1^3S_1` | 2.040 | 2.039 |    -0.9 | high |
| `2^3S_1` | 2.640 | 2.648 |    +7.6 | high |
| `1^3D_1` | 2.820 | 2.856 |   +36.2 | high |
| `1^3P_0` | 2.400 | 2.430 |   +30.4 | high |
| `1^1P_1` | 2.440 | 2.475 |   +35.3 | high |
| `1^3P_1` | 2.490 | 2.489 |    -1.4 | high |
| `1^3P_2` | 2.500 | 2.476 |   -23.7 | high |
| `1^3D_3` | 2.830 | 2.807 |   -23.3 | high |
| `1^3F_4` | 3.110 | 3.084 |   -25.9 | medium |
| `1^1S_0` | 1.980 | 1.962 |   -18.0 | high |
| `2^1S_0` | 2.670 | 2.667 |    -3.1 | high |
| `1^3S_1` | 2.130 | 2.125 |    -4.7 | high |
| `2^3S_1` | 2.730 | 2.737 |    +6.8 | high |
| `1^3D_1` | 2.900 | 2.913 |   +13.4 | high |
| `1^3P_0` | 2.480 | 2.505 |   +24.7 | high |
| `1^1P_1` | 2.530 | 2.564 |   +34.3 | high |
| `1^3P_1` | 2.570 | 2.572 |    +2.1 | high |
| `1^3P_2` | 2.590 | 2.572 |   -18.4 | high |
| `1^3D_3` | 2.920 | 2.906 |   -14.4 | high |
| `1^3F_4` | 3.190 | 3.183 |    -7.0 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.001 |  -141.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -141.5 | 1.860 |
| `2^1S_0` | 2.629 |   -54.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -54.0 | 2.575 |
| `1^3S_1` | 2.001 |   +37.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +37.8 | 2.039 |
| `2^3S_1` | 2.629 |   +18.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +18.8 | 2.648 |
| `1^3D_1` | 2.831 |    +0.0 |   -84.6 |  +117.5 |   +32.9 |    -7.5 |   +25.4 | 2.856 |
| `1^3P_0` | 2.475 |    +0.0 |  -148.8 |  +139.7 |    -9.1 |   -35.7 |   -44.8 | 2.430 |
| `1^1P_1` | 2.475 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.475 |
| `1^3P_1` | 2.475 |    +0.0 |   -74.4 |   +69.8 |    -4.6 |   +17.9 |   +13.3 | 2.489 |
| `1^3P_2` | 2.475 |    +0.0 |   +74.4 |   -69.8 |    +4.6 |    -3.6 |    +1.0 | 2.476 |
| `1^3D_3` | 2.831 |    +0.0 |   +56.4 |   -78.3 |   -22.0 |    -2.2 |   -24.1 | 2.807 |
| `1^3F_4` | 3.124 |    +0.0 |   +41.3 |   -79.8 |   -38.5 |    -1.3 |   -39.8 | 3.084 |
| `1^1S_0` | 2.091 |  -129.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -129.3 | 1.962 |
| `2^1S_0` | 2.719 |   -52.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -52.4 | 2.667 |
| `1^3S_1` | 2.091 |   +34.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +34.1 | 2.125 |
| `2^3S_1` | 2.719 |   +17.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +17.6 | 2.737 |
| `1^3D_1` | 2.913 |    +0.0 |   -51.5 |   +59.5 |    +7.9 |    -7.6 |    +0.3 | 2.913 |
| `1^3P_0` | 2.564 |    +0.0 |   -93.6 |   +71.6 |   -22.0 |   -37.6 |   -59.6 | 2.505 |
| `1^1P_1` | 2.564 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.564 |
| `1^3P_1` | 2.564 |    +0.0 |   -46.8 |   +35.8 |   -11.0 |   +18.8 |    +7.8 | 2.572 |
| `1^3P_2` | 2.564 |    +0.0 |   +46.8 |   -35.8 |   +11.0 |    -3.8 |    +7.2 | 2.572 |
| `1^3D_3` | 2.913 |    +0.0 |   +34.4 |   -39.6 |    -5.3 |    -2.2 |    -7.5 | 2.906 |
| `1^3F_4` | 3.200 |    +0.0 |   +25.0 |   -40.5 |   -15.5 |    -1.3 |   -16.8 | 3.183 |

## Fine-Structure Mass Convention (audit note)

Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 16.2 MeV.
Max absolute residual: 36.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.865 |    -5.8 |
| `1F` | 2 | 3.150 | 3.134 |   -16.5 |
| `1P` | 8 | 2.518 | 2.520 |    +2.3 |
| `1S` | 4 | 2.046 | 2.039 |    -6.9 |
| `2S` | 4 | 2.670 | 2.674 |    +4.4 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
