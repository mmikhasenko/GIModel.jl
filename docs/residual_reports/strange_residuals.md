# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.240 |  -229.8 | high |
| `2^1S_0` | 1.450 | 1.191 |  -258.9 | high |
| `3^1S_0` | 2.020 | 1.790 |  -230.3 | high |
| `1^3S_1` | 0.900 | 1.005 |  +105.4 | high |
| `2^3S_1` | 1.580 | 1.663 |   +83.2 | high |
| `1^3D_1` | 1.780 | 1.780 |    -0.3 | high |
| `3^3S_1` | 2.110 | 2.177 |   +67.1 | high |
| `2^3D_1` | 2.250 | 2.219 |   -31.4 | high |
| `1^3P_0` | 1.240 | 1.088 |  -151.8 | high |
| `2^3P_0` | 1.890 | 1.720 |  -169.9 | high |
| `1^1P_1` | 1.340 | 1.384 |   +43.7 | high |
| `1^3P_1` | 1.380 | 1.430 |   +50.1 | high |
| `2^1P_1` | 1.900 | 1.922 |   +21.7 | high |
| `2^3P_1` | 1.930 | 1.945 |   +15.3 | high |
| `1^3P_2` | 1.430 | 1.415 |   -15.0 | high |
| `2^3P_2` | 1.940 | 1.948 |    +7.8 | high |
| `1^3F_2` | 2.150 | 2.189 |   +38.8 | high |
| `1^1D_2` | 1.780 | 1.796 |   +16.1 | high |
| `1^3D_2` | 1.810 | 1.841 |   +31.4 | high |
| `2^1D_2` | 2.230 | 2.244 |   +14.3 | high |
| `2^3D_2` | 2.260 | 2.276 |   +16.5 | high |
| `1^3D_3` | 1.790 | 1.771 |   -19.2 | high |
| `2^3D_3` | 2.240 | 2.232 |    -7.6 | high |
| `1^3G_3` | 2.460 | 2.513 |   +52.9 | medium |
| `1^1F_3` | 2.120 | 2.132 |   +12.0 | high |
| `1^3F_3` | 2.150 | 2.169 |   +19.3 | high |
| `1^3F_4` | 2.110 | 2.071 |   -38.6 | high |
| `1^1G_4` | 2.410 | 2.423 |   +12.7 | medium |
| `1^3G_4` | 2.440 | 2.454 |   +13.7 | medium |
| `1^3G_5` | 2.390 | 2.340 |   -50.0 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.814 |  -573.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -573.9 | 0.240 |
| `2^1S_0` | 1.545 |  -354.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -354.1 | 1.191 |
| `3^1S_0` | 2.080 |  -290.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -290.6 | 1.790 |
| `1^3S_1` | 0.814 |  +191.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +191.3 | 1.005 |
| `2^3S_1` | 1.545 |  +118.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +118.0 | 1.663 |
| `1^3D_1` | 1.796 |    +0.0 |  -182.1 |  +203.8 |   +21.7 |   -38.0 |   -16.4 | 1.780 |
| `3^3S_1` | 2.080 |   +96.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +96.9 | 2.177 |
| `2^3D_1` | 2.244 |    +0.0 |  -159.8 |  +164.6 |    +4.8 |   -30.5 |   -25.7 | 2.219 |
| `1^3P_0` | 1.384 |    +0.0 |  -341.8 |  +240.4 |  -101.4 |  -194.1 |  -295.5 | 1.088 |
| `2^3P_0` | 1.922 |    +0.0 |  -247.2 |  +170.1 |   -77.1 |  -124.5 |  -201.6 | 1.720 |
| `1^1P_1` | 1.384 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.384 |
| `1^3P_1` | 1.384 |    +0.0 |  -170.9 |  +120.2 |   -50.7 |   +97.1 |   +46.4 | 1.430 |
| `2^1P_1` | 1.922 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.922 |
| `2^3P_1` | 1.922 |    +0.0 |  -123.6 |   +85.1 |   -38.6 |   +62.2 |   +23.7 | 1.945 |
| `1^3P_2` | 1.384 |    +0.0 |  +170.9 |  -120.2 |   +50.7 |   -19.4 |   +31.3 | 1.415 |
| `2^3P_2` | 1.922 |    +0.0 |  +123.6 |   -85.1 |   +38.6 |   -12.4 |   +26.1 | 1.948 |
| `1^3F_2` | 2.132 |    +0.0 |  -117.2 |  +189.5 |   +72.2 |   -15.4 |   +56.9 | 2.189 |
| `1^1D_2` | 1.796 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.796 |
| `1^3D_2` | 1.796 |    +0.0 |   -60.7 |   +67.9 |    +7.2 |   +38.0 |   +45.3 | 1.841 |
| `2^1D_2` | 2.244 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.244 |
| `2^3D_2` | 2.244 |    +0.0 |   -53.3 |   +54.9 |    +1.6 |   +30.5 |   +32.1 | 2.276 |
| `1^3D_3` | 1.796 |    +0.0 |  +121.4 |  -135.9 |   -14.5 |   -10.9 |   -25.3 | 1.771 |
| `2^3D_3` | 2.244 |    +0.0 |  +106.5 |  -109.8 |    -3.2 |    -8.7 |   -11.9 | 2.232 |
| `1^3G_3` | 2.423 |    +0.0 |   -84.5 |  +182.8 |   +98.3 |    -8.1 |   +90.1 | 2.513 |
| `1^1F_3` | 2.132 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.132 |
| `1^3F_3` | 2.132 |    +0.0 |   -29.3 |   +47.4 |   +18.1 |   +19.2 |   +37.3 | 2.169 |
| `1^3F_4` | 2.132 |    +0.0 |   +87.9 |  -142.1 |   -54.2 |    -6.4 |   -60.6 | 2.071 |
| `1^1G_4` | 2.423 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.423 |
| `1^3G_4` | 2.423 |    +0.0 |   -16.9 |   +36.6 |   +19.7 |   +11.4 |   +31.0 | 2.454 |
| `1^3G_5` | 2.423 |    +0.0 |   +67.6 |  -146.2 |   -78.6 |    -4.1 |   -82.7 | 2.340 |

## Fine-Structure Mass Convention (audit note)

Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `3^1S_0` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `3^3S_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3D_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3P_0` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1P_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3P_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3P_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1D_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1D_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3D_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3D_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3G_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1F_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1G_4` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3G_4` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3G_5` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 60.8 MeV.
Max absolute residual: 258.9 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.791 | 1.796 |    +5.1 |
| `2D` | 4 | 2.244 | 2.244 |    +0.3 |
| `1F` | 4 | 2.130 | 2.132 |    +2.3 |
| `1G` | 4 | 2.421 | 2.423 |    +1.6 |
| `1P` | 4 | 1.379 | 1.384 |    +4.6 |
| `2P` | 4 | 1.923 | 1.922 |    -1.7 |
| `1S` | 2 | 0.792 | 0.814 |   +21.6 |
| `2S` | 2 | 1.548 | 1.545 |    -2.3 |
| `3S` | 2 | 2.087 | 2.080 |    -7.2 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
