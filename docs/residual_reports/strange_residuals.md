# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.427 |   -43.0 | high |
| `2^1S_0` | 1.450 | 1.441 |    -8.6 | high |
| `3^1S_0` | 2.020 | 2.003 |   -17.2 | high |
| `1^3S_1` | 0.900 | 0.898 |    -2.3 | high |
| `2^3S_1` | 1.580 | 1.583 |    +3.0 | high |
| `1^3D_1` | 1.780 | 1.795 |   +15.3 | high |
| `3^3S_1` | 2.110 | 2.110 |    +0.3 | high |
| `2^3D_1` | 2.250 | 2.252 |    +2.4 | high |
| `1^3P_0` | 1.240 | 1.276 |   +35.9 | high |
| `2^3P_0` | 1.890 | 1.886 |    -3.8 | high |
| `1^1P_1` | 1.340 | 1.384 |   +43.7 | high |
| `1^3P_1` | 1.380 | 1.402 |   +22.5 | high |
| `2^1P_1` | 1.900 | 1.922 |   +21.7 | high |
| `2^3P_1` | 1.930 | 1.935 |    +5.3 | high |
| `1^3P_2` | 1.430 | 1.394 |   -36.0 | high |
| `2^3P_2` | 1.940 | 1.921 |   -19.5 | high |
| `1^3F_2` | 2.150 | 2.156 |    +6.3 | high |
| `1^1D_2` | 1.780 | 1.796 |   +16.1 | high |
| `1^3D_2` | 1.810 | 1.814 |    +3.7 | high |
| `2^1D_2` | 2.230 | 2.244 |   +14.3 | high |
| `2^3D_2` | 2.260 | 2.257 |    -2.6 | high |
| `1^3D_3` | 1.790 | 1.784 |    -6.1 | high |
| `2^3D_3` | 2.240 | 2.232 |    -8.5 | high |
| `1^3G_3` | 2.460 | 2.457 |    -3.5 | medium |
| `1^1F_3` | 2.120 | 2.132 |   +12.0 | high |
| `1^3F_3` | 2.150 | 2.146 |    -4.3 | high |
| `1^3F_4` | 2.110 | 2.108 |    -2.2 | high |
| `1^1G_4` | 2.410 | 2.423 |   +12.7 | medium |
| `1^3G_4` | 2.440 | 2.434 |    -6.4 | medium |
| `1^3G_5` | 2.390 | 2.392 |    +2.3 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.814 |  -387.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -387.1 | 0.427 |
| `2^1S_0` | 1.545 |  -103.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -103.8 | 1.441 |
| `3^1S_0` | 2.080 |   -77.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -77.4 | 2.003 |
| `1^3S_1` | 0.814 |   +83.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +83.6 | 0.898 |
| `2^3S_1` | 1.545 |   +37.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +37.8 | 1.583 |
| `1^3D_1` | 1.796 |    +0.0 |   -70.4 |   +83.0 |   +12.6 |   -13.4 |    -0.8 | 1.795 |
| `3^3S_1` | 2.080 |   +30.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +30.0 | 2.110 |
| `2^3D_1` | 2.244 |    +0.0 |   -45.6 |   +61.4 |   +15.8 |    -7.8 |    +8.0 | 2.252 |
| `1^3P_0` | 1.384 |    +0.0 |  -137.9 |  +102.7 |   -35.2 |   -72.7 |  -107.8 | 1.276 |
| `2^3P_0` | 1.922 |    +0.0 |   -67.7 |   +63.7 |    -4.1 |   -31.4 |   -35.5 | 1.886 |
| `1^1P_1` | 1.384 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.384 |
| `1^3P_1` | 1.384 |    +0.0 |   -68.9 |   +51.4 |   -17.6 |   +36.3 |   +18.7 | 1.402 |
| `2^1P_1` | 1.922 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.922 |
| `2^3P_1` | 1.922 |    +0.0 |   -33.9 |   +31.8 |    -2.0 |   +15.7 |   +13.7 | 1.935 |
| `1^3P_2` | 1.384 |    +0.0 |   +68.9 |   -51.4 |   +17.6 |    -7.3 |   +10.3 | 1.394 |
| `2^3P_2` | 1.922 |    +0.0 |   +33.9 |   -31.8 |    +2.0 |    -3.1 |    -1.1 | 1.921 |
| `1^3F_2` | 2.132 |    +0.0 |   -43.5 |   +72.9 |   +29.5 |    -5.1 |   +24.3 | 2.156 |
| `1^1D_2` | 1.796 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.796 |
| `1^3D_2` | 1.796 |    +0.0 |   -23.5 |   +27.7 |    +4.2 |   +13.4 |   +17.6 | 1.814 |
| `2^1D_2` | 2.244 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.244 |
| `2^3D_2` | 2.244 |    +0.0 |   -15.2 |   +20.5 |    +5.3 |    +7.8 |   +13.1 | 2.257 |
| `1^3D_3` | 1.796 |    +0.0 |   +47.0 |   -55.3 |    -8.4 |    -3.8 |   -12.2 | 1.784 |
| `2^3D_3` | 2.244 |    +0.0 |   +30.4 |   -40.9 |   -10.6 |    -2.2 |   -12.8 | 2.232 |
| `1^3G_3` | 2.423 |    +0.0 |   -30.0 |   +66.4 |   +36.4 |    -2.6 |   +33.8 | 2.457 |
| `1^1F_3` | 2.132 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.132 |
| `1^3F_3` | 2.132 |    +0.0 |   -10.9 |   +18.2 |    +7.4 |    +6.4 |   +13.8 | 2.146 |
| `1^3F_4` | 2.132 |    +0.0 |   +32.6 |   -54.7 |   -22.1 |    -2.1 |   -24.2 | 2.108 |
| `1^1G_4` | 2.423 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.423 |
| `1^3G_4` | 2.423 |    +0.0 |    -6.0 |   +13.3 |    +7.3 |    +3.6 |   +10.9 | 2.434 |
| `1^3G_5` | 2.423 |    +0.0 |   +24.0 |   -53.1 |   -29.1 |    -1.3 |   -30.4 | 2.392 |

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

Mean absolute residual: 12.7 MeV.
Max absolute residual: 43.7 MeV.

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
| `1S` | 2 | 0.792 | 0.780 |   -12.5 |
| `2S` | 2 | 1.548 | 1.548 |    +0.1 |
| `3S` | 2 | 2.087 | 2.083 |    -4.1 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
