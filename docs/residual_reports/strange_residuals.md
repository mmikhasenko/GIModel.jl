# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.427 |   -43.0 | high |
| `2^1S_0` | 1.450 | 1.441 |    -8.6 | high |
| `3^1S_0` | 2.020 | 2.003 |   -17.2 | high |
| `1^3S_1` | 0.900 | 0.898 |    -2.3 | high |
| `2^3S_1` | 1.580 | 1.583 |    +3.0 | high |
| `1^3D_1` | 1.780 | 1.781 |    +1.0 | high |
| `3^3S_1` | 2.110 | 2.110 |    +0.3 | high |
| `2^3D_1` | 2.250 | 2.242 |    -8.2 | high |
| `1^3P_0` | 1.240 | 1.259 |   +19.4 | high |
| `2^3P_0` | 1.890 | 1.876 |   -14.2 | high |
| `1^1P_1` | 1.340 | 1.382 |   +42.2 | high |
| `1^3P_1` | 1.380 | 1.396 |   +15.7 | high |
| `2^1P_1` | 1.900 | 1.918 |   +18.2 | high |
| `2^3P_1` | 1.930 | 1.934 |    +3.6 | high |
| `1^3P_2` | 1.430 | 1.402 |   -27.7 | high |
| `2^3P_2` | 1.940 | 1.926 |   -14.2 | high |
| `1^3F_2` | 2.150 | 2.143 |    -7.1 | high |
| `1^1D_2` | 1.780 | 1.787 |    +7.1 | high |
| `1^3D_2` | 1.810 | 1.818 |    +7.9 | high |
| `2^1D_2` | 2.230 | 2.236 |    +6.2 | high |
| `2^3D_2` | 2.260 | 2.262 |    +2.0 | high |
| `1^3D_3` | 1.790 | 1.793 |    +3.4 | high |
| `2^3D_3` | 2.240 | 2.239 |    -1.4 | high |
| `1^3G_3` | 2.460 | 2.444 |   -16.4 | medium |
| `1^1F_3` | 2.120 | 2.119 |    -1.1 | high |
| `1^3F_3` | 2.150 | 2.156 |    +5.5 | high |
| `1^3F_4` | 2.110 | 2.118 |    +7.8 | high |
| `1^1G_4` | 2.410 | 2.408 |    -2.4 | medium |
| `1^3G_4` | 2.440 | 2.446 |    +6.1 | medium |
| `1^3G_5` | 2.390 | 2.403 |   +12.7 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.814 |  -387.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -387.1 | 0.427 |
| `2^1S_0` | 1.545 |  -103.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -103.8 | 1.441 |
| `3^1S_0` | 2.080 |   -77.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -77.4 | 2.003 |
| `1^3S_1` | 0.814 |   +83.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +83.6 | 0.898 |
| `2^3S_1` | 1.545 |   +37.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +37.8 | 1.583 |
| `1^3D_1` | 1.796 |    +0.0 |   -70.4 |   +68.8 |    -1.7 |   -13.4 |    +0.0 |   -15.0 | 1.781 |
| `3^3S_1` | 2.080 |   +30.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +30.0 | 2.110 |
| `2^3D_1` | 2.244 |    +0.0 |   -45.6 |   +50.8 |    +5.3 |    -7.8 |    +0.0 |    -2.5 | 2.242 |
| `1^3P_0` | 1.384 |    +0.0 |  -137.9 |   +86.2 |   -51.7 |   -72.7 |    +0.0 |  -124.3 | 1.259 |
| `2^3P_0` | 1.922 |    +0.0 |   -67.7 |   +53.2 |   -14.5 |   -31.4 |    +0.0 |   -45.9 | 1.876 |
| `1^1P_1` | 1.384 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.382 |
| `1^3P_1` | 1.384 |    +0.0 |   -68.9 |   +43.1 |   -25.8 |   +36.3 |    +0.0 |   +10.5 | 1.396 |
| `2^1P_1` | 1.922 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.918 |
| `2^3P_1` | 1.922 |    +0.0 |   -33.9 |   +26.6 |    -7.2 |   +15.7 |    +0.0 |    +8.5 | 1.934 |
| `1^3P_2` | 1.384 |    +0.0 |   +68.9 |   -43.1 |   +25.8 |    -7.3 |    +0.0 |   +18.6 | 1.402 |
| `2^3P_2` | 1.922 |    +0.0 |   +33.9 |   -26.6 |    +7.2 |    -3.1 |    +0.0 |    +4.1 | 1.926 |
| `1^3F_2` | 2.132 |    +0.0 |   -43.5 |   +59.6 |   +16.1 |    -5.1 |    +0.0 |   +11.0 | 2.143 |
| `1^1D_2` | 1.796 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.787 |
| `1^3D_2` | 1.796 |    +0.0 |   -23.5 |   +22.9 |    -0.6 |   +13.4 |    +0.0 |   +12.8 | 1.818 |
| `2^1D_2` | 2.244 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.236 |
| `2^3D_2` | 2.244 |    +0.0 |   -15.2 |   +16.9 |    +1.8 |    +7.8 |    +0.0 |    +9.6 | 2.262 |
| `1^3D_3` | 1.796 |    +0.0 |   +47.0 |   -45.8 |    +1.1 |    -3.8 |    +0.0 |    -2.7 | 1.793 |
| `2^3D_3` | 2.244 |    +0.0 |   +30.4 |   -33.9 |    -3.5 |    -2.2 |    +0.0 |    -5.7 | 2.239 |
| `1^3G_3` | 2.423 |    +0.0 |   -30.0 |   +53.5 |   +23.5 |    -2.6 |    +0.0 |   +20.9 | 2.444 |
| `1^1F_3` | 2.132 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.119 |
| `1^3F_3` | 2.132 |    +0.0 |   -10.9 |   +14.9 |    +4.0 |    +6.4 |    +0.0 |   +10.4 | 2.156 |
| `1^3F_4` | 2.132 |    +0.0 |   +32.6 |   -44.7 |   -12.1 |    -2.1 |    +0.0 |   -14.2 | 2.118 |
| `1^1G_4` | 2.423 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.408 |
| `1^3G_4` | 2.423 |    +0.0 |    -6.0 |   +10.7 |    +4.7 |    +3.6 |    +0.0 |    +8.3 | 2.446 |
| `1^3G_5` | 2.423 |    +0.0 |   +24.0 |   -42.8 |   -18.8 |    -1.3 |    +0.0 |   -20.1 | 2.403 |

## Fine-Structure Mass Convention (audit note)

Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.

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
| `1^1P_1` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3P_1` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `2^1P_1` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `2^3P_1` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3P_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3P_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_2` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1D_2` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3D_2` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `2^1D_2` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `2^3D_2` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3D_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3D_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3G_3` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1F_3` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3F_3` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3F_4` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1G_4` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3G_4` | 0.220000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3G_5` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |

## Same-J Antisymmetric Spin-Orbit Mixing

Rows below use the mixed eigenvalues from the `(^1L_J, ^3L_J)` mass block. Components are ordered as singlet/triplet in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | theta deg | singlet component | triplet component |
|---|---:|---:|---:|---:|---:|---:|
| `1^1P_1` | 1.384 | 1.382 |    -4.3 |  +19.51 |  +0.943 |  +0.334 |
| `1^3P_1` | 1.394 | 1.396 |    -4.3 |  +19.51 |  +0.334 |  -0.943 |
| `2^1P_1` | 1.922 | 1.918 |    -6.5 |  +28.39 |  +0.880 |  +0.475 |
| `2^3P_1` | 1.930 | 1.934 |    -6.5 |  +28.39 |  +0.475 |  -0.880 |
| `1^1D_2` | 1.796 | 1.787 |   -14.0 |  +32.67 |  +0.842 |  +0.540 |
| `1^3D_2` | 1.809 | 1.818 |   -14.0 |  +32.67 |  +0.540 |  -0.842 |
| `2^1D_2` | 2.244 | 2.236 |   -12.0 |  +34.13 |  +0.828 |  +0.561 |
| `2^3D_2` | 2.254 | 2.262 |   -12.0 |  +34.13 |  +0.561 |  -0.828 |
| `1^1F_3` | 2.132 | 2.119 |   -17.6 |  +36.74 |  +0.801 |  +0.598 |
| `1^3F_3` | 2.142 | 2.156 |   -17.6 |  +36.74 |  +0.598 |  -0.801 |
| `1^1G_4` | 2.423 | 2.408 |   -18.8 |  +38.80 |  +0.779 |  +0.627 |
| `1^3G_4` | 2.431 | 2.446 |   -18.8 |  +38.80 |  +0.627 |  -0.779 |

Mean absolute residual: 10.9 MeV.
Max absolute residual: 43.0 MeV.

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
