# Residuals: strange (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.427 |   -43.0 | high |
| `2^1S_0` | 1.450 | 1.441 |    -8.6 | high |
| `3^1S_0` | 2.020 | 2.003 |   -17.2 | high |
| `1^3S_1` | 0.900 | 0.898 |    -2.3 | high |
| `2^3S_1` | 1.580 | 1.582 |    +2.0 | high |
| `1^3D_1` | 1.780 | 1.780 |    -0.0 | high |
| `3^3S_1` | 2.110 | 2.110 |    +0.1 | high |
| `2^3D_1` | 2.250 | 2.242 |    -8.1 | high |
| `1^3P_0` | 1.240 | 1.245 |    +5.3 | high |
| `2^3P_0` | 1.890 | 1.876 |   -13.6 | high |
| `1^1P_1` | 1.340 | 1.382 |   +42.1 | high |
| `1^3P_1` | 1.380 | 1.394 |   +14.1 | high |
| `2^1P_1` | 1.900 | 1.919 |   +18.5 | high |
| `2^3P_1` | 1.930 | 1.934 |    +4.0 | high |
| `1^3P_2` | 1.430 | 1.401 |   -29.5 | high |
| `2^3P_2` | 1.940 | 1.926 |   -14.2 | high |
| `1^3F_2` | 2.150 | 2.141 |    -8.5 | high |
| `1^1D_2` | 1.780 | 1.787 |    +7.0 | high |
| `1^3D_2` | 1.810 | 1.817 |    +7.3 | high |
| `2^1D_2` | 2.230 | 2.237 |    +6.5 | high |
| `2^3D_2` | 2.260 | 2.262 |    +2.2 | high |
| `1^3D_3` | 1.790 | 1.792 |    +2.1 | high |
| `2^3D_3` | 2.240 | 2.238 |    -1.7 | high |
| `1^3G_3` | 2.460 | 2.443 |   -17.4 | medium |
| `1^1F_3` | 2.120 | 2.119 |    -1.0 | high |
| `1^3F_3` | 2.150 | 2.155 |    +5.4 | high |
| `1^3F_4` | 2.110 | 2.117 |    +6.9 | high |
| `1^1G_4` | 2.410 | 2.408 |    -2.3 | medium |
| `1^3G_4` | 2.440 | 2.446 |    +6.0 | medium |
| `1^3G_5` | 2.390 | 2.402 |   +12.0 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.918 |  -490.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -490.6 | 0.427 |
| `2^1S_0` | 1.529 |   -87.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -87.3 | 1.441 |
| `3^1S_0` | 2.067 |   -64.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -64.2 | 2.003 |
| `1^3S_1` | 0.824 |   +73.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +73.6 | 0.898 |
| `2^3S_1` | 1.546 |   +37.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +37.1 | 1.582 |
| `1^3D_1` | 1.800 |    +0.0 |   -76.7 |   +67.4 |    -9.2 |   -14.2 |    +0.0 |   -23.4 | 1.780 |
| `3^3S_1` | 2.080 |   +30.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +30.2 | 2.110 |
| `2^3D_1` | 2.245 |    +0.0 |   -47.0 |   +50.9 |    +3.9 |    -8.0 |    +0.0 |    -4.1 | 2.242 |
| `1^3P_0` | 1.397 |    +0.0 |  -157.8 |   +85.3 |   -72.6 |   -79.2 |    +0.0 |  -151.8 | 1.245 |
| `2^3P_0` | 1.921 |    +0.0 |   -67.7 |   +54.3 |   -13.5 |   -31.2 |    +0.0 |   -44.7 | 1.876 |
| `1^1P_1` | 1.384 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.382 |
| `1^3P_1` | 1.385 |    +0.0 |   -71.6 |   +42.9 |   -28.6 |   +37.2 |    +0.0 |    +8.6 | 1.394 |
| `2^1P_1` | 1.922 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.919 |
| `2^3P_1` | 1.922 |    +0.0 |   -33.8 |   +26.7 |    -7.1 |   +15.6 |    +0.0 |    +8.5 | 1.934 |
| `1^3P_2` | 1.385 |    +0.0 |   +65.3 |   -43.2 |   +22.1 |    -7.0 |    +0.0 |   +15.1 | 1.401 |
| `2^3P_2` | 1.922 |    +0.0 |   +33.9 |   -26.5 |    +7.4 |    -3.2 |    +0.0 |    +4.2 | 1.926 |
| `1^3F_2` | 2.134 |    +0.0 |   -46.0 |   +58.4 |   +12.4 |    -5.3 |    +0.0 |    +7.0 | 2.141 |
| `1^1D_2` | 1.796 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.787 |
| `1^3D_2` | 1.796 |    +0.0 |   -23.9 |   +22.8 |    -1.0 |   +13.5 |    +0.0 |   +12.5 | 1.817 |
| `2^1D_2` | 2.244 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.237 |
| `2^3D_2` | 2.244 |    +0.0 |   -15.3 |   +16.9 |    +1.7 |    +7.8 |    +0.0 |    +9.5 | 2.262 |
| `1^3D_3` | 1.797 |    +0.0 |   +44.7 |   -46.3 |    -1.6 |    -3.7 |    +0.0 |    -5.3 | 1.792 |
| `2^3D_3` | 2.245 |    +0.0 |   +29.9 |   -33.9 |    -4.0 |    -2.2 |    +0.0 |    -6.2 | 2.238 |
| `1^3G_3` | 2.424 |    +0.0 |   -31.3 |   +52.5 |   +21.2 |    -2.6 |    +0.0 |   +18.6 | 2.443 |
| `1^1F_3` | 2.132 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.119 |
| `1^3F_3` | 2.132 |    +0.0 |   -11.0 |   +14.8 |    +3.9 |    +6.4 |    +0.0 |   +10.3 | 2.155 |
| `1^3F_4` | 2.133 |    +0.0 |   +31.4 |   -45.3 |   -13.9 |    -2.1 |    +0.0 |   -16.0 | 2.117 |
| `1^1G_4` | 2.423 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.408 |
| `1^3G_4` | 2.423 |    +0.0 |    -6.0 |   +10.7 |    +4.6 |    +3.6 |    +0.0 |    +8.2 | 2.446 |
| `1^3G_5` | 2.423 |    +0.0 |   +23.3 |   -43.4 |   -20.1 |    -1.3 |    +0.0 |   -21.4 | 2.402 |

## Fine-Structure Mass Convention (audit note)

Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 0.220000 | 0.419000 | `disabled` |
| `2^1S_0` | 0.220000 | 0.419000 | `disabled` |
| `3^1S_0` | 0.220000 | 0.419000 | `disabled` |
| `1^3S_1` | 0.220000 | 0.419000 | `disabled` |
| `2^3S_1` | 0.220000 | 0.419000 | `disabled` |
| `1^3D_1` | 0.220000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `3^3S_1` | 0.220000 | 0.419000 | `disabled` |
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
| `1^1P_1` | 1.384 | 1.382 |    -3.6 |  +18.52 |  +0.948 |  +0.318 |
| `1^3P_1` | 1.393 | 1.394 |    -3.6 |  +18.52 |  +0.318 |  -0.948 |
| `2^1P_1` | 1.922 | 1.919 |    -6.5 |  +28.39 |  -0.879 |  -0.475 |
| `2^3P_1` | 1.930 | 1.934 |    -6.5 |  +28.39 |  +0.475 |  -0.879 |
| `1^1D_2` | 1.796 | 1.787 |   -13.7 |  +32.63 |  +0.842 |  +0.539 |
| `1^3D_2` | 1.809 | 1.817 |   -13.7 |  +32.63 |  +0.539 |  -0.842 |
| `2^1D_2` | 2.244 | 2.237 |   -11.9 |  +34.12 |  -0.828 |  -0.561 |
| `2^3D_2` | 2.254 | 2.262 |   -11.9 |  +34.12 |  +0.561 |  -0.828 |
| `1^1F_3` | 2.132 | 2.119 |   -17.5 |  +36.73 |  +0.801 |  +0.598 |
| `1^3F_3` | 2.142 | 2.155 |   -17.5 |  +36.73 |  +0.598 |  -0.801 |
| `1^1G_4` | 2.423 | 2.408 |   -18.7 |  +38.79 |  +0.779 |  +0.627 |
| `1^3G_4` | 2.431 | 2.446 |   -18.7 |  +38.79 |  +0.627 |  -0.779 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 1.583 | 1.582 |   +13.3 |  -0.998 |  +0.067 |
| `1^3D_1` | 1.777 | 1.780 |   +13.3 |  +0.067 |  +0.996 |
| `3^3S_1` | 2.110 | 2.110 |    +3.9 |  +1.000 |  -0.030 |
| `2^3D_1` | 2.241 | 2.242 |    +3.9 |  -0.030 |  -0.999 |
| `2^3P_2` | 1.926 | 1.926 |    +3.3 |  -1.000 |  +0.015 |
| `1^3F_2` | 2.141 | 2.141 |    +3.3 |  +0.015 |  +1.000 |
| `2^3D_3` | 2.238 | 2.238 |    +1.0 |  -1.000 |  +0.005 |
| `1^3G_3` | 2.442 | 2.443 |    +1.0 |  +0.005 |  +1.000 |

Mean absolute residual: 10.3 MeV.
Max absolute residual: 43.0 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.791 | 1.795 |    +4.3 |
| `2D` | 4 | 2.244 | 2.244 |    +0.4 |
| `1F` | 4 | 2.130 | 2.131 |    +1.8 |
| `1G` | 4 | 2.421 | 2.422 |    +1.2 |
| `1P` | 4 | 1.379 | 1.381 |    +2.2 |
| `2P` | 4 | 1.923 | 1.922 |    -1.4 |
| `1S` | 2 | 0.792 | 0.780 |   -12.5 |
| `2S` | 2 | 1.548 | 1.547 |    -0.6 |
| `3S` | 2 | 2.087 | 2.083 |    -4.2 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
