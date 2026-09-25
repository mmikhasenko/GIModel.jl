# Residuals: strange (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.461 |    -9.1 | high |
| `2^1S_0` | 1.450 | 1.453 |    +3.2 | high |
| `3^1S_0` | 2.020 | 2.013 |    -6.7 | high |
| `1^3S_1` | 0.900 | 0.903 |    +2.7 | high |
| `2^3S_1` | 1.580 | 1.579 |    -1.4 | high |
| `1^3D_1` | 1.780 | 1.774 |    -5.6 | high |
| `3^3S_1` | 2.110 | 2.106 |    -4.2 | high |
| `2^3D_1` | 2.250 | 2.249 |    -1.5 | high |
| `1^3P_0` | 1.240 | 1.219 |   -21.4 | high |
| `2^3P_0` | 1.890 | 1.882 |    -8.5 | high |
| `1^1P_1` | 1.340 | 1.354 |   +13.8 | high |
| `1^3P_1` | 1.380 | 1.384 |    +3.6 | high |
| `2^1P_1` | 1.900 | 1.913 |   +12.6 | high |
| `2^3P_1` | 1.930 | 1.929 |    -0.8 | high |
| `1^3P_2` | 1.430 | 1.420 |    -9.6 | high |
| `2^3P_2` | 1.940 | 1.931 |    -8.7 | high |
| `1^3F_2` | 2.150 | 2.151 |    +0.8 | high |
| `1^1D_2` | 1.780 | 1.794 |   +14.5 | high |
| `1^3D_2` | 1.810 | 1.803 |    -6.9 | high |
| `2^1D_2` | 2.230 | 2.241 |   +10.5 | high |
| `2^3D_2` | 2.260 | 2.256 |    -3.6 | high |
| `1^3D_3` | 1.790 | 1.793 |    +2.7 | high |
| `2^3D_3` | 2.240 | 2.235 |    -4.8 | high |
| `1^3G_3` | 2.460 | 2.458 |    -2.0 | medium |
| `1^1F_3` | 2.120 | 2.131 |   +11.2 | high |
| `1^3F_3` | 2.150 | 2.143 |    -7.0 | high |
| `1^3F_4` | 2.110 | 2.108 |    -2.2 | high |
| `1^1G_4` | 2.410 | 2.422 |   +12.3 | medium |
| `1^3G_4` | 2.440 | 2.433 |    -6.8 | medium |
| `1^3G_5` | 2.390 | 2.388 |    -1.9 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.876 |  -415.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -415.5 | 0.461 |
| `2^1S_0` | 1.536 |   -82.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -82.7 | 1.453 |
| `3^1S_0` | 2.072 |   -59.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -59.0 | 2.013 |
| `1^3S_1` | 0.821 |   +81.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +81.4 | 0.903 |
| `2^3S_1` | 1.545 |   +34.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +34.0 | 1.579 |
| `1^3D_1` | 1.807 |    +0.0 |  -113.9 |   +89.1 |   -24.8 |    -8.7 |    +0.0 |   -33.6 | 1.774 |
| `3^3S_1` | 2.080 |   +26.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +26.1 | 2.106 |
| `2^3D_1` | 2.246 |    +0.0 |   -67.7 |   +74.4 |    +6.7 |    -4.8 |    +0.0 |    +1.9 | 2.249 |
| `1^3P_0` | 1.411 |    +0.0 |  -236.4 |   +93.0 |  -143.4 |   -48.5 |    +0.0 |  -191.9 | 1.219 |
| `2^3P_0` | 1.921 |    +0.0 |   -96.9 |   +76.1 |   -20.8 |   -18.6 |    +0.0 |   -39.4 | 1.882 |
| `1^1P_1` | 1.384 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.354 |
| `1^3P_1` | 1.388 |    +0.0 |  -106.7 |   +49.4 |   -57.2 |   +22.8 |    +0.0 |   -34.4 | 1.384 |
| `2^1P_1` | 1.922 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.913 |
| `2^3P_1` | 1.921 |    +0.0 |   -48.3 |   +37.4 |   -10.8 |    +9.3 |    +0.0 |    -1.6 | 1.929 |
| `1^3P_2` | 1.389 |    +0.0 |   +88.9 |   -53.9 |   +35.1 |    -4.0 |    +0.0 |   +31.1 | 1.420 |
| `2^3P_2` | 1.922 |    +0.0 |   +48.3 |   -36.5 |   +11.7 |    -1.9 |    +0.0 |    +9.8 | 1.931 |
| `1^3F_2` | 2.137 |    +0.0 |   -67.5 |   +84.0 |   +16.5 |    -3.3 |    +0.0 |   +13.3 | 2.151 |
| `1^1D_2` | 1.796 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.794 |
| `1^3D_2` | 1.797 |    +0.0 |   -34.6 |   +31.2 |    -3.3 |    +8.2 |    +0.0 |    +4.9 | 1.803 |
| `2^1D_2` | 2.244 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.241 |
| `2^3D_2` | 2.244 |    +0.0 |   -21.8 |   +25.0 |    +3.1 |    +4.7 |    +0.0 |    +7.8 | 2.256 |
| `1^3D_3` | 1.800 |    +0.0 |   +61.2 |   -66.3 |    -5.2 |    -2.1 |    +0.0 |    -7.3 | 1.793 |
| `2^3D_3` | 2.245 |    +0.0 |   +41.8 |   -50.4 |    -8.6 |    -1.3 |    +0.0 |    -9.9 | 2.235 |
| `1^3G_3` | 2.426 |    +0.0 |   -45.5 |   +79.1 |   +33.6 |    -1.6 |    +0.0 |   +32.0 | 2.458 |
| `1^1F_3` | 2.132 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.131 |
| `1^3F_3` | 2.132 |    +0.0 |   -15.7 |   +21.9 |    +6.2 |    +3.9 |    +0.0 |   +10.0 | 2.143 |
| `1^3F_4` | 2.135 |    +0.0 |   +43.2 |   -68.8 |   -25.6 |    -1.2 |    +0.0 |   -26.8 | 2.108 |
| `1^1G_4` | 2.423 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.422 |
| `1^3G_4` | 2.423 |    +0.0 |    -8.6 |   +16.4 |    +7.8 |    +2.2 |    +0.0 |    +9.9 | 2.433 |
| `1^3G_5` | 2.425 |    +0.0 |   +32.3 |   -68.1 |   -35.8 |    -0.7 |    +0.0 |   -36.6 | 2.388 |

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
| `1^1P_1` | 1.384 | 1.354 |    -1.3 |  +87.46 |  +0.044 |  +0.999 |
| `1^3P_1` | 1.354 | 1.384 |    -1.3 |  +87.46 |  +0.999 |  -0.044 |
| `2^1P_1` | 1.922 | 1.913 |    -8.3 |  +48.04 |  -0.669 |  -0.744 |
| `2^3P_1` | 1.920 | 1.929 |    -8.3 |  +48.04 |  +0.744 |  -0.668 |
| `1^1D_2` | 1.796 | 1.794 |    -3.2 |  +24.38 |  +0.911 |  +0.413 |
| `1^3D_2` | 1.802 | 1.803 |    -3.2 |  +24.38 |  +0.413 |  -0.911 |
| `2^1D_2` | 2.244 | 2.241 |    -6.8 |  +29.86 |  -0.867 |  -0.498 |
| `2^3D_2` | 2.252 | 2.256 |    -6.8 |  +29.86 |  +0.498 |  -0.867 |
| `1^1F_3` | 2.132 | 2.131 |    -2.9 |  +14.60 |  +0.968 |  +0.252 |
| `1^3F_3` | 2.142 | 2.143 |    -2.9 |  +14.60 |  +0.252 |  -0.968 |
| `1^1G_4` | 2.423 | 2.422 |    -2.1 |  +11.33 |  +0.980 |  +0.197 |
| `1^3G_4` | 2.433 | 2.433 |    -2.1 |  +11.33 |  +0.197 |  -0.980 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 1.579 | 1.579 |    +7.8 |  -0.999 |  +0.040 |
| `1^3D_1` | 1.773 | 1.774 |    +7.8 |  +0.040 |  +0.999 |
| `3^3S_1` | 2.106 | 2.106 |    +2.3 |  +1.000 |  -0.016 |
| `2^3D_1` | 2.248 | 2.249 |    +2.3 |  -0.016 |  -1.000 |
| `2^3P_2` | 1.931 | 1.931 |    +1.3 |  -1.000 |  +0.006 |
| `1^3F_2` | 2.151 | 2.151 |    +1.3 |  +0.006 |  +1.000 |
| `2^3D_3` | 2.235 | 2.235 |    +0.3 |  -1.000 |  +0.001 |
| `1^3G_3` | 2.458 | 2.458 |    +0.3 |  +0.001 |  +1.000 |

Mean absolute residual: 6.7 MeV.
Max absolute residual: 21.4 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.791 | 1.793 |    +2.0 |
| `2D` | 4 | 2.244 | 2.244 |    -0.2 |
| `1F` | 4 | 2.130 | 2.130 |    +0.5 |
| `1G` | 4 | 2.421 | 2.421 |    +0.4 |
| `1P` | 4 | 1.379 | 1.378 |    -1.4 |
| `2P` | 4 | 1.923 | 1.922 |    -1.4 |
| `1S` | 2 | 0.792 | 0.792 |    -0.2 |
| `2S` | 2 | 1.548 | 1.547 |    -0.2 |
| `3S` | 2 | 2.087 | 2.083 |    -4.9 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
