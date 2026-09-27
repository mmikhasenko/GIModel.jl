# Residuals: strange (GI-style)

Model: with GI momentum-sandwiched smeared contact hyperfine (every L), fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.461 |    -9.1 | high |
| `2^1S_0` | 1.450 | 1.453 |    +3.2 | high |
| `3^1S_0` | 2.020 | 2.013 |    -6.7 | high |
| `1^3S_1` | 0.900 | 0.903 |    +2.7 | high |
| `2^3S_1` | 1.580 | 1.579 |    -1.4 | high |
| `1^3D_1` | 1.780 | 1.777 |    -3.0 | high |
| `3^3S_1` | 2.110 | 2.106 |    -4.2 | high |
| `2^3D_1` | 2.250 | 2.251 |    +0.8 | high |
| `1^3P_0` | 1.240 | 1.233 |    -6.5 | high |
| `2^3P_0` | 1.890 | 1.890 |    -0.0 | high |
| `1^1P_1` | 1.340 | 1.352 |   +11.5 | high |
| `1^3P_1` | 1.380 | 1.366 |   -14.4 | high |
| `2^1P_1` | 1.900 | 1.895 |    -4.9 | high |
| `2^3P_1` | 1.930 | 1.930 |    +0.1 | high |
| `1^3P_2` | 1.430 | 1.428 |    -2.0 | high |
| `2^3P_2` | 1.940 | 1.938 |    -1.7 | high |
| `1^3F_2` | 2.150 | 2.151 |    +1.4 | high |
| `1^1D_2` | 1.780 | 1.790 |    +9.6 | high |
| `1^3D_2` | 1.810 | 1.804 |    -5.6 | high |
| `2^1D_2` | 2.230 | 2.236 |    +5.9 | high |
| `2^3D_2` | 2.260 | 2.257 |    -3.0 | high |
| `1^3D_3` | 1.790 | 1.794 |    +4.1 | high |
| `2^3D_3` | 2.240 | 2.237 |    -3.0 | high |
| `1^3G_3` | 2.460 | 2.458 |    -1.8 | medium |
| `1^1F_3` | 2.120 | 2.130 |    +9.9 | high |
| `1^3F_3` | 2.150 | 2.143 |    -6.6 | high |
| `1^3F_4` | 2.110 | 2.108 |    -1.8 | high |
| `1^1G_4` | 2.410 | 2.422 |   +11.9 | medium |
| `1^3G_4` | 2.440 | 2.433 |    -6.7 | medium |
| `1^3G_5` | 2.390 | 2.388 |    -1.8 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 0.876 |  -415.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -415.5 | 0.461 |
| `2^1S_0` | 1.536 |   -82.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -82.7 | 1.453 |
| `3^1S_0` | 2.072 |   -59.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -59.0 | 2.013 |
| `1^3S_1` | 0.821 |   +81.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +81.4 | 0.903 |
| `2^3S_1` | 1.545 |   +34.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +34.0 | 1.579 |
| `1^3D_1` | 1.806 |    +2.5 |  -113.3 |   +89.3 |   -24.0 |    -8.7 |    +0.0 |   -30.2 | 1.777 |
| `3^3S_1` | 2.080 |   +26.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +26.1 | 2.106 |
| `2^3D_1` | 2.246 |    +2.3 |   -67.4 |   +74.4 |    +7.0 |    -4.8 |    +0.0 |    +4.5 | 2.251 |
| `1^3P_0` | 1.405 |   +14.5 |  -232.3 |   +93.9 |  -138.4 |   -48.1 |    +0.0 |  -171.9 | 1.233 |
| `2^3P_0` | 1.920 |    +8.3 |   -95.9 |   +75.8 |   -20.1 |   -18.5 |    +0.0 |   -30.2 | 1.890 |
| `1^1P_1` | 1.386 |   -34.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -34.6 | 1.352 |
| `1^3P_1` | 1.387 |   +11.4 |  -104.6 |   +49.8 |   -54.8 |   +22.5 |    +0.0 |   -20.9 | 1.366 |
| `2^1P_1` | 1.923 |   -26.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -26.1 | 1.895 |
| `2^3P_1` | 1.921 |    +7.8 |   -47.7 |   +37.3 |   -10.4 |    +9.2 |    +0.0 |    +6.6 | 1.930 |
| `1^3P_2` | 1.391 |    +7.4 |   +87.3 |   -54.1 |   +33.2 |    -4.0 |    +0.0 |   +36.7 | 1.428 |
| `2^3P_2` | 1.922 |    +6.9 |   +47.6 |   -36.5 |   +11.1 |    -1.9 |    +0.0 |   +16.1 | 1.938 |
| `1^3F_2` | 2.137 |    +0.6 |   -67.3 |   +84.1 |   +16.7 |    -3.3 |    +0.0 |   +14.1 | 2.151 |
| `1^1D_2` | 1.796 |    -5.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -5.7 | 1.790 |
| `1^3D_2` | 1.797 |    +2.0 |   -34.4 |   +31.3 |    -3.1 |    +8.2 |    +0.0 |    +7.0 | 1.804 |
| `2^1D_2` | 2.244 |    -6.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -6.2 | 2.236 |
| `2^3D_2` | 2.244 |    +2.0 |   -21.7 |   +25.0 |    +3.2 |    +4.7 |    +0.0 |   +10.0 | 2.257 |
| `1^3D_3` | 1.800 |    +1.4 |   +60.9 |   -66.4 |    -5.5 |    -2.1 |    +0.0 |    -6.2 | 1.794 |
| `2^3D_3` | 2.245 |    +1.8 |   +41.6 |   -50.4 |    -8.8 |    -1.3 |    +0.0 |    -8.3 | 2.237 |
| `1^3G_3` | 2.426 |    +0.2 |   -45.5 |   +79.1 |   +33.6 |    -1.6 |    +0.0 |   +32.2 | 2.458 |
| `1^1F_3` | 2.132 |    -1.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -1.4 | 2.130 |
| `1^3F_3` | 2.132 |    +0.5 |   -15.7 |   +21.9 |    +6.2 |    +3.9 |    +0.0 |   +10.6 | 2.143 |
| `1^3F_4` | 2.135 |    +0.4 |   +43.2 |   -68.8 |   -25.7 |    -1.2 |    +0.0 |   -26.5 | 2.108 |
| `1^1G_4` | 2.423 |    -0.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    -0.4 | 2.422 |
| `1^3G_4` | 2.423 |    +0.1 |    -8.6 |   +16.4 |    +7.8 |    +2.2 |    +0.0 |   +10.1 | 2.433 |
| `1^3G_5` | 2.425 |    +0.1 |   +32.3 |   -68.1 |   -35.9 |    -0.7 |    +0.0 |   -36.5 | 2.388 |

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
| `1^1P_1` | 1.352 | 1.352 |    -1.1 |   +4.31 |  +0.997 |  +0.075 |
| `1^3P_1` | 1.366 | 1.366 |    -1.1 |   +4.31 |  +0.075 |  -0.997 |
| `2^1P_1` | 1.897 | 1.895 |    -8.3 |  +14.18 |  -0.969 |  -0.245 |
| `2^3P_1` | 1.928 | 1.930 |    -8.3 |  +14.18 |  +0.245 |  -0.969 |
| `1^1D_2` | 1.790 | 1.790 |    -3.2 |  +12.84 |  +0.975 |  +0.222 |
| `1^3D_2` | 1.804 | 1.804 |    -3.2 |  +12.84 |  +0.222 |  -0.975 |
| `2^1D_2` | 2.238 | 2.236 |    -6.8 |  +20.14 |  -0.939 |  -0.344 |
| `2^3D_2` | 2.254 | 2.257 |    -6.8 |  +20.14 |  +0.344 |  -0.939 |
| `1^1F_3` | 2.131 | 2.130 |    -2.9 |  +12.60 |  +0.976 |  +0.218 |
| `1^3F_3` | 2.143 | 2.143 |    -2.9 |  +12.60 |  +0.218 |  -0.976 |
| `1^1G_4` | 2.422 | 2.422 |    -2.1 |  +10.82 |  +0.982 |  +0.188 |
| `1^3G_4` | 2.433 | 2.433 |    -2.1 |  +10.82 |  +0.188 |  -0.982 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 1.579 | 1.579 |    +7.8 |  -0.999 |  +0.039 |
| `1^3D_1` | 1.776 | 1.777 |    +7.8 |  +0.039 |  +0.999 |
| `3^3S_1` | 2.106 | 2.106 |    +2.3 |  +1.000 |  -0.016 |
| `2^3D_1` | 2.251 | 2.251 |    +2.3 |  -0.016 |  -1.000 |
| `2^3P_2` | 1.938 | 1.938 |    +1.2 |  -1.000 |  +0.006 |
| `1^3F_2` | 2.151 | 2.151 |    +1.2 |  +0.006 |  +1.000 |
| `2^3D_3` | 2.237 | 2.237 |    +0.2 |  -1.000 |  +0.001 |
| `1^3G_3` | 2.458 | 2.458 |    +0.2 |  +0.001 |  +1.000 |

Mean absolute residual: 4.8 MeV.
Max absolute residual: 14.4 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.791 | 1.793 |    +2.0 |
| `2D` | 4 | 2.244 | 2.244 |    -0.2 |
| `1F` | 4 | 2.130 | 2.130 |    +0.5 |
| `1G` | 4 | 2.421 | 2.421 |    +0.4 |
| `1P` | 4 | 1.379 | 1.377 |    -2.1 |
| `2P` | 4 | 1.923 | 1.921 |    -1.9 |
| `1S` | 2 | 0.792 | 0.792 |    -0.2 |
| `2S` | 2 | 1.548 | 1.547 |    -0.2 |
| `3S` | 2 | 2.087 | 2.083 |    -4.9 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
