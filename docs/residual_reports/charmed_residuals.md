# Residuals: charmed (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia package **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.860 |   -20.2 | high |
| `2^1S_0` | 2.580 | 2.575 |    -5.2 | high |
| `1^3S_1` | 2.040 | 2.039 |    -0.9 | high |
| `2^3S_1` | 2.640 | 2.647 |    +7.4 | high |
| `1^3D_1` | 2.820 | 2.840 |   +19.9 | high |
| `1^3P_0` | 2.400 | 2.412 |   +11.8 | high |
| `1^1P_1` | 2.440 | 2.474 |   +34.2 | high |
| `1^3P_1` | 2.490 | 2.480 |    -9.7 | high |
| `1^3P_2` | 2.500 | 2.486 |   -14.4 | high |
| `1^3D_3` | 2.830 | 2.818 |   -12.2 | high |
| `1^3F_4` | 3.110 | 3.096 |   -14.1 | medium |
| `1^1S_0` | 1.980 | 1.962 |   -18.0 | high |
| `2^1S_0` | 2.670 | 2.667 |    -3.1 | high |
| `1^3S_1` | 2.130 | 2.125 |    -4.7 | high |
| `2^3S_1` | 2.730 | 2.737 |    +6.6 | high |
| `1^3D_1` | 2.900 | 2.908 |    +7.9 | high |
| `1^3P_0` | 2.480 | 2.498 |   +18.2 | high |
| `1^1P_1` | 2.530 | 2.564 |   +34.0 | high |
| `1^3P_1` | 2.570 | 2.569 |    -0.8 | high |
| `1^3P_2` | 2.590 | 2.575 |   -15.2 | high |
| `1^3D_3` | 2.920 | 2.909 |   -10.7 | high |
| `1^3F_4` | 3.190 | 3.187 |    -3.0 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.001 |  -141.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -141.5 | 1.860 |
| `2^1S_0` | 2.629 |   -54.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -54.0 | 2.575 |
| `1^3S_1` | 2.001 |   +37.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +37.8 | 2.039 |
| `2^3S_1` | 2.629 |   +18.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +18.8 | 2.647 |
| `1^3D_1` | 2.831 |    +0.0 |   -84.6 |  +100.9 |   +16.3 |    -7.5 |    +0.0 |    +8.8 | 2.840 |
| `1^3P_0` | 2.475 |    +0.0 |  -148.8 |  +121.0 |   -27.8 |   -35.7 |    +0.0 |   -63.5 | 2.412 |
| `1^1P_1` | 2.475 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.474 |
| `1^3P_1` | 2.475 |    +0.0 |   -74.4 |   +60.5 |   -13.9 |   +17.9 |    +0.0 |    +4.0 | 2.480 |
| `1^3P_2` | 2.475 |    +0.0 |   +74.4 |   -60.5 |   +13.9 |    -3.6 |    +0.0 |   +10.3 | 2.486 |
| `1^3D_3` | 2.831 |    +0.0 |   +56.4 |   -67.3 |   -10.9 |    -2.2 |    +0.0 |   -13.0 | 2.818 |
| `1^3F_4` | 3.124 |    +0.0 |   +41.3 |   -67.9 |   -26.6 |    -1.3 |    +0.0 |   -27.9 | 3.096 |
| `1^1S_0` | 2.091 |  -129.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -129.3 | 1.962 |
| `2^1S_0` | 2.719 |   -52.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -52.4 | 2.667 |
| `1^3S_1` | 2.091 |   +34.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +34.1 | 2.125 |
| `2^3S_1` | 2.719 |   +17.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +17.6 | 2.737 |
| `1^3D_1` | 2.913 |    +0.0 |   -51.5 |   +53.8 |    +2.3 |    -7.6 |    +0.0 |    -5.3 | 2.908 |
| `1^3P_0` | 2.564 |    +0.0 |   -93.6 |   +65.1 |   -28.5 |   -37.6 |    +0.0 |   -66.1 | 2.498 |
| `1^1P_1` | 2.564 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.564 |
| `1^3P_1` | 2.564 |    +0.0 |   -46.8 |   +32.6 |   -14.3 |   +18.8 |    +0.0 |    +4.5 | 2.569 |
| `1^3P_2` | 2.564 |    +0.0 |   +46.8 |   -32.6 |   +14.3 |    -3.8 |    +0.0 |   +10.5 | 2.575 |
| `1^3D_3` | 2.913 |    +0.0 |   +34.4 |   -35.9 |    -1.5 |    -2.2 |    +0.0 |    -3.7 | 2.909 |
| `1^3F_4` | 3.200 |    +0.0 |   +25.0 |   -36.5 |   -11.5 |    -1.3 |    +0.0 |   -12.8 | 3.187 |

## Fine-Structure Mass Convention (audit note)

Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.220000 | `unequal_mass_same_j_mixed` |
| `1^3P_1` | 1.628000 | 0.220000 | `unequal_mass_same_j_mixed` |
| `1^3P_2` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3P_1` | 1.628000 | 0.419000 | `unequal_mass_same_j_mixed` |
| `1^3P_2` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |

## Same-J Antisymmetric Spin-Orbit Mixing

Rows below use the mixed eigenvalues from the `(^1L_J, ^3L_J)` mass block. Components are ordered as singlet/triplet in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | theta deg | singlet component | triplet component |
|---|---:|---:|---:|---:|---:|---:|
| `1^1P_1` | 2.475 | 2.474 |    +2.3 |  -24.82 |  +0.908 |  -0.420 |
| `1^3P_1` | 2.479 | 2.480 |    +2.3 |  -24.82 |  +0.420 |  +0.908 |
| `1^1P_1` | 2.564 | 2.564 |    +1.2 |  -14.15 |  +0.970 |  -0.244 |
| `1^3P_1` | 2.569 | 2.569 |    +1.2 |  -14.15 |  +0.244 |  +0.970 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 2.648 | 2.647 |    -6.9 |  +0.999 |  +0.036 |
| `1^3D_1` | 2.840 | 2.840 |    -6.9 |  +0.036 |  -0.999 |
| `2^3S_1` | 2.737 | 2.737 |    +5.9 |  +0.999 |  -0.035 |
| `1^3D_1` | 2.908 | 2.908 |    +5.9 |  +0.035 |  +0.999 |

Mean absolute residual: 12.4 MeV.
Max absolute residual: 34.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.867 |    -3.8 |
| `1F` | 2 | 3.150 | 3.141 |    -8.5 |
| `1P` | 8 | 2.518 | 2.520 |    +2.3 |
| `1S` | 4 | 2.046 | 2.039 |    -6.9 |
| `2S` | 4 | 2.670 | 2.674 |    +4.2 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
