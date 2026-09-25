# Residuals: charmed (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.873 |    -6.6 | high |
| `2^1S_0` | 2.580 | 2.582 |    +2.4 | high |
| `1^3S_1` | 2.040 | 2.038 |    -2.0 | high |
| `2^3S_1` | 2.640 | 2.645 |    +4.9 | high |
| `1^3D_1` | 2.820 | 2.815 |    -5.1 | high |
| `1^3P_0` | 2.400 | 2.387 |   -13.0 | high |
| `1^1P_1` | 2.440 | 2.456 |   +16.0 | high |
| `1^3P_1` | 2.490 | 2.476 |   -13.8 | high |
| `1^3P_2` | 2.500 | 2.497 |    -2.7 | high |
| `1^3D_3` | 2.830 | 2.832 |    +1.9 | high |
| `1^3F_4` | 3.110 | 3.113 |    +2.7 | medium |
| `1^1S_0` | 1.980 | 1.974 |    -6.0 | high |
| `2^1S_0` | 2.670 | 2.675 |    +4.6 | high |
| `1^3S_1` | 2.130 | 2.125 |    -5.3 | high |
| `2^3S_1` | 2.730 | 2.734 |    +4.4 | high |
| `1^3D_1` | 2.900 | 2.897 |    -2.7 | high |
| `1^3P_0` | 2.480 | 2.474 |    -5.7 | high |
| `1^1P_1` | 2.530 | 2.546 |   +15.7 | high |
| `1^3P_1` | 2.570 | 2.565 |    -4.9 | high |
| `1^3P_2` | 2.590 | 2.587 |    -2.7 | high |
| `1^3D_3` | 2.920 | 2.915 |    -4.8 | high |
| `1^3F_4` | 3.190 | 3.190 |    -0.4 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.014 |  -141.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -141.0 | 1.873 |
| `2^1S_0` | 2.627 |   -44.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -44.4 | 2.582 |
| `1^3S_1` | 2.003 |   +35.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +35.3 | 2.038 |
| `2^3S_1` | 2.629 |   +16.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +16.2 | 2.645 |
| `1^3D_1` | 2.836 |    +0.0 |   -70.7 |   +54.0 |   -16.7 |    -4.9 |    +0.0 |   -21.6 | 2.815 |
| `1^3P_0` | 2.488 |    +0.0 |  -135.2 |   +57.7 |   -77.5 |   -23.8 |    +0.0 |  -101.3 | 2.387 |
| `1^1P_1` | 2.475 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.456 |
| `1^3P_1` | 2.477 |    +0.0 |   -61.9 |   +30.2 |   -31.7 |   +11.2 |    +0.0 |   -20.6 | 2.476 |
| `1^3P_2` | 2.478 |    +0.0 |   +53.7 |   -32.2 |   +21.4 |    -2.0 |    +0.0 |   +19.4 | 2.497 |
| `1^3D_3` | 2.833 |    +0.0 |   +39.2 |   -39.0 |    +0.3 |    -1.2 |    +0.0 |    -0.9 | 2.832 |
| `1^3F_4` | 3.125 |    +0.0 |   +28.6 |   -40.3 |   -11.6 |    -0.7 |    +0.0 |   -12.4 | 3.113 |
| `1^1S_0` | 2.104 |  -130.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -130.3 | 1.974 |
| `2^1S_0` | 2.718 |   -43.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -43.3 | 2.675 |
| `1^3S_1` | 2.093 |   +32.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +32.1 | 2.125 |
| `2^3S_1` | 2.719 |   +15.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +15.2 | 2.734 |
| `1^3D_1` | 2.917 |    +0.0 |   -64.0 |   +48.7 |   -15.3 |    -5.0 |    +0.0 |   -20.3 | 2.897 |
| `1^3P_0` | 2.577 |    +0.0 |  -126.6 |   +49.3 |   -77.4 |   -25.7 |    +0.0 |  -103.1 | 2.474 |
| `1^1P_1` | 2.564 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.546 |
| `1^3P_1` | 2.566 |    +0.0 |   -56.7 |   +25.2 |   -31.5 |   +11.8 |    +0.0 |   -19.6 | 2.565 |
| `1^3P_2` | 2.566 |    +0.0 |   +48.8 |   -25.9 |   +22.9 |    -2.1 |    +0.0 |   +20.8 | 2.587 |
| `1^3D_3` | 2.914 |    +0.0 |   +35.6 |   -33.6 |    +2.0 |    -1.2 |    +0.0 |    +0.8 | 2.915 |
| `1^3F_4` | 3.201 |    +0.0 |   +26.1 |   -36.4 |   -10.3 |    -0.7 |    +0.0 |   -11.1 | 3.190 |

## Fine-Structure Mass Convention (audit note)

Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 1.628000 | 0.220000 | `disabled` |
| `2^1S_0` | 1.628000 | 0.220000 | `disabled` |
| `1^3S_1` | 1.628000 | 0.220000 | `disabled` |
| `2^3S_1` | 1.628000 | 0.220000 | `disabled` |
| `1^3D_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.220000 | `unequal_mass_same_j_mixed` |
| `1^3P_1` | 1.628000 | 0.220000 | `unequal_mass_same_j_mixed` |
| `1^3P_2` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 1.628000 | 0.419000 | `disabled` |
| `2^1S_0` | 1.628000 | 0.419000 | `disabled` |
| `1^3S_1` | 1.628000 | 0.419000 | `disabled` |
| `2^3S_1` | 1.628000 | 0.419000 | `disabled` |
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
| `1^1P_1` | 2.475 | 2.456 |    +4.2 |  -77.73 |  +0.212 |  -0.977 |
| `1^3P_1` | 2.457 | 2.476 |    +4.2 |  -77.73 |  +0.977 |  +0.212 |
| `1^1P_1` | 2.564 | 2.546 |    +3.7 |  -78.79 |  +0.194 |  -0.981 |
| `1^3P_1` | 2.546 | 2.565 |    +3.7 |  -78.79 |  +0.981 |  +0.194 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 2.645 | 2.645 |    +3.3 |  -1.000 |  +0.020 |
| `1^3D_1` | 2.815 | 2.815 |    +3.3 |  +0.020 |  +1.000 |
| `2^3S_1` | 2.734 | 2.734 |    +2.7 |  -1.000 |  +0.017 |
| `1^3D_1` | 2.897 | 2.897 |    +2.7 |  +0.017 |  +1.000 |

Mean absolute residual: 5.8 MeV.
Max absolute residual: 16.0 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.868 |    -2.2 |
| `1F` | 2 | 3.150 | 3.151 |    +1.1 |
| `1P` | 8 | 2.518 | 2.517 |    -0.3 |
| `1S` | 4 | 2.046 | 2.042 |    -4.3 |
| `2S` | 4 | 2.670 | 2.674 |    +4.3 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
