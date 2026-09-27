# Residuals: charmed (GI-style)

Model: with GI momentum-sandwiched smeared contact hyperfine (every L), fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.873 |    -6.6 | high |
| `2^1S_0` | 2.580 | 2.582 |    +2.4 | high |
| `1^3S_1` | 2.040 | 2.038 |    -2.0 | high |
| `2^3S_1` | 2.640 | 2.645 |    +4.9 | high |
| `1^3D_1` | 2.820 | 2.816 |    -3.5 | high |
| `1^3P_0` | 2.400 | 2.395 |    -5.2 | high |
| `1^1P_1` | 2.440 | 2.456 |   +15.6 | high |
| `1^3P_1` | 2.490 | 2.465 |   -24.8 | high |
| `1^3P_2` | 2.500 | 2.502 |    +2.0 | high |
| `1^3D_3` | 2.830 | 2.833 |    +2.9 | high |
| `1^3F_4` | 3.110 | 3.113 |    +3.0 | medium |
| `1^1S_0` | 1.980 | 1.974 |    -6.0 | high |
| `2^1S_0` | 2.670 | 2.675 |    +4.6 | high |
| `1^3S_1` | 2.130 | 2.125 |    -5.3 | high |
| `2^3S_1` | 2.730 | 2.734 |    +4.4 | high |
| `1^3D_1` | 2.900 | 2.899 |    -1.5 | high |
| `1^3P_0` | 2.480 | 2.481 |    +0.6 | high |
| `1^1P_1` | 2.530 | 2.547 |   +17.4 | high |
| `1^3P_1` | 2.570 | 2.554 |   -15.6 | high |
| `1^3P_2` | 2.590 | 2.591 |    +1.1 | high |
| `1^3D_3` | 2.920 | 2.916 |    -3.9 | high |
| `1^3F_4` | 3.190 | 3.190 |    -0.2 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.014 |  -141.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -141.0 | 1.873 |
| `2^1S_0` | 2.627 |   -44.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -44.4 | 2.582 |
| `1^3S_1` | 2.003 |   +35.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +35.3 | 2.038 |
| `2^3S_1` | 2.629 |   +16.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +16.2 | 2.645 |
| `1^3D_1` | 2.836 |    +1.6 |   -70.4 |   +54.1 |   -16.3 |    -4.9 |    +0.0 |   -19.6 | 2.816 |
| `1^3P_0` | 2.486 |    +7.7 |  -133.4 |   +58.0 |   -75.4 |   -23.6 |    +0.0 |   -91.3 | 2.395 |
| `1^1P_1` | 2.476 |   -18.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -18.8 | 2.456 |
| `1^3P_1` | 2.477 |    +6.3 |   -61.0 |   +30.3 |   -30.7 |   +11.1 |    +0.0 |   -13.3 | 2.465 |
| `1^3P_2` | 2.479 |    +4.6 |   +52.9 |   -32.4 |   +20.6 |    -2.0 |    +0.0 |   +23.2 | 2.502 |
| `1^3D_3` | 2.833 |    +1.0 |   +39.1 |   -39.0 |    +0.1 |    -1.2 |    +0.0 |    -0.1 | 2.833 |
| `1^3F_4` | 3.125 |    +0.3 |   +28.6 |   -40.3 |   -11.7 |    -0.7 |    +0.0 |   -12.1 | 3.113 |
| `1^1S_0` | 2.104 |  -130.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -130.3 | 1.974 |
| `2^1S_0` | 2.718 |   -43.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -43.3 | 2.675 |
| `1^3S_1` | 2.093 |   +32.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +32.1 | 2.125 |
| `2^3S_1` | 2.719 |   +15.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +15.2 | 2.734 |
| `1^3D_1` | 2.917 |    +1.2 |   -63.7 |   +48.7 |   -15.0 |    -5.0 |    +0.0 |   -18.8 | 2.899 |
| `1^3P_0` | 2.576 |    +6.2 |  -125.0 |   +49.4 |   -75.6 |   -25.5 |    +0.0 |   -94.9 | 2.481 |
| `1^1P_1` | 2.565 |   -14.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -14.7 | 2.547 |
| `1^3P_1` | 2.566 |    +5.0 |   -55.9 |   +25.3 |   -30.7 |   +11.7 |    +0.0 |   -13.9 | 2.554 |
| `1^3P_2` | 2.567 |    +3.7 |   +48.2 |   -25.9 |   +22.3 |    -2.1 |    +0.0 |   +23.9 | 2.591 |
| `1^3D_3` | 2.915 |    +0.8 |   +35.5 |   -33.6 |    +1.9 |    -1.2 |    +0.0 |    +1.5 | 2.916 |
| `1^3F_4` | 3.201 |    +0.3 |   +26.1 |   -36.4 |   -10.4 |    -0.7 |    +0.0 |   -10.8 | 3.190 |

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
| `1^1P_1` | 2.457 | 2.456 |    +3.7 |  -25.77 |  +0.901 |  -0.435 |
| `1^3P_1` | 2.463 | 2.465 |    +3.7 |  -25.77 |  +0.435 |  +0.901 |
| `1^1P_1` | 2.550 | 2.547 |    +3.5 |  -39.57 |  +0.771 |  -0.637 |
| `1^3P_1` | 2.552 | 2.554 |    +3.5 |  -39.57 |  +0.637 |  +0.771 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 2.645 | 2.645 |    +3.4 |  -1.000 |  +0.020 |
| `1^3D_1` | 2.816 | 2.816 |    +3.4 |  +0.020 |  +1.000 |
| `2^3S_1` | 2.734 | 2.734 |    +2.8 |  -1.000 |  +0.017 |
| `1^3D_1` | 2.898 | 2.899 |    +2.8 |  +0.017 |  +1.000 |

Mean absolute residual: 6.1 MeV.
Max absolute residual: 24.8 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.869 |    -1.1 |
| `1F` | 2 | 3.150 | 3.151 |    +1.4 |
| `1P` | 8 | 2.518 | 2.517 |    -0.5 |
| `1S` | 4 | 2.046 | 2.042 |    -4.3 |
| `2S` | 4 | 2.670 | 2.674 |    +4.3 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
