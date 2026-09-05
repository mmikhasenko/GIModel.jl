# Residuals: charmed (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05333), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.860 |   -20.2 | high |
| `2^1S_0` | 2.580 | 2.575 |    -5.2 | high |
| `1^3S_1` | 2.040 | 2.039 |    -0.9 | high |
| `2^3S_1` | 2.640 | 2.647 |    +7.5 | high |
| `1^3D_1` | 2.820 | 2.835 |   +15.4 | high |
| `1^3P_0` | 2.400 | 2.398 |    -1.7 | high |
| `1^1P_1` | 2.440 | 2.475 |   +35.2 | high |
| `1^3P_1` | 2.490 | 2.478 |   -12.1 | high |
| `1^3P_2` | 2.500 | 2.484 |   -16.4 | high |
| `1^3D_3` | 2.830 | 2.816 |   -13.7 | high |
| `1^3F_4` | 3.110 | 3.095 |   -15.0 | medium |
| `1^1S_0` | 1.980 | 1.962 |   -18.0 | high |
| `2^1S_0` | 2.670 | 2.667 |    -3.1 | high |
| `1^3S_1` | 2.130 | 2.125 |    -4.7 | high |
| `2^3S_1` | 2.730 | 2.737 |    +6.7 | high |
| `1^3D_1` | 2.900 | 2.906 |    +6.4 | high |
| `1^3P_0` | 2.480 | 2.490 |   +10.2 | high |
| `1^1P_1` | 2.530 | 2.564 |   +34.2 | high |
| `1^3P_1` | 2.570 | 2.569 |    -1.3 | high |
| `1^3P_2` | 2.590 | 2.574 |   -16.0 | high |
| `1^3D_3` | 2.920 | 2.909 |   -11.2 | high |
| `1^3F_4` | 3.190 | 3.187 |    -3.3 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.023 |  -163.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -163.2 | 1.860 |
| `2^1S_0` | 2.625 |   -50.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -50.7 | 2.575 |
| `1^3S_1` | 2.003 |   +35.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +35.6 | 2.039 |
| `2^3S_1` | 2.629 |   +18.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +18.8 | 2.647 |
| `1^3D_1` | 2.836 |    +0.0 |   -95.7 |  +102.6 |    +6.9 |    -8.3 |    +0.0 |    -1.4 | 2.835 |
| `1^3P_0` | 2.489 |    +0.0 |  -176.4 |  +126.0 |   -50.3 |   -40.0 |    +0.0 |   -90.4 | 2.398 |
| `1^1P_1` | 2.475 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.475 |
| `1^3P_1` | 2.477 |    +0.0 |   -79.0 |   +61.4 |   -17.6 |   +18.6 |    +0.0 |    +1.0 | 2.478 |
| `1^3P_2` | 2.477 |    +0.0 |   +69.1 |   -59.4 |    +9.7 |    -3.4 |    +0.0 |    +6.3 | 2.484 |
| `1^3D_3` | 2.832 |    +0.0 |   +52.6 |   -66.6 |   -14.0 |    -2.0 |    +0.0 |   -16.1 | 2.816 |
| `1^3F_4` | 3.125 |    +0.0 |   +39.1 |   -67.7 |   -28.6 |    -1.2 |    +0.0 |   -29.8 | 3.095 |
| `1^1S_0` | 2.113 |  -150.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -150.8 | 1.962 |
| `2^1S_0` | 2.717 |   -50.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -50.4 | 2.667 |
| `1^3S_1` | 2.093 |   +32.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +32.1 | 2.125 |
| `2^3S_1` | 2.719 |   +17.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +17.4 | 2.737 |
| `1^3D_1` | 2.915 |    +0.0 |   -56.1 |   +54.9 |    -1.2 |    -8.2 |    +0.0 |    -9.4 | 2.906 |
| `1^3P_0` | 2.573 |    +0.0 |  -108.9 |   +68.7 |   -40.2 |   -42.1 |    +0.0 |   -82.3 | 2.490 |
| `1^1P_1` | 2.564 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.564 |
| `1^3P_1` | 2.565 |    +0.0 |   -48.4 |   +33.0 |   -15.5 |   +19.3 |    +0.0 |    +3.8 | 2.569 |
| `1^3P_2` | 2.565 |    +0.0 |   +44.4 |   -32.0 |   +12.5 |    -3.6 |    +0.0 |    +8.9 | 2.574 |
| `1^3D_3` | 2.914 |    +0.0 |   +32.9 |   -35.5 |    -2.6 |    -2.1 |    +0.0 |    -4.7 | 2.909 |
| `1^3F_4` | 3.200 |    +0.0 |   +24.2 |   -36.3 |   -12.1 |    -1.3 |    +0.0 |   -13.3 | 3.187 |

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
| `1^1P_1` | 2.475 | 2.475 |    +0.5 |  -10.46 |  +0.983 |  -0.181 |
| `1^3P_1` | 2.478 | 2.478 |    +0.5 |  -10.46 |  +0.181 |  +0.983 |
| `1^1P_1` | 2.564 | 2.564 |    +0.8 |  -10.45 |  +0.983 |  -0.181 |
| `1^3P_1` | 2.569 | 2.569 |    +0.8 |  -10.45 |  +0.181 |  +0.983 |

## Same-J Tensor Mixing

Rows below use the mixed eigenvalues from triplet `L=J-1` / `L=J+1` tensor blocks. Components are ordered as lower-`L`/higher-`L` in the unmixed basis.

| state | unmixed GeV | mixed GeV | offdiag MeV | lower-L component | higher-L component |
|---|---:|---:|---:|---:|---:|
| `2^3S_1` | 2.648 | 2.647 |    +5.4 |  -1.000 |  +0.029 |
| `1^3D_1` | 2.835 | 2.835 |    +5.4 |  +0.029 |  +0.999 |
| `2^3S_1` | 2.737 | 2.737 |    +4.7 |  -1.000 |  +0.028 |
| `1^3D_1` | 2.906 | 2.906 |    +4.7 |  +0.028 |  +0.999 |

Mean absolute residual: 11.8 MeV.
Max absolute residual: 35.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.865 |    -5.4 |
| `1F` | 2 | 3.150 | 3.141 |    -9.2 |
| `1P` | 8 | 2.518 | 2.518 |    +0.6 |
| `1S` | 4 | 2.046 | 2.039 |    -6.9 |
| `2S` | 4 | 2.670 | 2.674 |    +4.3 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
