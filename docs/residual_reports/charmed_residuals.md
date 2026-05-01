# Residuals: charmed (GI-style)

Model: finite-difference + `relativistic` kinetic, with GI momentum-sandwiched smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, see `src/GIModel.jl` and sibling sources under `src/` (Julia **GIModel**, `Project.toml`).

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.800 |   -79.5 | high |
| `2^1S_0` | 2.580 | 2.504 |   -76.0 | high |
| `1^3S_1` | 2.040 | 2.068 |   +28.2 | high |
| `2^3S_1` | 2.640 | 2.670 |   +30.4 | high |
| `1^3D_1` | 2.820 | 2.878 |   +57.8 | high |
| `1^3P_0` | 2.400 | 2.375 |   -24.8 | high |
| `1^1P_1` | 2.440 | 2.475 |   +35.3 | high |
| `1^3P_1` | 2.490 | 2.501 |   +11.4 | high |
| `1^3P_2` | 2.500 | 2.480 |   -20.4 | high |
| `1^3D_3` | 2.830 | 2.783 |   -46.6 | high |
| `1^3F_4` | 3.110 | 3.040 |   -69.5 | medium |
| `1^1S_0` | 1.980 | 1.933 |   -47.0 | high |
| `2^1S_0` | 2.670 | 2.622 |   -48.3 | high |
| `1^3S_1` | 2.130 | 2.144 |   +14.0 | high |
| `2^3S_1` | 2.730 | 2.752 |   +21.7 | high |
| `1^3D_1` | 2.900 | 2.910 |   +10.4 | high |
| `1^3P_0` | 2.480 | 2.463 |   -16.9 | high |
| `1^1P_1` | 2.530 | 2.564 |   +34.3 | high |
| `1^3P_1` | 2.570 | 2.577 |    +6.6 | high |
| `1^3P_2` | 2.590 | 2.577 |   -12.8 | high |
| `1^3D_3` | 2.920 | 2.902 |   -17.7 | high |
| `1^3F_4` | 3.190 | 3.173 |   -17.0 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.001 |  -200.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -200.8 | 1.800 |
| `2^1S_0` | 2.629 |  -124.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -124.8 | 2.504 |
| `1^3S_1` | 2.001 |   +66.9 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +66.9 | 2.068 |
| `2^3S_1` | 2.629 |   +41.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +41.6 | 2.670 |
| `1^3D_1` | 2.831 |    +0.0 |  -176.1 |  +240.0 |   +63.9 |   -17.0 |   +46.9 | 2.878 |
| `1^3P_0` | 2.475 |    +0.0 |  -299.2 |  +275.3 |   -23.9 |   -76.2 |  -100.1 | 2.375 |
| `1^1P_1` | 2.475 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.475 |
| `1^3P_1` | 2.475 |    +0.0 |  -149.6 |  +137.6 |   -12.0 |   +38.1 |   +26.1 | 2.501 |
| `1^3P_2` | 2.475 |    +0.0 |  +149.6 |  -137.6 |   +12.0 |    -7.6 |    +4.3 | 2.480 |
| `1^3D_3` | 2.831 |    +0.0 |  +117.4 |  -160.0 |   -42.6 |    -4.8 |   -47.4 | 2.783 |
| `1^3F_4` | 3.124 |    +0.0 |   +88.4 |  -168.7 |   -80.4 |    -3.0 |   -83.4 | 3.040 |
| `1^1S_0` | 2.091 |  -158.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -158.2 | 1.933 |
| `2^1S_0` | 2.719 |   -97.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -97.5 | 2.622 |
| `1^3S_1` | 2.091 |   +52.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +52.7 | 2.144 |
| `2^3S_1` | 2.719 |   +32.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +32.5 | 2.752 |
| `1^3D_1` | 2.913 |    +0.0 |   -84.5 |   +95.0 |   +10.5 |   -13.1 |    -2.7 | 2.910 |
| `1^3P_0` | 2.564 |    +0.0 |  -150.9 |  +112.6 |   -38.3 |   -62.9 |  -101.3 | 2.463 |
| `1^1P_1` | 2.564 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.564 |
| `1^3P_1` | 2.564 |    +0.0 |   -75.5 |   +56.3 |   -19.2 |   +31.5 |   +12.3 | 2.577 |
| `1^3P_2` | 2.564 |    +0.0 |   +75.5 |   -56.3 |   +19.2 |    -6.3 |   +12.9 | 2.577 |
| `1^3D_3` | 2.913 |    +0.0 |   +56.3 |   -63.3 |    -7.0 |    -3.7 |   -10.7 | 2.902 |
| `1^3F_4` | 3.200 |    +0.0 |   +41.7 |   -66.2 |   -24.5 |    -2.3 |   -26.7 | 3.173 |

## Fine-Structure Mass Convention (audit note)

Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_1` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_0` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1P_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_1` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 1.628000 | 0.419000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 33.0 MeV.
Max absolute residual: 79.5 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.858 |   -12.3 |
| `1F` | 2 | 3.150 | 3.107 |   -43.2 |
| `1P` | 8 | 2.518 | 2.520 |    +2.3 |
| `1S` | 4 | 2.046 | 2.046 |    +0.0 |
| `2S` | 4 | 2.670 | 2.674 |    +4.0 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
