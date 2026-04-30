# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.385 |   -84.6 | high |
| `2^1S_0` | 1.450 | 0.992 |  -457.5 | high |
| `3^1S_0` | 2.020 | 1.517 |  -503.2 | high |
| `1^3S_1` | 0.900 | 1.258 |  +357.7 | high |
| `2^3S_1` | 1.580 | 1.974 |  +394.2 | high |
| `1^3D_1` | 1.780 | 1.905 |  +125.1 | high |
| `3^3S_1` | 2.110 | 2.479 |  +368.9 | high |
| `2^3D_1` | 2.250 | 2.198 |   -52.0 | high |
| `1^3P_0` | 1.240 | 0.748 |  -492.1 | high |
| `2^3P_0` | 1.890 | 0.887 | -1002.7 | high |
| `1^1P_1` | 1.340 | 1.552 |  +212.0 | high |
| `1^3P_1` | 1.380 | 1.651 |  +271.1 | high |
| `2^1P_1` | 1.900 | 2.079 |  +179.0 | high |
| `2^3P_1` | 1.930 | 2.104 |  +174.5 | high |
| `1^3P_2` | 1.430 | 1.653 |  +223.4 | high |
| `2^3P_2` | 1.940 | 2.302 |  +362.0 | high |
| `1^3F_2` | 2.150 | 2.414 |  +263.9 | high |
| `1^1D_2` | 1.780 | 1.930 |  +150.3 | high |
| `1^3D_2` | 1.810 | 2.049 |  +239.4 | high |
| `2^1D_2` | 2.230 | 2.376 |  +145.7 | high |
| `2^3D_2` | 2.260 | 2.472 |  +211.9 | high |
| `1^3D_3` | 1.790 | 1.856 |   +66.0 | high |
| `2^3D_3` | 2.240 | 2.383 |  +143.2 | high |
| `1^3G_3` | 2.460 | 2.791 |  +331.0 | medium |
| `1^1F_3` | 2.120 | 2.246 |  +126.3 | high |
| `1^3F_3` | 2.150 | 2.351 |  +201.1 | high |
| `1^3F_4` | 2.110 | 2.072 |   -38.3 | high |
| `1^1G_4` | 2.410 | 2.524 |  +114.1 | medium |
| `1^3G_4` | 2.440 | 2.616 |  +175.7 | medium |
| `1^3G_5` | 2.390 | 2.279 |  -110.7 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 1.040 |  -654.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -654.2 | 0.385 |
| `2^1S_0` | 1.729 |  -736.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -736.3 | 0.992 |
| `3^1S_0` | 2.238 |  -721.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -721.6 | 1.517 |
| `1^3S_1` | 1.040 |  +218.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +218.1 | 1.258 |
| `2^3S_1` | 1.729 |  +245.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +245.4 | 1.974 |
| `1^3D_1` | 1.930 |    +0.0 |  -404.9 |  +475.4 |   +70.5 |   -95.6 |   -25.1 | 1.905 |
| `3^3S_1` | 2.238 |  +240.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +240.5 | 2.479 |
| `2^3D_1` | 2.376 |    +0.0 |  -545.3 |  +484.2 |   -61.2 |  -116.6 |  -177.7 | 2.198 |
| `1^3P_0` | 1.552 |    +0.0 |  -828.2 |  +525.2 |  -302.9 |  -501.2 |  -804.1 | 0.748 |
| `2^3P_0` | 2.079 |    +0.0 | -1113.8 |  +543.5 |  -570.4 |  -621.4 | -1191.7 | 0.887 |
| `1^1P_1` | 1.552 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.552 |
| `1^3P_1` | 1.552 |    +0.0 |  -414.1 |  +262.6 |  -151.5 |  +250.6 |   +99.1 | 1.651 |
| `2^1P_1` | 2.079 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.079 |
| `2^3P_1` | 2.079 |    +0.0 |  -556.9 |  +271.7 |  -285.2 |  +310.7 |   +25.5 | 2.104 |
| `1^3P_2` | 1.552 |    +0.0 |  +414.1 |  -262.6 |  +151.5 |   -50.1 |  +101.4 | 1.653 |
| `2^3P_2` | 2.079 |    +0.0 |  +556.9 |  -271.7 |  +285.2 |   -62.1 |  +223.0 | 2.302 |
| `1^3F_2` | 2.246 |    +0.0 |  -270.2 |  +479.7 |  +209.5 |   -41.9 |  +167.6 | 2.414 |
| `1^1D_2` | 1.930 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.930 |
| `1^3D_2` | 1.930 |    +0.0 |  -135.0 |  +158.5 |   +23.5 |   +95.6 |  +119.1 | 2.049 |
| `2^1D_2` | 2.376 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.376 |
| `2^3D_2` | 2.376 |    +0.0 |  -181.8 |  +161.4 |   -20.4 |  +116.6 |   +96.2 | 2.472 |
| `1^3D_3` | 1.930 |    +0.0 |  +269.9 |  -316.9 |   -47.0 |   -27.3 |   -74.3 | 1.856 |
| `2^3D_3` | 2.376 |    +0.0 |  +363.5 |  -322.8 |   +40.8 |   -33.3 |    +7.5 | 2.383 |
| `1^3G_3` | 2.524 |    +0.0 |  -205.9 |  +496.7 |  +290.8 |   -23.9 |  +267.0 | 2.791 |
| `1^1F_3` | 2.246 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.246 |
| `1^3F_3` | 2.246 |    +0.0 |   -67.5 |  +119.9 |   +52.4 |   +52.4 |  +104.8 | 2.351 |
| `1^3F_4` | 2.246 |    +0.0 |  +202.6 |  -359.8 |  -157.1 |   -17.5 |  -174.6 | 2.072 |
| `1^1G_4` | 2.524 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.524 |
| `1^3G_4` | 2.524 |    +0.0 |   -41.2 |   +99.3 |   +58.2 |   +33.4 |   +91.6 | 2.616 |
| `1^3G_5` | 2.524 |    +0.0 |  +164.7 |  -397.4 |  -232.7 |   -12.2 |  -244.8 | 2.279 |

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

Mean absolute residual: 252.6 MeV.
Max absolute residual: 1002.7 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.791 | 1.930 |  +139.3 |
| `2D` | 4 | 2.244 | 2.376 |  +131.7 |
| `1F` | 4 | 2.130 | 2.246 |  +116.7 |
| `1G` | 4 | 2.421 | 2.524 |  +103.0 |
| `1P` | 4 | 1.379 | 1.552 |  +172.8 |
| `2P` | 4 | 1.923 | 2.079 |  +155.6 |
| `1S` | 2 | 0.792 | 1.040 |  +247.1 |
| `2S` | 2 | 1.548 | 1.729 |  +181.3 |
| `3S` | 2 | 2.087 | 2.238 |  +150.9 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
