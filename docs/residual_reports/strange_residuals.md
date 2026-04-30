# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.385 |   -84.6 | high |
| `2^1S_0` | 1.450 | 0.992 |  -457.5 | high |
| `3^1S_0` | 2.020 | 1.517 |  -503.2 | high |
| `1^3S_1` | 0.900 | 1.258 |  +357.7 | high |
| `2^3S_1` | 1.580 | 1.974 |  +394.2 | high |
| `1^3D_1` | 1.780 | 1.880 |   +99.9 | high |
| `3^3S_1` | 2.110 | 2.479 |  +368.9 | high |
| `2^3D_1` | 2.250 | 2.184 |   -65.6 | high |
| `1^3P_0` | 1.240 | 0.629 |  -611.3 | high |
| `2^3P_0` | 1.890 | 0.739 | -1150.7 | high |
| `1^1P_1` | 1.340 | 1.552 |  +212.0 | high |
| `1^3P_1` | 1.380 | 1.785 |  +405.4 | high |
| `2^1P_1` | 1.900 | 2.079 |  +179.0 | high |
| `2^3P_1` | 1.930 | 2.344 |  +413.9 | high |
| `1^3P_2` | 1.430 | 1.597 |  +166.6 | high |
| `2^3P_2` | 1.940 | 2.188 |  +247.9 | high |
| `1^3F_2` | 2.150 | 2.391 |  +241.1 | high |
| `1^1D_2` | 1.780 | 1.930 |  +150.3 | high |
| `1^3D_2` | 1.810 | 2.065 |  +254.5 | high |
| `2^1D_2` | 2.230 | 2.376 |  +145.7 | high |
| `2^3D_2` | 2.260 | 2.515 |  +255.4 | high |
| `1^3D_3` | 1.790 | 1.856 |   +66.0 | high |
| `2^3D_3` | 2.240 | 2.358 |  +118.0 | high |
| `1^3G_3` | 2.460 | 2.771 |  +310.9 | medium |
| `1^1F_3` | 2.120 | 2.246 |  +126.3 | high |
| `1^3F_3` | 2.150 | 2.351 |  +200.5 | high |
| `1^3F_4` | 2.110 | 2.085 |   -25.2 | high |
| `1^1G_4` | 2.410 | 2.524 |  +114.1 | medium |
| `1^3G_4` | 2.440 | 2.613 |  +172.9 | medium |
| `1^3G_5` | 2.390 | 2.294 |   -95.7 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 1.040 |  -654.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -654.2 | 0.385 |
| `2^1S_0` | 1.729 |  -736.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -736.3 | 0.992 |
| `3^1S_0` | 2.238 |  -721.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -721.6 | 1.517 |
| `1^3S_1` | 1.040 |  +218.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +218.1 | 1.258 |
| `2^3S_1` | 1.729 |  +245.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +245.4 | 1.974 |
| `1^3D_1` | 1.930 |    +0.0 |  -306.6 |  +369.5 |   +62.9 |  -113.3 |   -50.4 | 1.880 |
| `3^3S_1` | 2.238 |  +240.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  +240.5 | 2.479 |
| `2^3D_1` | 2.376 |    +0.0 |  -388.7 |  +349.9 |   -38.7 |  -152.6 |  -191.3 | 2.184 |
| `1^3P_0` | 1.552 |    +0.0 |  -560.0 |  +331.8 |  -228.3 |  -695.1 |  -923.3 | 0.629 |
| `2^3P_0` | 2.079 |    +0.0 |  -704.8 |  +300.0 |  -404.9 |  -934.8 | -1339.6 | 0.739 |
| `1^1P_1` | 1.552 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.552 |
| `1^3P_1` | 1.552 |    +0.0 |  -280.0 |  +165.9 |  -114.1 |  +347.5 |  +233.4 | 1.785 |
| `2^1P_1` | 2.079 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.079 |
| `2^3P_1` | 2.079 |    +0.0 |  -352.4 |  +150.0 |  -202.4 |  +467.4 |  +265.0 | 2.344 |
| `1^3P_2` | 1.552 |    +0.0 |  +280.0 |  -165.9 |  +114.1 |   -69.5 |   +44.6 | 1.597 |
| `2^3P_2` | 2.079 |    +0.0 |  +352.4 |  -150.0 |  +202.4 |   -93.5 |  +109.0 | 2.188 |
| `1^3F_2` | 2.246 |    +0.0 |  -215.2 |  +405.3 |  +190.1 |   -45.3 |  +144.8 | 2.391 |
| `1^1D_2` | 1.930 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 1.930 |
| `1^3D_2` | 1.930 |    +0.0 |  -102.2 |  +123.2 |   +21.0 |  +113.3 |  +134.2 | 2.065 |
| `2^1D_2` | 2.376 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.376 |
| `2^3D_2` | 2.376 |    +0.0 |  -129.6 |  +116.6 |   -12.9 |  +152.6 |  +139.6 | 2.515 |
| `1^3D_3` | 1.930 |    +0.0 |  +204.4 |  -246.3 |   -41.9 |   -32.4 |   -74.3 | 1.856 |
| `2^3D_3` | 2.376 |    +0.0 |  +259.1 |  -233.3 |   +25.8 |   -43.6 |   -17.8 | 2.358 |
| `1^3G_3` | 2.524 |    +0.0 |  -167.4 |  +438.9 |  +271.5 |   -24.7 |  +246.8 | 2.771 |
| `1^1F_3` | 2.246 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.246 |
| `1^3F_3` | 2.246 |    +0.0 |   -53.8 |  +101.3 |   +47.5 |   +56.7 |  +104.2 | 2.351 |
| `1^3F_4` | 2.246 |    +0.0 |  +161.4 |  -304.0 |  -142.6 |   -18.9 |  -161.5 | 2.085 |
| `1^1G_4` | 2.524 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.524 |
| `1^3G_4` | 2.524 |    +0.0 |   -33.5 |   +87.8 |   +54.3 |   +34.6 |   +88.9 | 2.613 |
| `1^3G_5` | 2.524 |    +0.0 |  +133.9 |  -351.1 |  -217.2 |   -12.6 |  -229.8 | 2.294 |

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

Mean absolute residual: 266.5 MeV.
Max absolute residual: 1150.7 MeV.

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
