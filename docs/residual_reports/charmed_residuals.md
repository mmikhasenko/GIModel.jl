# Residuals: charmed (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 1.880 | 1.843 |   -37.2 | high |
| `2^1S_0` | 2.580 | 2.501 |   -79.2 | high |
| `1^3S_1` | 2.040 | 2.165 |  +124.9 | high |
| `2^3S_1` | 2.640 | 2.789 |  +149.1 | high |
| `1^3D_1` | 2.820 | 2.848 |   +28.0 | high |
| `1^3P_0` | 2.400 | 2.008 |  -392.3 | high |
| `1^1P_1` | 2.440 | 2.561 |  +121.0 | high |
| `1^3P_1` | 2.490 | 2.476 |   -13.9 | high |
| `1^3P_2` | 2.500 | 2.723 |  +222.5 | high |
| `1^3D_3` | 2.830 | 2.908 |   +78.0 | high |
| `1^3F_4` | 3.110 | 3.093 |   -16.7 | medium |
| `1^1S_0` | 1.980 | 2.039 |   +58.9 | high |
| `2^1S_0` | 2.670 | 2.688 |   +17.5 | high |
| `1^3S_1` | 2.130 | 2.241 |  +111.1 | high |
| `2^3S_1` | 2.730 | 2.864 |  +134.4 | high |
| `1^3D_1` | 2.900 | 2.957 |   +56.8 | high |
| `1^3P_0` | 2.480 | 2.425 |   -55.2 | high |
| `1^1P_1` | 2.530 | 2.654 |  +124.3 | high |
| `1^3P_1` | 2.570 | 2.650 |   +79.7 | high |
| `1^3P_2` | 2.590 | 2.703 |  +113.0 | high |
| `1^3D_3` | 2.920 | 2.990 |   +70.2 | high |
| `1^3F_4` | 3.190 | 3.239 |   +48.6 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 2.084 |  -241.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -241.6 | 1.843 |
| `2^1S_0` | 2.717 |  -216.3 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -216.3 | 2.501 |
| `1^3S_1` | 2.084 |   +80.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +80.5 | 2.165 |
| `2^3S_1` | 2.717 |   +72.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +72.1 | 2.789 |
| `1^3D_1` | 2.904 |    +0.0 |  -361.5 |  +340.8 |   -20.7 |   -35.7 |   -56.4 | 2.848 |
| `1^3P_0` | 2.561 |    +0.0 |  -671.1 |  +309.6 |  -361.5 |  -191.8 |  -553.3 | 2.008 |
| `1^1P_1` | 2.561 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.561 |
| `1^3P_1` | 2.561 |    +0.0 |  -335.6 |  +154.8 |  -180.7 |   +95.9 |   -84.8 | 2.476 |
| `1^3P_2` | 2.561 |    +0.0 |  +335.6 |  -154.8 |  +180.7 |   -19.2 |  +161.6 | 2.723 |
| `1^3D_3` | 2.904 |    +0.0 |  +241.0 |  -227.2 |   +13.8 |   -10.2 |    +3.6 | 2.908 |
| `1^3F_4` | 3.189 |    +0.0 |  +188.6 |  -277.9 |   -89.3 |    -6.5 |   -95.8 | 3.093 |
| `1^1S_0` | 2.191 |  -151.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -151.6 | 2.039 |
| `2^1S_0` | 2.820 |  -132.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -132.6 | 2.688 |
| `1^3S_1` | 2.191 |   +50.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +50.5 | 2.241 |
| `2^3S_1` | 2.820 |   +44.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +44.2 | 2.864 |
| `1^3D_1` | 2.988 |    +0.0 |  -112.8 |  +101.4 |   -11.4 |   -20.2 |   -31.6 | 2.957 |
| `1^3P_0` | 2.654 |    +0.0 |  -212.2 |   +92.8 |  -119.4 |  -110.1 |  -229.5 | 2.425 |
| `1^1P_1` | 2.654 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 | 2.654 |
| `1^3P_1` | 2.654 |    +0.0 |  -106.1 |   +46.4 |   -59.7 |   +55.1 |    -4.6 | 2.650 |
| `1^3P_2` | 2.654 |    +0.0 |  +106.1 |   -46.4 |   +59.7 |   -11.0 |   +48.7 | 2.703 |
| `1^3D_3` | 2.988 |    +0.0 |   +75.2 |   -67.6 |    +7.6 |    -5.8 |    +1.8 | 2.990 |
| `1^3F_4` | 3.266 |    +0.0 |   +58.3 |   -82.3 |   -24.0 |    -3.6 |   -27.6 | 3.239 |

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

Mean absolute residual: 96.9 MeV.
Max absolute residual: 392.3 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 2.870 | 2.935 |   +64.6 |
| `1F` | 2 | 3.150 | 3.166 |   +15.9 |
| `1P` | 8 | 2.518 | 2.608 |   +90.2 |
| `1S` | 4 | 2.046 | 2.137 |   +91.2 |
| `2S` | 4 | 2.670 | 2.769 |   +98.6 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
