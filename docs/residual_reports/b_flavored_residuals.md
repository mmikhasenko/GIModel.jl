# Residuals: b_flavored (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.236 |   -73.8 | high |
| `2^1S_0` | 5.900 | 5.838 |   -61.5 | high |
| `1^3S_1` | 5.370 | 5.386 |   +15.7 | high |
| `2^3S_1` | 5.930 | 5.967 |   +36.6 | high |
| `1^3P_2` | 5.800 | 6.094 |  +293.6 | high |
| `1^3D_3` | 6.110 | 6.208 |   +98.4 | medium |
| `1^3F_4` | 6.360 | 6.349 |   -11.1 | medium |
| `1^1S_0` | 5.390 | 5.373 |   -17.3 | high |
| `2^1S_0` | 5.980 | 5.970 |   -10.2 | high |
| `1^3S_1` | 5.450 | 5.470 |   +20.3 | high |
| `2^3S_1` | 6.010 | 6.052 |   +41.9 | high |
| `1^3P_2` | 5.880 | 5.982 |  +102.2 | high |
| `1^3D_3` | 6.180 | 6.235 |   +55.4 | medium |
| `1^3F_4` | 6.430 | 6.457 |   +27.2 | medium |
| `1^1S_0` | 6.270 | 6.335 |   +64.8 | high |
| `2^1S_0` | 6.850 | 6.910 |   +59.6 | high |
| `1^3S_1` | 6.340 | 6.392 |   +52.2 | high |
| `2^3S_1` | 6.890 | 6.950 |   +59.9 | high |
| `1^3P_2` | 6.770 | 6.805 |   +35.1 | high |
| `1^3D_3` | 7.040 | 7.079 |   +38.5 | medium |
| `1^3F_4` | 7.270 | 7.302 |   +31.8 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 5.348 |  -112.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |  -112.2 | 5.236 |
| `2^1S_0` | 5.935 |   -96.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -96.1 | 5.838 |
| `1^3S_1` | 5.348 |   +37.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +37.4 | 5.386 |
| `2^3S_1` | 5.935 |   +32.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +32.0 | 5.967 |
| `1^3P_2` | 5.811 |    +0.0 |  +469.4 |  -173.3 |  +296.2 |   -13.2 |  +283.0 | 6.094 |
| `1^3D_3` | 6.134 |    +0.0 |  +332.8 |  -252.8 |   +80.0 |    -6.0 |   +74.0 | 6.208 |
| `1^3F_4` | 6.399 |    +0.0 |  +261.9 |  -308.7 |   -46.8 |    -3.4 |   -50.2 | 6.349 |
| `1^1S_0` | 5.446 |   -73.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -73.2 | 5.373 |
| `2^1S_0` | 6.031 |   -61.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -61.6 | 5.970 |
| `1^3S_1` | 5.446 |   +24.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +24.4 | 5.470 |
| `2^3S_1` | 6.031 |   +20.5 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +20.5 | 6.052 |
| `1^3P_2` | 5.896 |    +0.0 |  +143.8 |   -50.1 |   +93.7 |    -7.7 |   +86.0 | 5.982 |
| `1^3D_3` | 6.211 |    +0.0 |  +100.3 |   -72.4 |   +27.9 |    -3.4 |   +24.5 | 6.235 |
| `1^3F_4` | 6.469 |    +0.0 |   +78.1 |   -88.0 |    -9.9 |    -1.9 |   -11.8 | 6.457 |
| `1^1S_0` | 6.378 |   -43.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -43.1 | 6.335 |
| `2^1S_0` | 6.940 |   -30.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -30.2 | 6.910 |
| `1^3S_1` | 6.378 |   +14.4 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +14.4 | 6.392 |
| `2^3S_1` | 6.940 |   +10.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +10.1 | 6.950 |
| `1^3P_2` | 6.795 |    +0.0 |   +18.3 |    -4.6 |   +13.7 |    -3.6 |   +10.1 | 6.805 |
| `1^3D_3` | 7.074 |    +0.0 |   +12.2 |    -6.5 |    +5.8 |    -1.5 |    +4.2 | 7.079 |
| `1^3F_4` | 7.301 |    +0.0 |    +9.2 |    -7.7 |    +1.5 |    -0.8 |    +0.7 | 7.302 |

## Fine-Structure Mass Convention (audit note)

Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `2^1S_0` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3S_1` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `2^3S_1` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3P_2` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 57.5 MeV.
Max absolute residual: 293.6 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.507 |   +64.1 |
| `1F` | 3 | 6.687 | 6.703 |   +16.0 |
| `1P` | 3 | 6.150 | 6.294 |  +143.6 |
| `1S` | 6 | 5.704 | 5.724 |   +19.9 |
| `2S` | 6 | 6.268 | 6.302 |   +33.6 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
