# Residuals: b_flavored (GI-style)

Model: with GI momentum-sandwiched smeared S-wave contact hyperfine, fixed-sector L·S (vector+Thomas) and OGE-tensor diagonalization with GI momentum-factor sandwiches and smeared-G/S derivative kernels; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` in the solver's native p² representation, see the **GIModel** package (`src/` at the repository root) for the computation and **GIPaper** (`GIPaper/src/`) for the comparison layer.

Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.310 |    -0.2 | high |
| `2^1S_0` | 5.900 | 5.904 |    +4.4 | high |
| `1^3S_1` | 5.370 | 5.369 |    -0.7 | high |
| `2^3S_1` | 5.930 | 5.933 |    +3.5 | high |
| `1^3P_2` | 5.800 | 5.795 |    -5.2 | high |
| `1^3D_3` | 6.110 | 6.105 |    -5.2 | medium |
| `1^3F_4` | 6.360 | 6.364 |    +3.9 | medium |
| `1^1S_0` | 5.390 | 5.391 |    +0.7 | high |
| `2^1S_0` | 5.980 | 5.985 |    +4.6 | high |
| `1^3S_1` | 5.450 | 5.447 |    -2.6 | high |
| `2^3S_1` | 6.010 | 6.013 |    +3.3 | high |
| `1^3P_2` | 5.880 | 5.874 |    -6.2 | high |
| `1^3D_3` | 6.180 | 6.178 |    -2.0 | medium |
| `1^3F_4` | 6.430 | 6.431 |    +1.1 | medium |
| `1^1S_0` | 6.270 | 6.265 |    -5.2 | high |
| `2^1S_0` | 6.850 | 6.855 |    +5.2 | high |
| `1^3S_1` | 6.340 | 6.333 |    -7.1 | high |
| `2^3S_1` | 6.890 | 6.889 |    -1.4 | high |
| `1^3P_2` | 6.770 | 6.765 |    -4.9 | high |
| `1^3D_3` | 7.040 | 7.044 |    +4.4 | medium |
| `1^3F_4` | 7.270 | 7.271 |    +0.5 | medium |

## Contribution Breakdown (diagnostic)

The central, contact, spin-orbit, and tensor columns are operator contributions in the same fully diagonalized fixed-sector eigenstate. They are not differences from a separate central-only solve.

| state | central contribution GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | annihilation MeV | spin + annihilation MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 5.358 |   -48.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -48.0 | 5.310 |
| `2^1S_0` | 5.926 |   -21.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -21.7 | 5.904 |
| `1^3S_1` | 5.356 |   +13.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.8 | 5.369 |
| `2^3S_1` | 5.926 |    +7.2 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +7.2 | 5.933 |
| `1^3P_2` | 5.788 |    +0.0 |   +34.6 |   -26.5 |    +8.0 |    -0.8 |    +0.0 |    +7.2 | 5.795 |
| `1^3D_3` | 6.111 |    +0.0 |   +25.5 |   -30.7 |    -5.2 |    -0.5 |    +0.0 |    -5.7 | 6.105 |
| `1^3F_4` | 6.376 |    +0.0 |   +18.7 |   -30.6 |   -11.9 |    -0.3 |    +0.0 |   -12.2 | 6.364 |
| `1^1S_0` | 5.437 |   -46.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -46.0 | 5.391 |
| `2^1S_0` | 6.006 |   -21.7 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -21.7 | 5.985 |
| `1^3S_1` | 5.434 |   +13.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +13.0 | 5.447 |
| `2^3S_1` | 6.006 |    +7.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +7.1 | 6.013 |
| `1^3P_2` | 5.865 |    +0.0 |   +30.4 |   -20.9 |    +9.5 |    -0.9 |    +0.0 |    +8.6 | 5.874 |
| `1^3D_3` | 6.182 |    +0.0 |   +22.6 |   -26.2 |    -3.7 |    -0.5 |    +0.0 |    -4.2 | 6.178 |
| `1^3F_4` | 6.442 |    +0.0 |   +16.7 |   -27.7 |   -11.0 |    -0.3 |    +0.0 |   -11.3 | 6.431 |
| `1^1S_0` | 6.324 |   -59.1 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -59.1 | 6.265 |
| `2^1S_0` | 6.882 |   -26.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   -26.6 | 6.855 |
| `1^3S_1` | 6.318 |   +14.6 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |   +14.6 | 6.333 |
| `2^3S_1` | 6.881 |    +7.8 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +0.0 |    +7.8 | 6.889 |
| `1^3P_2` | 6.752 |    +0.0 |   +21.1 |    -6.6 |   +14.5 |    -1.0 |    +0.0 |   +13.5 | 6.765 |
| `1^3D_3` | 7.039 |    +0.0 |   +14.7 |    -9.2 |    +5.6 |    -0.5 |    +0.0 |    +5.1 | 7.044 |
| `1^3F_4` | 7.270 |    +0.0 |   +11.0 |   -10.6 |    +0.4 |    -0.3 |    +0.0 |    +0.1 | 7.271 |

## Fine-Structure Mass Convention (audit note)

Unequal-mass diagonal fine structure still uses the symmetric `L·S` contraction. Same-`J` `^1L_J`/`^3L_J` rows are then corrected by the antisymmetric spin-orbit block when both partner rows are present; other unequal-mass rows remain under the symmetric convention.

| state | m1 GeV | m2 GeV | convention |
|---|---:|---:|---|
| `1^1S_0` | 4.977000 | 0.220000 | `disabled` |
| `2^1S_0` | 4.977000 | 0.220000 | `disabled` |
| `1^3S_1` | 4.977000 | 0.220000 | `disabled` |
| `2^3S_1` | 4.977000 | 0.220000 | `disabled` |
| `1^3P_2` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.220000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 0.419000 | `disabled` |
| `2^1S_0` | 4.977000 | 0.419000 | `disabled` |
| `1^3S_1` | 4.977000 | 0.419000 | `disabled` |
| `2^3S_1` | 4.977000 | 0.419000 | `disabled` |
| `1^3P_2` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 0.419000 | `unequal_mass_equal_share_LdotS` |
| `1^1S_0` | 4.977000 | 1.628000 | `disabled` |
| `2^1S_0` | 4.977000 | 1.628000 | `disabled` |
| `1^3S_1` | 4.977000 | 1.628000 | `disabled` |
| `2^3S_1` | 4.977000 | 1.628000 | `disabled` |
| `1^3P_2` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3D_3` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |
| `1^3F_4` | 4.977000 | 1.628000 | `unequal_mass_equal_share_LdotS` |

Mean absolute residual: 3.4 MeV.
Max absolute residual: 7.1 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.442 |    -0.9 |
| `1F` | 3 | 6.687 | 6.689 |    +1.9 |
| `1P` | 3 | 6.150 | 6.145 |    -5.4 |
| `1S` | 6 | 5.704 | 5.701 |    -3.0 |
| `2S` | 6 | 6.268 | 6.271 |    +2.5 |

The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.
