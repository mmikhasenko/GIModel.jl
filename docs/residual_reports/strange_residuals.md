# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.385 |   -84.6 | high |
| `2^1S_0` | 1.450 | 0.992 |  -457.5 | high |
| `3^1S_0` | 2.020 | 1.517 |  -503.2 | high |
| `1^3S_1` | 0.900 | 1.258 |  +357.7 | high |
| `2^3S_1` | 1.580 | 1.974 |  +394.2 | high |
| `1^3D_1` | 1.780 | 1.965 |  +184.8 | high |
| `3^3S_1` | 2.110 | 2.479 |  +368.9 | high |
| `2^3D_1` | 2.250 | 2.299 |   +48.8 | high |
| `1^3P_0` | 1.240 | 1.150 |   -90.0 | high |
| `2^3P_0` | 1.890 | 1.440 |  -449.6 | high |
| `1^1P_1` | 1.340 | 1.552 |  +212.0 | high |
| `1^3P_1` | 1.380 | 1.525 |  +144.8 | high |
| `2^1P_1` | 1.900 | 2.079 |  +179.0 | high |
| `2^3P_1` | 1.930 | 1.993 |   +63.4 | high |
| `1^3P_2` | 1.430 | 1.649 |  +218.8 | high |
| `2^3P_2` | 1.940 | 2.258 |  +318.0 | high |
| `1^3F_2` | 2.150 | 2.425 |  +275.1 | high |
| `1^1D_2` | 1.780 | 1.930 |  +150.3 | high |
| `1^3D_2` | 1.810 | 1.980 |  +169.6 | high |
| `2^1D_2` | 2.230 | 2.376 |  +145.7 | high |
| `2^3D_2` | 2.260 | 2.401 |  +140.9 | high |
| `1^3D_3` | 1.790 | 1.880 |   +90.3 | high |
| `2^3D_3` | 2.240 | 2.391 |  +150.7 | high |
| `1^3G_3` | 2.460 | 2.789 |  +329.4 | medium |
| `1^1F_3` | 2.120 | 2.246 |  +126.3 | high |
| `1^3F_3` | 2.150 | 2.308 |  +158.0 | high |
| `1^3F_4` | 2.110 | 2.099 |   -11.0 | high |
| `1^1G_4` | 2.410 | 2.524 |  +114.1 | medium |
| `1^3G_4` | 2.440 | 2.587 |  +147.0 | medium |
| `1^3G_5` | 2.390 | 2.304 |   -86.3 | medium |

Mean absolute residual: 205.7 MeV.
Max absolute residual: 503.2 MeV.

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
