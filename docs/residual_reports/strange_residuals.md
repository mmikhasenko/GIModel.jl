# Residuals: strange (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A smearing), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 0.470 | 0.385 |   -84.6 | high |
| `2^1S_0` | 1.450 | 0.992 |  -457.5 | high |
| `3^1S_0` | 2.020 | 1.517 |  -503.2 | high |
| `1^3S_1` | 0.900 | 1.258 |  +357.7 | high |
| `2^3S_1` | 1.580 | 1.974 |  +394.2 | high |
| `1^3D_1` | 1.780 | 1.654 |  -126.3 | high |
| `3^3S_1` | 2.110 | 2.479 |  +368.9 | high |
| `2^3D_1` | 2.250 | 2.145 |  -104.8 | high |
| `1^3P_0` | 1.240 | 1.334 |   +94.2 | high |
| `2^3P_0` | 1.890 | 1.946 |   +55.9 | high |
| `1^1P_1` | 1.340 | 1.552 |  +212.0 | high |
| `1^3P_1` | 1.380 | 1.498 |  +117.6 | high |
| `2^1P_1` | 1.900 | 2.079 |  +179.0 | high |
| `2^3P_1` | 1.930 | 2.063 |  +133.5 | high |
| `1^3P_2` | 1.430 | 1.621 |  +190.9 | high |
| `2^3P_2` | 1.940 | 2.108 |  +168.1 | high |
| `1^3F_2` | 2.150 | 1.915 |  -235.1 | high |
| `1^1D_2` | 1.780 | 1.930 |  +150.3 | high |
| `1^3D_2` | 1.810 | 1.846 |   +35.8 | high |
| `2^1D_2` | 2.230 | 2.376 |  +145.7 | high |
| `2^3D_2` | 2.260 | 2.307 |   +47.4 | high |
| `1^3D_3` | 1.790 | 2.108 |  +317.5 | high |
| `2^3D_3` | 2.240 | 2.522 |  +281.6 | high |
| `1^3G_3` | 2.460 | 2.142 |  -318.3 | medium |
| `1^1F_3` | 2.120 | 2.246 |  +126.3 | high |
| `1^3F_3` | 2.150 | 2.162 |   +11.8 | high |
| `1^3F_4` | 2.110 | 2.497 |  +387.1 | high |
| `1^1G_4` | 2.410 | 2.524 |  +114.1 | medium |
| `1^3G_4` | 2.440 | 2.447 |    +7.1 | medium |
| `1^3G_5` | 2.390 | 2.831 |  +440.7 | medium |

Mean absolute residual: 205.6 MeV.
Max absolute residual: 503.2 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 4 | 1.791 | 1.930 |  +138.7 |
| `2D` | 4 | 2.244 | 2.375 |  +131.1 |
| `1F` | 4 | 2.130 | 2.247 |  +117.0 |
| `1G` | 4 | 2.421 | 2.524 |  +103.1 |
| `1P` | 4 | 1.379 | 1.549 |  +169.8 |
| `2P` | 4 | 1.923 | 2.076 |  +152.8 |
| `1S` | 2 | 0.792 | 1.040 |  +247.1 |
| `2S` | 2 | 1.548 | 1.729 |  +181.3 |
| `3S` | 2 | 2.087 | 2.238 |  +150.9 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
