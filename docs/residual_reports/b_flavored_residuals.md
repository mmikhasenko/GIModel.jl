# Residuals: b_flavored (GI-style)

Model: finite-difference + `relativistic` kinetic, with smeared S-wave contact hyperfine, first-order L·S (vector+Thomas) and OGE-tensor; GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`; pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), see `src/GIModel/`.

| state | reference GeV | baseline GeV | residual MeV | confidence |
|---|---:|---:|---:|---|
| `1^1S_0` | 5.310 | 5.236 |   -73.8 | high |
| `2^1S_0` | 5.900 | 5.838 |   -61.5 | high |
| `1^3S_1` | 5.370 | 5.386 |   +15.7 | high |
| `2^3S_1` | 5.930 | 5.967 |   +36.6 | high |
| `1^3P_2` | 5.800 | 6.104 |  +303.5 | high |
| `1^3D_3` | 6.110 | 6.213 |  +102.8 | medium |
| `1^3F_4` | 6.360 | 6.351 |    -8.6 | medium |
| `1^1S_0` | 5.390 | 5.373 |   -17.3 | high |
| `2^1S_0` | 5.980 | 5.970 |   -10.2 | high |
| `1^3S_1` | 5.450 | 5.470 |   +20.3 | high |
| `2^3S_1` | 6.010 | 6.052 |   +41.9 | high |
| `1^3P_2` | 5.880 | 5.988 |  +108.0 | high |
| `1^3D_3` | 6.180 | 6.238 |   +58.0 | medium |
| `1^3F_4` | 6.430 | 6.459 |   +28.7 | medium |
| `1^1S_0` | 6.270 | 6.335 |   +64.8 | high |
| `2^1S_0` | 6.850 | 6.910 |   +59.6 | high |
| `1^3S_1` | 6.340 | 6.392 |   +52.2 | high |
| `2^3S_1` | 6.890 | 6.950 |   +59.9 | high |
| `1^3P_2` | 6.770 | 6.808 |   +37.8 | high |
| `1^3D_3` | 7.040 | 7.080 |   +39.7 | medium |
| `1^3F_4` | 7.270 | 7.302 |   +32.4 | medium |

## Contribution Breakdown (diagnostic)

All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).

| state | central GeV | contact MeV | L·S MeV | tensor MeV | total shift MeV | predicted GeV |
|---|---:|---:|---:|---:|---:|---:|
| `1^1S_0` | 5.348 |  -112.2 |    +0.0 |    +0.0 |  -112.2 | 5.236 |
| `2^1S_0` | 5.935 |   -96.1 |    +0.0 |    +0.0 |   -96.1 | 5.838 |
| `1^3S_1` | 5.348 |   +37.4 |    +0.0 |    +0.0 |   +37.4 | 5.386 |
| `2^3S_1` | 5.935 |   +32.0 |    +0.0 |    +0.0 |   +32.0 | 5.967 |
| `1^3P_2` | 5.811 |    +0.0 |  +296.2 |    -3.3 |  +292.9 | 6.104 |
| `1^3D_3` | 6.134 |    +0.0 |   +80.0 |    -1.5 |   +78.5 | 6.213 |
| `1^3F_4` | 6.399 |    +0.0 |   -46.8 |    -0.9 |   -47.7 | 6.351 |
| `1^1S_0` | 5.446 |   -73.2 |    +0.0 |    +0.0 |   -73.2 | 5.373 |
| `2^1S_0` | 6.031 |   -61.6 |    +0.0 |    +0.0 |   -61.6 | 5.970 |
| `1^3S_1` | 5.446 |   +24.4 |    +0.0 |    +0.0 |   +24.4 | 5.470 |
| `2^3S_1` | 6.031 |   +20.5 |    +0.0 |    +0.0 |   +20.5 | 6.052 |
| `1^3P_2` | 5.896 |    +0.0 |   +93.7 |    -1.9 |   +91.8 | 5.988 |
| `1^3D_3` | 6.211 |    +0.0 |   +27.9 |    -0.9 |   +27.0 | 6.238 |
| `1^3F_4` | 6.469 |    +0.0 |    -9.9 |    -0.5 |   -10.3 | 6.459 |
| `1^1S_0` | 6.378 |   -43.1 |    +0.0 |    +0.0 |   -43.1 | 6.335 |
| `2^1S_0` | 6.940 |   -30.2 |    +0.0 |    +0.0 |   -30.2 | 6.910 |
| `1^3S_1` | 6.378 |   +14.4 |    +0.0 |    +0.0 |   +14.4 | 6.392 |
| `2^3S_1` | 6.940 |   +10.1 |    +0.0 |    +0.0 |   +10.1 | 6.950 |
| `1^3P_2` | 6.795 |    +0.0 |   +13.7 |    -0.9 |   +12.8 | 6.808 |
| `1^3D_3` | 7.074 |    +0.0 |    +5.8 |    -0.4 |    +5.4 | 7.080 |
| `1^3F_4` | 7.301 |    +0.0 |    +1.5 |    -0.2 |    +1.3 | 7.302 |

Mean absolute residual: 58.7 MeV.
Max absolute residual: 303.5 MeV.

## Spin-Averaged Diagnostics

Weighted by `2J+1` within each available `(n, L)` group.

| multiplet | states | reference GeV | baseline GeV | residual MeV |
|---|---:|---:|---:|---:|
| `1D` | 3 | 6.443 | 6.510 |   +66.8 |
| `1F` | 3 | 6.687 | 6.704 |   +17.5 |
| `1P` | 3 | 6.150 | 6.300 |  +149.8 |
| `1S` | 6 | 5.704 | 5.724 |   +19.9 |
| `2S` | 6 | 6.268 | 6.302 |   +33.6 |

This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.
