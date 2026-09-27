# Table VI independent solver comparison

Compares the complete native-HO run with the independent FD run at
ngrid=700, rmax=24 GeV^-1. This is an observable cross-check at the stated
resolution, not a new FD precision certification. Both runs use identical
canonical states, charges, kinematic inputs, and shared phase conventions.

| Multipole | Rows | Median relative difference | Largest absolute difference |
|---|---:|---:|---:|
| M1 | 42 | 0.099% | 0.00880238 |
| E1 | 35 | 0.047% | 0.00348363 |
| M2 | 2 | 0.094% | 0.000552937 |

Absolute differences are in μN for M1 and MeV^(1/2) for E1/M2.
The following rows expose radial mixing and cancellation sensitivity.

| Parent -> daughter | HO | FD | Paper |
|---|---:|---:|---:|
| eta_r -> omega | +0.0875515 | +0.0900856 | -0.18 |
| eta_r -> rho | +0.297925 | +0.305952 | +0.57 |
| eta_r -> phi | +0.502568 | +0.507555 | +0.37 |
| etaprime_r -> omega | -0.0722687 | -0.070018 | +0.01 |
| etaprime_r -> rho | -0.218603 | -0.212131 | +0.008 |
| etaprime_r -> phi | -0.026182 | -0.0300775 | -0.29 |
| Upsilondoubleprime -> etaprime_b | +0.00713656 | +0.00716937 | +0.007 |
| Upsilondoubleprime -> eta_b | -0.00373863 | -0.00375418 | -0.004 |
| Upsilon -> etaprime | -3.27426e-05 | -3.27729e-05 | +3e-06 |
| Upsilondoubleprime -> chi_0b | +0.00197291 | +0.00191614 | -0.002 |

The excited eta magnitude residuals survive the change of numerical
representation. They are not resolved by replacing the old central P
waves or fixing the duplicated perfect-mixing coefficient. Near-zero
amplitudes are shown in absolute units so a large relative error does not
conceal the cancellation scale. Canonical paper values are never replaced
by either solver's answer.
