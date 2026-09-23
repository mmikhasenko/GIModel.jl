# Phase III completion record

Phase III ends at a symbolic boundary. It derives every helicity and
partial-wave coefficient as a linear combination of named spatial integrals;
it does not assign numerical values to those integrals. Native Eq. (19)
integration and normalization conversion belong to Phase IV.

## One representation

`_AngularCoefficientDecomposition` is the sole assembled algebra result:

- columns are a shared tuple of spatial-integral labels;
- `helicity_coefficients` maps those columns to the complete compressed
  Appendix-C helicity vector;
- `partial_wave_coefficients` is exactly the Appendix-C projection matrix
  multiplied by `helicity_coefficients`.

Thus S/D or P/F waves never own or recompute separate integral collections.
The spin, flavor, and topology term records remain internal provenance used to
construct this matrix. No Phase-III type or function is exported.

## Evidence

The tests cover:

- all Table-XI signs and magnitudes for parent `J=0:5`;
- independent hand-reduced coefficient matrices for `rho -> pi pi` and
  `A1 -> rho pi`;
- the two S/D waves of `A1 -> rho pi` over one identical column basis;
- direct and recoil Eq. (19) pieces and quark/antiquark routing;
- an unequal-mass heavy-light flavor transition without a special angular
  implementation;
- exact spectator/OZI zeros;
- pure-`ss` `f' -> pi pi` zero and a nonzero physical mixing contribution;
- coherent mixed external states, overall-phase covariance, and interference;
- identical-daughter selection and its single normalization factor.

The legacy Table-IV/V backend remains frozen as the equal-beta numerical
oracle. Reducing the named spatial integrals to its `A,S,D,P,...` functions is
the analytic-SHO cross-check in Phase IV, because it depends on the spatial
operator rather than on Phase-III recoupling alone.

## Remaining Phase-IV inputs

- define `q'` from constituent masses for each emitting topology;
- implement the derivative and plane-wave spatial integrals;
- supply `g,h` and the supported complex-momentum domain;
- convert Appendix-C amplitudes to canonical relativistic normalization;
- compare analytic SHO, numerical SHO, native HO, and native FD evaluations.
