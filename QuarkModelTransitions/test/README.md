# Numerical test coverage

Run from the repository root:

```sh
julia --project=QuarkModelTransitions QuarkModelTransitions/test/runtests.jl
```

Tests protect calculated amplitudes and widths, their normalization and units,
physical selection rules, coherent phases, and numerical accuracy. They do not
freeze export lists, constructor layouts, exception types, rejected keywords,
printed output, or compatibility with older calling conventions.

| Calculation | Independent check |
| --- | --- |
| Flavor and spin | Appendix-B signs, forbidden matrix elements, tensor Hermiticity |
| Angular decomposition | Hand-reduced Eq. (19) coefficients and Appendix-C projection matrices |
| Strong emission | Analytic Gaussian integrals, S/D ratio, interference, identical-particle normalization |
| Photon emission | M1/E1/M2 kernel values, recoil, angular factors, momentum powers, closed-channel widths |
| Annihilation | Eq. (17) gluonic prefactors, electromagnetic charge coefficients, coherent sums, physical width conversions |
| Mass corrections | Analytic mass/momentum powers and recomputation at changed external masses with fixed waves |
| Wave integration | Sampled-mesh versus analytic oscillator integrals, independently solved HO/FD waves, mesh refinement |
| Reference formulas | Table-IV/V reduced amplitudes and selected Table-VII observable benchmarks |

Use exact comparisons for exact selection-rule zeros. Algebraic identities use
floating-point tolerances; analytic integral comparisons specify tolerances near
1e-10, while sampled-mesh integrals allow discretization errors near 2e-5.
Solved strong-emission HO/FD overlaps use 0.1–0.3% tolerances, with the looser
bound reserved for node-sensitive S-wave overlaps. The smeared-origin comparison
uses 1%. These are numerical agreements, not estimates of model accuracy.

The wider windows around historical paper values are model benchmarks with
stated wave and mass prescriptions. They are not high-precision numerical
reference values. Do not tighten them by recording the current implementation's
output as an independent expected value.

Keep solver eigenvalue and mixing-diagonalization tests in GIModel. This suite
checks the transition integrals that consume those waves. Keep paper data loading
and report validation in GIPaper. The documentation audit remains a separate
final check: it executes copyable examples and validates help links dynamically,
without a frozen list of public names.

For a new regression, prefer one case that distinguishes the incorrect numerical
result from the correct one. Reuse existing fixtures and remove a superseded
check instead of accumulating another layer of historical API tests.
