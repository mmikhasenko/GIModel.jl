# Formula Map

Every implemented formula must be mapped to the original Godfrey-Isgur paper by
equation, paragraph, page, or table.

No model term should be introduced in code without an entry here.
# Formula Map

This file maps implemented local code to the Godfrey-Isgur paper. It is also a
guardrail: when the code is only a diagnostic approximation, say so here.

## Implemented Baseline

- `src/GIModel/GIModel.jl`: semirelativistic kinetic operator
  `sqrt(p^2 + m_1^2) + sqrt(p^2 + m_2^2)`.
  - Paper anchor: Eq. (1b), PDF page 2.
  - Numerical implementation: finite-difference radial `p^2` operator and
    dense matrix square root.

- `src/GIModel/GIModel.jl`: central spin-independent potential
  `b r - 4 alpha_s(r) / (3 r) + c`.
  - Paper anchor: nonrelativistic orientation around Eqs. (2)-(3), PDF pages
    2-3.
  - Parameters: `data/parameters.provisional.toml`, copied from Table II.

- `src/GIModel/GIModel.jl`: running Coulomb ansatz
  `alpha_s(r) = sum_k alpha_k erf(gamma_k r)`.
  - Paper anchor: Eq. (12), Eq. (13), and Fig. 2 caption, PDF page 3.
  - Coefficients used from Fig. 2 caption:
    `alpha_k = (0.25, 0.15, 0.20)` and momentum-space denominators
    `(1, 10, 1000) GeV^2`, giving
    `gamma_k = (1/2, sqrt(10)/2, sqrt(1000)/2) GeV`.

- `src/GIModel/GIModel.jl`: smeared S-wave contact hyperfine shift.
  - Paper anchor: color hyperfine term around Eq. (4), PDF page 2, and
    smearing discussion in Appendix A, PDF pages 36-37.
  - Current status: first diagnostic implementation. It uses the Table II
    `sigma0` and `s` smearing parameters with the standard GI mass-dependent
    Gaussian width form. The full momentum-dependent relativization factors
    are not yet included.

## Not Yet Implemented

- Full GI smeared Coulomb/confinement effective potentials from Appendix A.
- Momentum-dependent relativization factors for contact, tensor, vector
  spin-orbit, and scalar spin-orbit interactions.
- Tensor and spin-orbit fine structure.
- Unequal-mass antisymmetric spin-orbit mixing.
- Isoscalar annihilation mixing.
