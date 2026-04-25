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
  - Parameters: `data/parameters.provisional.toml`, copied from Table II
    (Table II $\epsilon_i$ from (A10) are loaded as `relativistic_factors` but
    not yet used in the baseline Hamiltonian beyond parameter storage).

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

- `src/GIModel/GIModel.jl`: **experimental** 3D isotropic Gaussian smearing
  of pointwise $G(r)$ and $S(r)$ as in the structure of (A7)–(A8) (spherical
  shell integral with a normalized $(2\pi\sigma^2)^{-3/2} e^{-r^2/(2\sigma^2)}$
  factor; $R>0$ via the standard difference-of-Gaussians radial reduction). The
  quadrature extends the $r$ mesh by $8\sigma$ in $r$ to cover the kernel tail.
  - Paper: (A7)–(A8) with $\sigma$ from (A9), Table II; PDF p. 36.
  - **Not** the expanded derivative forms (A12)–(A13) that the paper actually uses
  for the spin-independent part in the HO diagonalization, and not equivalent to
  convolving with the *same* width as the contact hyperfine when $\sigma$ is
  large: a naive 3D blur of the pointwise $-4\alpha_s/(3r)$ piece can move the
  small-$r$ potential toward less binding on a fixed radial line. The flag
  `appendix_a_smearing` in the parameters file is **off** by default; keep it off
  until the (A12)–(A13) structure (or a momentum/HO-basis path) is implemented.

## Not Yet Implemented

- Full GI effective spin-independent smearing from (A12)–(A13) and/or the paper’s
  HO-basis smearing, replacing the separate experimental (A7)–(A8) convolution
  when `appendix_a_smearing` is enabled.
- Momentum-dependent relativization factors for contact, tensor, vector
  spin-orbit, and scalar spin-orbit interactions.
- Tensor and spin-orbit fine structure.
- Unequal-mass antisymmetric spin-orbit mixing.
- Isoscalar annihilation mixing.
