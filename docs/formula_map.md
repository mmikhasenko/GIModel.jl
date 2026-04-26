# Formula Map

This file maps implemented local code to the Godfrey-Isgur paper. It is also a
guardrail: when the code is only a diagnostic approximation, say so here.

## Implemented Baseline

- `src/GIModel/GIModel.jl`: semirelativistic kinetic operator
  `sqrt(p^2 + m_1^2) + sqrt(p^2 + m_2^2)`.
  - Paper anchor: Eq. (1b), PDF page 2.
  - Numerical implementation: `radial_grid` (uniform $r$ mesh), `p2_operator`
    (FD Laplacian with $L(L+1)/r^2$), and `sqrt_kinetic_matrix` on the $p^2$
    eigenbasis; see `GIModel.jl` for the discrete stencil and matrix square
    root.

- `src/GIModel/GIModel.jl`: central spin-independent potential
  `b r - 4 alpha_s(r) / (3 r) + c`.
  - Paper anchor: nonrelativistic orientation around Eqs. (2)-(3), PDF pages
    2-3.
  - Parameters: `data/parameters.provisional.toml`, copied from Table II
    (`data/table_ii_parameters.csv`). `scripts/verify_table_ii_toml.py` checks
    they still agree on mapped entries. Table II `relativistic_factors`
    $\epsilon_i$ from (A10) are applied in the first-order fine-structure and
    contact terms as noted below.

- `src/GIModel/GIModel.jl`: running Coulomb ansatz
  `alpha_s(r) = sum_k alpha_k erf(gamma_k r)`.
  - Paper anchor: Eq. (12), Eq. (13), and Fig. 2 caption, PDF page 3.
  - Coefficients used from Fig. 2 caption:
    `alpha_k = (0.25, 0.15, 0.20)` and momentum-space denominators
    `(1, 10, 1000) GeV^2`, giving
    `gamma_k = (1/2, sqrt(10)/2, sqrt(1000)/2) GeV`.
  - Numerical implementation: `erf_approx` is used (no external dependencies);
    derivatives used by the Coulomb spin–orbit piece (`dV/dr`) are computed via
    the analytic derivative of the same approximation (`erf_approx_prime`) to
    avoid mixed conventions.

- `src/GIModel/GIModel.jl`: smeared S-wave contact hyperfine shift.
  - Paper anchor: color hyperfine term around Eq. (4), PDF page 2, and
    smearing discussion in Appendix A, PDF pages 36-37.
  - Current status: first diagnostic implementation. It uses the Table II
    `sigma0` and `s` smearing parameters with the standard GI mass-dependent
    Gaussian width form. The full momentum-dependent relativization factors
    are not yet included.
  - Normalization convention: FD eigenvectors are treated as reduced radial
    wavefunctions $u(r)$ with $\int |u|^2\,dr = 1$; the smeared 3D delta kernel
    is normalized to $\int d^3r\,\delta_\sigma(r)=1$, so the S-wave radial
    expectation uses $\int |u|^2 \delta_\sigma(r)\,dr$ (no extra $4\pi$ factor).

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
  - Code toggle: `GIParameters.appendix_a_smearing` gates `potential_diagonal` →
    `smeared_central_values`.

- `src/GIModel/masses_from_content.jl`: map `sector` + first `quark_content`
  segment to constituent $(m_1, m_2)$ for `compare_sector` (unequal-mass channels).
  - Paper: Sec. II flavor content; Table II masses.

- `src/GIModel/spin_fine_structure.jl`: first-order color-magnetic + Thomas
  (scalar confinement) spin–orbit and an OGE-style tensor term on the FD radial
  mesh; unsmeared $\alpha_s(r)/r^3$ proxy for the tensor piece; Table II $\epsilon_t$,
  $\epsilon_{\rm so(v)}$, $\epsilon_{\rm so(s)}$; global `k_spin_orbit`,
  `k_tensor` in `[fine_structure]`.
  - Paper: spin-dependent structure around Eqs. (3)–(7) (text), and (A10) (Appendix
    A) for the $\epsilon$ factors. Tensor angular factors for triplet
    $J=L-1,L,L+1$ states use the closed forms
    $-2(L+1)/(2L-1)$, $2$, and $-2L/(2L+3)$, whose $(2J+1)$-weighted
    average vanishes across the triplet multiplet. The spin-orbit radial
    coefficient uses the standard color-magnetic minus Thomas structure
    proportional to $3G'(r)-S'(r)$. The tensor radial term is deliberately
    unsmeared for now because applying the broad contact width to $1/r^3$
    overdamps the P/D splittings; the correct next replacement is the
    derivative of the Appendix A smeared $G(r)$. Radial expectation values treat
    FD eigenvectors as reduced radial functions $u(r)$ with $\int |u|^2\,dr=1$
    (uniform-mesh proxy $\sum |u_i|^2 h = 1$; no additional $4\pi$ factor).
    Scales $k$ bridge the small FD basis to the large HO result and are
    diagnostic, not paper refits.
  - Regression tests in `test/runtests.jl`: triplet tensor/L·S angular sum rules
    and Coulomb $d\alpha_s/dr$ consistency; `fine_structure_split` checks that
    S-waves and P singlets have zero first-order fine-structure shift and that
    the $1P$ triplet $J=0,1,2$ splittings are not all identical (finite
    $r$-space on the diagnostic mesh).

## Appendix A: paper method vs this codebase (completion target)

The 1985 paper does **not** stop at a pointwise $V(r)$ on a radial line. Appendix
A explains (i) **relativistic smearing** so spin-dependent singularities become
well-defined operators, and (ii) a **harmonic-oscillator (HO) diagonalization**
path for the spin-independent problem, with expanded effective forms of the
potentials. In the prose (e.g. around the discussion of Eq. (2)), Godfrey and
Isgur state that the coordinate $r$ is smeared at the scale of inverse quark
masses and that detailed smearing is “relegated to Appendix A.”

**What we implement today**

- A **finite-difference** radial mesh with semirelativistic $\sqrt{p^2+m^2}$
  kinetics and a **pointwise** spin-independent $V$ (plus the experimental
  `appendix_a_smearing` branch that 3D-blurs pointwise $G$ and $S$ in the spirit
  of (A7)–(A8), but **not** the (A12)–(A13) derivative expansion the HO
  solution actually uses).
- Smeared **contact** hyperfine and first-order **fine structure** on the same
  $u(r)$, with (A10) $\epsilon$ factors and **diagnostic** global $k$ scales.

**What “done” should look like for the spin-independent sector**

- Replace or strictly **bracket** the pointwise central potential with the
  **paper’s Appendix A** effective spin-independent operator: either implement the
  **(A12)–(A13)** structure (or equivalent) on the FD mesh, or reproduce the
  paper’s **HO-basis** construction and map to observables we can compare to
  Fig. 6 / Fig. 8. Until then, the **common mass offset** seen in
  `heavy_quarkonium_diagnostics.md` is expected to be dominated by this gap, not
  by retuning `k_spin_orbit` / `k_tensor`.
- Spin-dependent operators should eventually use the **same** smeared $G(r)$ and
  confinement $S(r)$ as the central sector (the paper ties this together in
  Appendix A and in the discussion of (A10)–(A13)).

**Reference row lock-in:** Table II inputs are checked against
`data/table_ii_parameters.csv` via `scripts/verify_table_ii_toml.py`. Reference
spectrum rows are checked for schema via `scripts/validate_reference_spectra.py`.

## Not Yet Implemented

- Full GI effective spin-independent smearing from (A12)–(A13) and/or the paper’s
  HO-basis smearing, replacing the separate experimental (A7)–(A8) convolution
  when `appendix_a_smearing` is enabled.
- Full $E_i/m_i$ or $(p^2{+}m_i^2)^{1/2}$ momentum dependence in the spin
  couplings (paper’s relativization beyond constant $\epsilon$).
- Momentum-dependent relativization factors for the contact and tensor
  *operators* (beyond the $\epsilon$ factors already applied in the
  `fine_structure` block).
- Unequal-mass antisymmetric spin–orbit and tensor off-diagonal mixing
  (perturbative in the text).
- Isoscalar annihilation and explicit $n\bar n$—$s\bar s$ large mixings
  (Table III).
