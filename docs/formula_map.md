# Formula Map

> **Detail layer.** *Which* code realizes each paper unit — with its status,
> tests, and page image — is the machine-checked manifest/dashboard
> (`docs/paper_manifest/*.toml`; `cd docs && make dashboard`). This file is the
> complementary **how-it's-computed** ledger: the actual numerical conventions,
> closed-form kernels, normalization choices, and diagnostic-vs-active guardrails
> behind those code pointers. Nothing here is duplicated by the manifest, and the
> manifest's `code` fields point back into the same `src/` symbols documented below.

This file maps implemented local code to the Godfrey-Isgur paper. It is also a
guardrail: when the code is only a diagnostic approximation, say so here.

Bundling: Julia package **GIModel** in `Project.toml`, module [`src/GIModel.jl`](../src/GIModel.jl), with `include`'d helpers in `src/`.

## Implemented Baseline

- `src/GIModel.jl`: semirelativistic kinetic operator
  `sqrt(p^2 + m_1^2) + sqrt(p^2 + m_2^2)`.
  - Paper anchor: Eq. (1b), PDF page 2.
  - Numerical implementation: `radial_grid` (uniform $r$ mesh), `p2_operator`
    (FD Laplacian with $L(L+1)/r^2$), and `sqrt_kinetic_matrix` on the $p^2$
    eigenbasis; see `GIModel.jl` for the discrete stencil and matrix square
    root.

- `src/GIModel.jl`: central spin-independent potential
  `b r - 4 alpha_s(r) / (3 r) + c`.
  - Paper anchor: nonrelativistic orientation around Eqs. (2)-(3), PDF pages
    2-3.
  - Parameters: `data/parameters.provisional.toml`, copied from Table II
    (`data/table_ii_parameters.csv`). `test/data_validation.jl` checks
    they still agree on mapped entries. Table II `relativistic_factors`
    $\epsilon$ values are loaded into `GIParameters`; the active
    contact and fine-structure paths use GI-style Hermitian momentum-factor
    sandwiches when `contact_momentum_sandwich` and
    `fine_structure_momentum_sandwich` are enabled. Legacy scalar
    `(1 + epsilon_i)` paths remain available for diagnostics and regression
    comparisons.

- `src/GIModel.jl`: running Coulomb ansatz
  `alpha_s(r) = sum_k alpha_k erf(gamma_k r)`.
  - Paper anchor: Eq. (12), Eq. (13), and Fig. 2 caption, PDF page 3.
  - Coefficients used from Fig. 2 caption:
    `alpha_k = (0.25, 0.15, 0.20)` and momentum-space denominators
    `(1, 10, 1000) GeV^2`, giving
    `gamma_k = (1/2, sqrt(10)/2, sqrt(1000)/2) GeV`.
  - Numerical implementation: `SpecialFunctions.erf` in `alpha_s_r`
    (`src/running_coupling.jl`); derivatives used by the Coulomb spin–orbit
    piece (`dV/dr`) are computed via the analytic derivative of the same profile
    (`erf_prime`)
    to avoid mixed conventions. In code, for the Coulomb central piece
    $V_G(r)=-4\alpha_s(r)/(3r)$ we use
    $$
      \frac{dV_G}{dr}=\frac{4\alpha_s(r)}{3r^2}-\frac{4\alpha_s'(r)}{3r},
    $$
    implemented as `GIModel.dV_coul_central_dr` in
    `src/spin_fine_structure.jl` and regression-tested against a
    finite-difference derivative in `test/spin_kernels.jl`.

- `src/contact_hyperfine.jl`: smeared contact hyperfine, in every `L`.
  - Paper anchor: color hyperfine term around Eq. (4), PDF page 2, and
    smearing discussion in Appendix A, PDF pages 36-37.
  - Current status: active GI-style implementation. It uses the Table II
    `sigma0` and `s` smearing parameters with the standard GI mass-dependent
    Gaussian width form, and applies the contact momentum factor as a Hermitian
    sandwich around the smeared contact kernel:
    $B_c(p)V_c(r)B_c(p)$ with
    $B_c=(m_1m_2/E_1E_2)^{1/2+\epsilon_c}$. The legacy diagonal
    `(1+epsilon_c)` implementation remains available as
    `contact_hyperfine_shift`.
  - Ordering convention: in the finite-difference path, `compare` now follows
    the paper's first diagonalization more closely for S waves by diagonalizing
    the central S-wave Hamiltonian plus the contact operator in fixed
    multiplicity sectors. The reported `contact_shift_GeV` is therefore the
    nonperturbative level displacement relative to the spin-independent central
    level. The operator is also part of every L>0 fixed-sector Hamiltonian
    (FD and HO): the smeared kernel is nonzero at r>0, and GI Eqs. (23)-(26)
    keep its S in all P-wave levels. Until 2026-09 it was S-wave-only, which
    spoiled the singlet-triplet gaps behind the published mixing angles; see
    `GIPaper/docs/investigations/mixing_composition_investigation.md`.
  - Smearing width: Appendix A (A9) is implemented as
    $$
      \sigma^2(m_1,m_2)=\sigma_0^2\left(\tfrac12+\tfrac12\left(\frac{4m_1m_2}{(m_1+m_2)^2}\right)^4\right)
      + s^2\left(\frac{2m_1m_2}{m_1+m_2}\right)^2,
    $$
    via `contact_smearing_sigma` (Table II gives $\sigma_0$ and $s$).
  - Normalization convention: FD eigenvectors are treated as reduced radial
    wavefunctions $u(r)$ with $\int |u|^2\,dr = 1$; the smeared 3D delta kernel
    is normalized to $\int d^3r\,\delta_\sigma(r)=1$ with
    $\delta_\sigma(r)=(\sigma^3/\pi^{3/2})e^{-\sigma^2 r^2}$ (Table II $\sigma$
    has units GeV, so $r$ is treated as GeV$^{-1}$). In code this kernel is
    `GIModel.delta_sigma_3d(r, σ)` and is regression-tested to satisfy
    $4\pi\int r^2 \delta_\sigma(r)\,dr=1$. The S-wave radial expectation then
    uses $\int |u|^2 \alpha_s(r)\,\delta_\sigma(r)\,dr$ (no extra $4\pi$ factor).
    Guardrail: `GIModel.physical_u_norm(r, h, u)` validates that the supplied
    `h` matches a **uniformly spaced** `r` mesh and throws on mismatch, so
    expectation-value conventions cannot silently drift.

- `src/GIModel.jl`: **experimental** 3D isotropic Gaussian smearing
  of pointwise $G(r)$ and $S(r)$ as in the structure of (A7)–(A8) (spherical
  shell integral with a normalized $\rho(\Delta r)=(\sigma^3/\pi^{3/2})
  e^{-\sigma^2 \Delta r^2}$ kernel; $R>0$ via the standard
  difference-of-Gaussians radial reduction). The quadrature extends the $r$
  mesh by $8/\sigma$ in $r$ (so the tail is $\sim e^{-64}$) to cover the kernel
  support.
  - Mesh note: the default FD grid from `radial_grid` starts at `r=h` (not
    `r=0`), so the explicit $R\to 0$ limiting branch inside `smear_3d_radial`
    is only exercised if a caller provides a mesh including `r=0` (kept as a
    numerical guardrail for future Appendix A work).
  - Paper: (A7)–(A8) with $\sigma$ from (A9), Table II; PDF p. 36.
  - **Not** the closed-form (A12)–(A14) branch used by the active central path:
  a naive 3D blur of the pointwise $-4\alpha_s/(3r)$ piece can move the small-$r$
  potential toward less binding on a fixed radial line. This is an explicitly selected diagnostic comparator. The shipped preset
  selects `appendix_a_momentum_sandwich`.
  - Code toggle: `central = AppendixASmearing3D()` selects `potential_diagonal` →
    `smeared_central_values`.
  - Regression test: `test/central_potentials.jl` checks that
    `smear_3d_radial` preserves a
    constant function once the mesh is extended by $8/\sigma$ (the same tail
    coverage used in `smeared_central_values`), so any stray $4\pi$ or kernel
    prefactor drift becomes a loud failure.

- `src/radial_1d_coulomb_smear.jl`: **diagnostic** 1D Gaussian
  renormalization of the pointwise Coulomb piece $G(r)$ on the radial mesh, with
  $S(r)=br+c$ kept pointwise.
  - Paper context: shares the Table II $\sigma_0$, $s$ mass-dependent width family
    (A9) with the contact term, but this is **not** the 3D convolution (A7)–(A8),
    and **not** the paper’s spin-independent effective form (A12)–(A13).
  - Numerical implementation: `convolve_1d_gaussian_same_length` computes a
    local Gaussian-weighted average in $r$ using weights $\propto
    e^{-\sigma^2(r-r')^2}$ (Table II $\sigma$ as in the 3D kernel), and explicitly
    renormalizes by the 1D weight sum (so constants are preserved on a uniform
    grid). There are no $4\pi$ volume factors in this construction; it is a
    mesh-local proxy used to bracket “pointwise” vs “smeared” sensitivity in
    heavy-quarkonium diagnostics.
  - Code toggle: `central = Coulomb1DSmearing()` selects `potential_diagonal` →
    `coulomb_1d_smeared_central_values`. Central methods are mutually exclusive typed dispatch choices.

- `src/appendix_a_derivative_potential.jl`: **Appendix-A derivative
  proxy** for the Coulomb block on the same FD mesh.
  - Paper context: the Gaussian smearing operator motivates
    $\exp(\nabla^2/(4\sigma^2))G$. The implemented comparator keeps only the
    first correction, $G_{\rm eff}\approx G+\nabla^2G/(4\sigma^2)$, and keeps
    $S(r)=br+c$ pointwise, matching the pragmatic simplification noted around
    (A12)–(A14). This is now an older comparator: the closed-form (A12)–(A14)
    branch is the primary Appendix-A central source, and this first-term
    expansion remains useful only as a sensitivity check.
  - Numerical implementation: `radial_laplacian_values` uses the uniform FD
    radial Laplacian `f''+2f'/r`; `appendix_a_derivative_central_values` applies
    it to `static_coulomb_G` using the same Table II width from
    `contact_smearing_sigma`.
  - Code toggle: `central = AppendixADerivativeG()` selects
    `potential_diagonal` through the named dispatcher `central_potential_values`.
    This diagnostic must be selected explicitly; there is no dispatcher precedence.

- `src/GIModel.jl`: **closed-form Appendix-A central candidates**
  from the checked Appendix-A markdown source.
  - `:appendix_a_closed_form` evaluates the analytic Gaussian-smearing forms
    from (A12)-(A14):
    $\tilde G(r)=-\sum_k 4\alpha_k\,\mathrm{erf}(\tau_k r)/(3r)$ with
    $\tau_k^{-2}=\sigma^{-2}+\gamma_k^{-2}$, and the corresponding closed-form
    smeared linear $\tilde S(r)$.
  - `:appendix_a_momentum_sandwich` uses the same $\tilde G,\tilde S$ but builds
    the central Coulomb factor as a matrix on the existing FD $p^2$ eigenbasis:
    $G'=A(p)\tilde G(r)A(p)$ with
    $A(p)=\sqrt{1+p^2/(E_1E_2)}$ and $E_i=\sqrt{p^2+m_i^2}$.
  - Both FD and paper-order HO solvers implement this central prescription in
    their native momentum-squared representation.
  - Select `AppendixAClosedForm()` or `AppendixAMomentumSandwich()` through
    the immutable `GIParameters` constructor’s `central` keyword; the latter
    is the shipped preset. Methods do not form a precedence chain.

- `GIPaper/src/reference_state.jl`: `ReferenceState` and
  `load_reference_spectrum` own CSV IO.
- `GIPaper/src/reference_meson.jl`: `reference_meson` maps supported paper
  sector/content labels to core `Meson` objects and fails loudly otherwise.
  - Paper: Sec. II flavor content; Table II masses.
- `src/sector_solver.jl`: `RadialChannelKey`, `ChannelRadialSolution`,
  `SectorComputation`, and `solve_sector` own the native radial-solution cache.
- `src/spectrum.jl`: `central_spectrum`, `fixed_spectrum`,
  `add_intra_meson_mixing`, and `compute_spectrum` implement the typed spectrum
  pipeline.
- `GIPaper/src/comparison.jl` and `residual_report.jl`: `compare_reference`
  maps paper rows to core spectra and `write_residual_report` emits markdown.
  Structural overview: **`docs/code_architecture.md`**.

- `src/contact_hyperfine.jl` and `src/spin_fine_structure.jl`: literal A15-A16
  contact, vector/scalar spin–orbit, and tensor operators. Each `11`, `22`, or
  `12` term has its own smeared kernel, mass denominator, and two-sided
  $(m_\alpha m_\beta/E_\alpha E_\beta)^{1/2+\epsilon_i}$ factor. HO builds these
  as native matrices; FD builds the same decomposition on its native grid.
  - The contact density is the analytic Laplacian of A12,
    $\sum_k\alpha_k\delta_{\tau_k}(r)$.
  - Diagonal triplets use
    $\langle\mathbf S_1\!\cdot\!\mathbf L\rangle=
    \langle\mathbf S_2\!\cdot\!\mathbf L\rangle=
    \langle\mathbf S\!\cdot\!\mathbf L\rangle/2$; unequal-mass same-$J$
    mixing uses the corresponding antisymmetric difference.
  - A16 contains only the scalar-confinement `11` and `22` derivatives.
  - `tensor_triplet_LJ` stores conventional Pauli-$S_{12}$ angular elements;
    the A15 spin bracket is therefore $S_{12}/12$.
  - No fitted spin-strength multiplier exists. `FineStructure` is only a master
    switch; Table-II $\epsilon_c$, $\epsilon_t$, $\epsilon_{so(v)}$, and
    $\epsilon_{so(s)}$ supply the paper's relativization.
  - `fine_structure_components(...)` exposes `I_vector_11`, `I_vector_22`,
    `I_vector_12`, `I_scalar_11`, `I_scalar_22`, and `I_tk` together with the
    assembled contributions. Radial expectation values use the shared
    `RadialWave` normalization contract.
  - Regression tests in `test/fine_structure.jl` and
    `test/spin_kernels.jl`: triplet tensor/L·S angular sum rules
    and Coulomb $d\alpha_s/dr$ consistency; the OGE tensor kernel
    $K(r)=(1/r)\,dG/dr-d^2G/dr^2$ is also checked against finite-difference
    derivatives of $G(r)=-4\alpha_s(r)/(3r)$; `fine_structure_split` checks that
    S-waves and P singlets have zero first-order fine-structure shift and that
    the $1P$ triplet $J=0,1,2$ splittings are not all identical (finite
    $r$-space on the diagnostic mesh). Both fine-structure and smeared-contact
    expectation values are tested to be **scale-invariant** under rescaling of
    the FD eigenvector, enforcing the documented convention that the code treats
    solver eigenvectors as reduced radial functions $u(r)$ that are normalized
    via $\int |u|^2 dr=1$ before computing expectation values; the contact test
    also checks the expected singlet/triplet ratio from `spin_dot`.
  - Unequal-mass scope note: diagonal fine structure still contracts the
    spin-orbit operator into total `L·S` with a symmetric mass prefactor. The
    antisymmetric part is now applied afterward for open-flavor same-`J`
    `^1L_J`/`^3L_J` pairs when both partner rows are present.
  - Same-`J` mixing path: `spin_orbit_mixing_components(...)` computes the
    antisymmetric spin–orbit off-diagonal matrix element for `^1L_L`/`^3L_L`,
    `same_j_mixing(...)` diagonalizes the resulting `2x2` mass matrix, and
    `compare(...)` assigns mixed eigenvalues while reporting the off-diagonal
    element, angle, and singlet/triplet components. Equal-mass channels are
    regression-tested to give zero off-diagonal mixing.

## Appendix A: paper method vs this codebase (completion target)

The 1985 paper does **not** stop at a pointwise $V(r)$ on a radial line. Appendix
A explains (i) **relativistic smearing** so spin-dependent singularities become
well-defined operators, and (ii) a **harmonic-oscillator (HO) diagonalization**
path for the spin-independent problem, with expanded effective forms of the
potentials. In the prose (e.g. around the discussion of Eq. (2)), Godfrey and
Isgur state that the coordinate $r$ is smeared at the scale of inverse quark
masses and that detailed smearing is “relegated to Appendix A.”

**What we implement today**

- A native **harmonic-oscillator** implementation of the paper route and an
  independent **finite-difference** comparator, both with semirelativistic
  $\sqrt{p^2+m^2}$ kinetics and the closed-form Appendix-A central operator
  `A(p)G~A(p)+S~` in their own $p^2$ representations.
- Smeared contact, symmetric spin-orbit, and diagonal tensor matrices assembled
  with the central operator in one fixed-$(L,S,J)$ Hamiltonian before
  diagonalization. The Table II $\epsilon$ values enter through native
  momentum-factor sandwiches; later tensor and antisymmetric spin-orbit blocks
  consume the resulting native waves.
- An adaptive HO beta/basis controller with a recorded convergence certificate.
  Pointwise and alternate-smearing central methods remain explicitly named
  standalone comparators.

**What “done” should look like for the spin-independent sector**

- Use converged native HO as the paper path and FD as an independent comparator.
  The focused FD/HO central comparison is complete in
  `docs/residual_reports/appendix_a_ho_comparison.md`: the active central
  operator agrees between bases at the sub-MeV level over the audited channels.
  PA-18 completed the A15-A16 normalization and paper-path certification;
  `fd_comparator_convergence.md` separately certifies the modern FD route.
- The equal-mass spin-dependent operators now use the same closed-form smeared
  $G(r)$ and confinement $S(r)$ derivatives as the active central path. The
  open-flavor antisymmetric spin-orbit block and triplet same-`J` tensor blocks
  dispatch on both native wave representations with no fitted bridge factor.

**Reference row lock-in:** Table II inputs are checked against
`data/table_ii_parameters.csv` via `test/data_validation.jl`. Reference
spectrum rows are checked for schema by the same command.
Formula-audit checkpoint: Table II `epsilon_so_scalar` was rechecked against
`paper/vision_ocr/page_images/page-005.png`; the paper value is
`epsilon_so(S)=+0.055`, now reflected in both the CSV and TOML inputs.

**Paper navigation for Appendix A:** see `docs/appendix_a_from_paper.md` and
`docs/appendix_a_equation_audit.md` (equation labels, local source status, and
next stages). The runtime flag for which central path is active is summarized
by `GIModel.central_potential_path` in `src/appendix_a_status.jl`.

## Formula Audit Ledger

Treat `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` as searchable. For
Appendix A, use the checked split markdown pages and
`docs/appendix_a_equation_audit.md` as the local equation-source ledger.

| Paper item | Status | Implementation/readout |
| --- | --- | --- |
| Eq. (1a)-(1b), semirelativistic Hamiltonian | FD and HO analogues implemented | `relativistic_hamiltonian` uses `sqrt(p^2+m_1^2)+sqrt(p^2+m_2^2)` on the finite-difference `p^2` operator. Passing an `OscillatorSolver` selects the finite oscillator-basis analogue; native HO is the paper path and FD is the independent modern comparator. |
| Eq. (3), spin-independent color Coulomb plus linear confinement | exact color-singlet sign/normalization for pointwise limit | `central_potential = b*r - 4*alpha_s(r)/(3r) + c`; the color-singlet factor turns Eq. (3) into this form. |
| Eq. (4), contact hyperfine | reproduced in the Appendix-A prescription | `smeared_contact_kernel` is the analytic Laplacian of A12, `sum(alpha_k delta_tau_k)`, enclosed by the post-A14 `m/E` sandwich and diagonalized with the central Hamiltonian. |
| Eq. (4), tensor hyperfine | reproduced through A15 | Native HO and FD use derivatives of closed-form smeared `G~`, the literal `S12/12` contraction, and complete same-`J` `L=J-1`/`L=J+1` radial mixing blocks. Angular factors and matrix symmetry are tested. |
| Eq. (6), vector spin-orbit | reproduced through A15 | Separate `G11`, `G22`, and `G12` sandwiches produce diagonal symmetric and unequal-mass antisymmetric contractions. |
| Eq. (7), Thomas/scalar spin-orbit | reproduced through A16 | Separate scalar-confinement `S11` and `S22` sandwiches produce diagonal symmetric and unequal-mass antisymmetric contractions; no Coulomb derivative is duplicated here. |
| Eq. (12)-(13) and Fig. 2, running `alpha_s` | exact for fitted GI ansatz | `alpha_s(Q^2)` coefficients map to `alpha_s(r)=sum alpha_k erf(gamma_k*r)` with `gamma=(1,sqrt(10),sqrt(1000))/2` GeV. Derivatives are regression-tested. |
| Eq. (14), staged diagonalization | implemented and certified | One production algorithm assembles and diagonalizes complete fixed-sector Hamiltonians, with native matrices and waves supplied by solver dispatch. HO refines beta/basis to a recorded 0.1 MeV certificate. Generic `MixingBlock`/`MixingResult` stages apply antisymmetric spin-orbit, tensor, and flavor-annihilation transformations and expose recursively composed physical states to all consumers. |
| Eq. (16), general annihilation matrix element | implemented and integrated | OCR/page pass confirms the `4*pi*(2L+1)`, `alpha_s(M_i^2)alpha_s(M_j^2)`, `S_L(Psi_j)S_L(Psi_i)`, and `1/(m_i m_j)` factors. `isoscalar_general_annihilation_solution` implements the general channel block with the `(alpha_i alpha_j/pi^2)^(n/2)` bracket (`n=2`/`3` for `C=+`/`-`). `add_isoscalar_annihilation` stores the resulting shared `MixingResult` in a final two-channel `Spectrum`, and `physical_components` resolves its flavor-tagged native waves. GIPaper retains only reference-row assignment. |
| Eq. (17), `S_L(Psi)` wavefunction factor | implemented for general `L` | OCR/page pass confirms the momentum-space normalized-wavefunction factor `(2*pi)^(-3/2) * integral d^3p/sqrt(4*pi) * Phi_i(p) * (p/E_i)^L * (m_i/E_i)`. `wavefunction_origin_smearing` evaluates this for general `L` through the spectrum wave's native momentum interface. The independent q/s/c/b six-level audit finds converged HO/FD differences no larger than 0.7%; the obsolete coordinate-origin proxy has been removed. |
| Eq. (18a), pseudoscalar P1 replacement | implemented as `:p1` / `PaperP1Annihilation` | OCR/page pass confirms the bracket replacement `A_np*exp(-(m_i^2+m_j^2)/m_eta^2) + (2*pi/3)(ln2-1) alpha_s(M_j^2)alpha_s(M_i^2)/pi^2`. Table III constants are loaded from `[annihilation]` in `data/parameters.provisional.toml`. |
| Eq. (18b), pseudoscalar P2 replacement | implemented as `:p2` / `PaperP2Annihilation` | OCR/page pass confirms the mass-dependent sign-changing term `A_np*(1-(M/M0)^4)*exp(-(m_i^2+m_j^2)/M0^2 - M^4/(4M0^4))` plus the perturbative `(alpha_s(M^2)/pi)^2` term. The implementation solves each pole as a mass-dependent fixed point, so the vectors are not treated as orthogonal. |
| Table II parameters | audited for active solver inputs | CSV/TOML sync is tested. The latest audited correction is `epsilon_so(S)=+0.055`. Remaining low-confidence labels should still be promoted only after image/PDF checks. |
| Table III isoscalar mixings | visible page-11 rows promoted; reproduced by `:table_iii` | The visible page-11 pseudoscalar P1/P2 rows plus `1^3S_1` and `1^3P_2` rows are promoted into clean `mixings.csv`. The general Eq. (16) blocks reproduce the non-pseudoscalar eigenvectors; current residual statistics live only in `table_iii_mixing_audit.md`. Channels not listed in Table III are ideally mixed, as the caption prescribes. |
| Appendix A (A7)-(A9), Gaussian smearing and sigma | exact for contact width and diagnostic convolution | A9 width is implemented and tested; A7-A8 3D convolution is available as a diagnostic path, not the paper's main spin-independent calculation. |
| Appendix A post-A14 spin-dependent momentum factors | native HO and FD | Side exponent is `1/2+epsilon_i`, so `epsilon=0` replaces `1/(m_1m_2)` by `1/(E_1E_2)` after the two-sided sandwich. |
| Appendix A (A12)-(A14), closed-form smeared `G~`, `S~`, `tau_k` | native HO and FD | Closed-form `G~` and `S~` are tested against quadrature/derivatives and enter the HO matrix-element workflow directly. |
| Appendix A (A15)-(A17), effective operators and HO matrix elements | reproduced and certified | Native HO assembles literal contact, vector/scalar spin-orbit, and tensor matrices with exact A17 factorization, diagonalizes complete fixed sectors, then applies full radial mixing blocks. `ho_convergence.md`, `w6_ho_order_validation.md`, and the Table VI/VII reports gate the path. |

## Remaining physics comparisons

- Literal pseudoscalar P1/P2 spectral scoring and the residual eta-prime
  amplitude discrepancy. With HO wavefunctions in Eq. (17) and the corrected
  Eq. (16) gluon-power bracket, the earlier under-mixing diagnosis is
  resolved (see `docs/residual_reports/table_iii_mixing_audit.md`); the
  calibrated P1 control still carries the headline mass scoring until the
  literal modes also score the digitized Fig. 5 masses.
