# Formula Map

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
    (`data/table_ii_parameters.csv`). `scripts/verify_table_ii_toml.py` checks
    they still agree on mapped entries. Table II `relativistic_factors`
    $\epsilon$ values from (A10) are currently applied in the **diagnostic**
    code as **scalar multipliers** on the corresponding terms:
    `epsilon_c`, `epsilon_t`, `epsilon_so_vector`, `epsilon_so_scalar` enter as
    `(1 + epsilon_i)` prefactors on the contact / tensor / vector spin–orbit /
    scalar (Thomas) spin–orbit contributions. This is **not** yet the paper’s
    full relativization in Eq. (A10), which replaces simple `1/m` factors by
    operator factors involving `E_i = sqrt(p^2 + m_i^2)` (and introduces the
    small $\epsilon_i$ as exponents).

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
    finite-difference derivative in `test/runtests.jl`.

- `src/GIModel.jl`: smeared S-wave contact hyperfine shift.
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
    level. Non-S waves and non-FD basis diagnostics keep the perturbative
    expectation path.
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
  - **Not** the expanded derivative forms (A12)–(A13) that the paper actually uses
  for the spin-independent part in the HO diagonalization, and not equivalent to
  convolving with the *same* width as the contact hyperfine when $\sigma$ is
  large: a naive 3D blur of the pointwise $-4\alpha_s/(3r)$ piece can move the
  small-$r$ potential toward less binding on a fixed radial line. The flag
  `appendix_a_smearing` in the parameters file is **off** by default; keep it off
  until the (A12)–(A13) structure (or a momentum/HO-basis path) is implemented.
  - Code toggle: `GIParameters.appendix_a_smearing` gates `potential_diagonal` →
    `smeared_central_values`.
  - Regression test: `test/runtests.jl` checks that `smear_3d_radial` preserves a
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
  - Code toggle: `GIParameters.coulomb_1d_smear` gates `potential_diagonal` →
    `coulomb_1d_smeared_central_values`. This is an older comparator path below
    the closed-form and derivative modes in dispatcher precedence.

- `src/appendix_a_derivative_potential.jl`: **Appendix-A derivative
  proxy** for the Coulomb block on the same FD mesh.
  - Paper context: the Gaussian smearing operator motivates
    $\exp(\nabla^2/(4\sigma^2))G$. The implemented comparator keeps only the
    first correction, $G_{\rm eff}\approx G+\nabla^2G/(4\sigma^2)$, and keeps
    $S(r)=br+c$ pointwise, matching the pragmatic simplification noted around
    (A12)–(A14). This is **not** yet the exact (A12)–(A13) transcription; those
    coefficients still require PDF audit because the OCR extraction is fragile.
  - Numerical implementation: `radial_laplacian_values` uses the uniform FD
    radial Laplacian `f''+2f'/r`; `appendix_a_derivative_central_values` applies
    it to `static_coulomb_G` using the same Table II width from
    `contact_smearing_sigma`.
  - Code toggle: `GIParameters.appendix_a_derivative_g` gates
    `potential_diagonal` through the named dispatcher `central_potential_values`.
    This is now an older comparator path below the closed-form modes, but above
    the raw 3D/1D smearing diagnostics.

- `src/GIModel.jl`: **closed-form Appendix-A central candidates**
  from the expanded web/literature research trail.
  - `:appendix_a_closed_form` evaluates the analytic Gaussian-smearing forms
    used in later GI/MGI implementations:
    $\tilde G(r)=-\sum_k 4\alpha_k\,\mathrm{erf}(\tau_k r)/(3r)$ with
    $\tau_k^{-2}=\sigma^{-2}+\gamma_k^{-2}$, and the corresponding closed-form
    smeared linear $\tilde S(r)$.
  - `:appendix_a_momentum_sandwich` uses the same $\tilde G,\tilde S$ but builds
    the central Coulomb factor as a matrix on the existing FD $p^2$ eigenbasis:
    $G'=A(p)\tilde G(r)A(p)$ with
    $A(p)=\sqrt{1+p^2/(E_1E_2)}$ and $E_i=\sqrt{p^2+m_i^2}$.
  - This is the first implementation path that follows the cross-source
    consensus from later GI/MGI papers. It is still an FD-basis analogue rather
    than the original HO-basis code.
  - Code toggles: `GIParameters.appendix_a_closed_form` and
    `GIParameters.appendix_a_momentum_sandwich`. Precedence is now
    `appendix_a_momentum_sandwich`, then `appendix_a_closed_form`, then the older
    comparator modes.

- `src/reference_state.jl`: **`ReferenceState`**, **`ReferenceStateWithMasses`** (+ **`FineStructureMultiplet`** overloads from CSV rows),
  **`load_reference_spectrum`** (CSV IO).
  Loaded with IO **after** **`sector_comparison.jl`** in **`GIModel.jl`**; **`compute_sector`** / **`compare`** duck-type annotated rows as **`AbstractVector`** with `.constituent_masses` and `.state`.
- `src/masses_from_content.jl`: **`parse_quark_masses`**, **`resolve_constituent_masses`**, **`attach_constituent_masses`** — map string `sector` +
  `quark_content` to **`ConstituentMasses`**; **`attach_constituent_masses`** pairs **[ReferenceState](@ref)** rows from CSV with masses via **`resolve_constituent_masses`**.
  - Paper: Sec. II flavor content; Table II masses.

- `src/sector_solver.jl`: **`RadialChannelKey`**, **`ChannelRadialSolution`**, **`SectorComputation`**, **`solve_sector`**.
- `src/sector_comparison.jl`: **`compute_sector`** batches **`channel_solution`** calls per distinct
  **`RadialChannelKey`** and fills **`SectorComputation.channel_cache`**; **`compare`** maps reference rows to
  cached channels and builds residual **`NamedTuple`** rows; **`write_residual_report`** emits markdown under
  **`docs/residual_reports/`**.
  Structural overview: **`docs/code_architecture.md`**.

- `src/spin_fine_structure.jl`: first-order color-magnetic + Thomas
  (scalar confinement) spin–orbit and an OGE-style tensor term on the FD radial
  mesh; unsmeared Coulomb-kernel tensor term derived from $G(r)=-4\alpha_s(r)/(3r)$
  via $(1/r)\,dG/dr-d^2G/dr^2$ (so running $\alpha_s$ contributes $\alpha_s'(r)$
  and $\alpha_s''(r)$ pieces); Table II
  $\epsilon_t$, $\epsilon_{\rm so(v)}$, $\epsilon_{\rm so(s)}$ are now applied
  as GI-style Hermitian momentum-factor sandwiches around the tensor, vector
  spin–orbit, and scalar Thomas radial kernels; global `k_spin_orbit`,
  `k_tensor` in `[fine_structure]`.
  - Paper: spin-dependent structure around Eqs. (3)–(7) (text), and (A10) (Appendix
    A) for the $\epsilon$ factors. Tensor angular factors for triplet
    $J=L-1,L,L+1$ states use the closed forms
    $-2(L+1)/(2L-1)$, $2$, and $-2L/(2L+3)$, whose $(2J+1)$-weighted
    average vanishes across the triplet multiplet. The spin–orbit implementation
    now follows the *paper text* Eqs. (6)–(7) directly: the color-magnetic piece
    uses $\alpha_s(r)/r^3$ (Eq. (6), **no** $\alpha_s'(r)$ term), while the
    Thomas-precession piece uses $(1/2r)\,dH^{\rm conf}/dr$ (Eq. (7), which
    **does** include $\alpha_s'(r)$ through $d/dr[-\alpha_s(r)/r]$ when running
    $\alpha_s(r)$ is inserted). For unequal masses the exact operator contains
    both symmetric and antisymmetric spin–orbit structures; the current
    diagnostic code keeps only the symmetric contraction into total $L\!\cdot\!S$
    (exact in the equal-mass validation sectors `ccbar`/`bbbar`).
    The tensor radial kernel is still built from the derivative of the pointwise
    running-Coulomb $G(r)$ rather than derivatives of $\tilde G(r)$; the momentum
    sandwich supplies the dominant GI relativization. The next more exact
    replacement is to derive tensor/spin-orbit kernels from the Appendix-A
    smeared $\tilde G(r)$. In the present step we keep the **Coulomb-limit color
    factor** implied by
    $G(r)=-4\alpha_s/(3r)$ and evaluate the full unsmeared running-$\alpha_s$
    kernel $(1/r)\,dG/dr-d^2G/dr^2$, so the tensor prefactor reduces to
    $(4/(3m_1m_2))\langle\alpha_s(r)/r^3\rangle$ in the constant-$\alpha_s$
    limit and includes the corresponding $\alpha_s'(r)$ and $\alpha_s''(r)$
    contributions for the actual running ansatz.
    Radial expectation values treat
    FD eigenvectors as reduced radial functions $u(r)$ with $\int |u|^2\,dr=1$
    (uniform-mesh proxy $\sum |u_i|^2 h = 1$; no additional $4\pi$ factor).
    Scales $k$ bridge the small FD basis to the large HO result and are
    diagnostic, not paper refits.
  - Audit hook: `fine_structure_components(...)` returns the `spin_orbit` and
    `tensor` contributions (and their sum), splits the spin–orbit term into
    `spin_orbit_vector` (color-magnetic) and `spin_orbit_thomas` (scalar
    confinement / Thomas), and also exposes the underlying radial expectation
    values `I_cm`, `I_tp`, and `I_tk` so reports can separate “radial integral”
    effects from angular/mass prefactors.
  - Regression tests in `test/runtests.jl`: triplet tensor/L·S angular sum rules
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
  - Unequal-mass scope note: the current implementation contracts the spin–orbit
    operator into total `L·S` with a symmetric mass prefactor (an “equal-share”
    convention). This is exact for equal-mass $q\bar q$ (the primary validation
    target `ccbar`/`bbbar`), but is only diagnostic for unequal-mass channels
    until the paper’s antisymmetric spin–orbit term and same-`J` mixing
    machinery are implemented.
  - Diagnostic same-`J` mixing hook: `spin_orbit_mixing_components(...)` now
    computes the antisymmetric spin–orbit off-diagonal matrix element for
    `^1L_L`/`^3L_L`, and `same_j_mixing(...)` diagonalizes the resulting `2x2`
    mass matrix. This is deliberately reported as a convention-sensitive
    diagnostic layer rather than folded into `compare` by default; equal-mass
    channels are regression-tested to give zero off-diagonal mixing, and the
    heavy-quarkonium diagnostic report prints the `1P` result for `ccbar`,
    `bbbar`, and the Fig. 9 `bcbar` panel.

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
  solution actually uses). The named central dispatcher also exposes
  `:coulomb_1d` and `:appendix_a_derivative_g` comparison modes so these choices
  can be measured without editing solver internals.
- Smeared **contact** hyperfine and first-order **fine structure** on the same
  $u(r)$. With `contact_momentum_sandwich=true` and
  `fine_structure_momentum_sandwich=true`, the Table II $\epsilon$ values enter
  as GI-style momentum-factor sandwiches around the radial operators. With
  `fine_structure_smeared_kernels=true`, the vector spin-orbit and OGE tensor
  kernels use derivatives of the closed-form smeared Coulomb $\tilde G(r)$, and
  the Thomas term uses derivatives of $\tilde G(r)+\tilde S(r)$.

**What “done” should look like for the spin-independent sector**

- Replace or strictly **bracket** the pointwise central potential with the
  **paper’s Appendix A** effective spin-independent operator: either implement the
  **(A12)–(A13)** structure (or equivalent) on the FD mesh, or reproduce the
  paper’s **HO-basis** construction and map to observables we can compare to
  Fig. 6 / Fig. 8. Until then, the **common mass offset** seen in
  `heavy_quarkonium_diagnostics.md` is expected to be dominated by this gap, not
  by retuning `k_spin_orbit` / `k_tensor`.
- The equal-mass spin-dependent operators now use the same closed-form smeared
  $G(r)$ and confinement $S(r)$ derivatives as the active central path. The
  remaining spin-side gaps are unequal-mass antisymmetric spin-orbit terms,
  same-`J` tensor mixing, and comparison against the paper's perturbative
  ordering in the HO basis.

**Reference row lock-in:** Table II inputs are checked against
`data/table_ii_parameters.csv` via `scripts/verify_table_ii_toml.py`. Reference
spectrum rows are checked for schema via `scripts/validate_reference_spectra.py`.

**Paper navigation for Appendix A:** see `docs/appendix_a_from_paper.md` (equation
labels and PDF pages). The runtime flag for which central path is active is
summarized by `GIModel.central_potential_path` in `src/appendix_a_status.jl`.

## Not Yet Implemented

- Full GI effective spin-independent smearing from (A12)–(A13) and/or the paper’s
  HO-basis smearing, replacing the separate experimental (A7)–(A8) convolution
  and the first-term `appendix_a_derivative_g` proxy when those comparison modes
  are enabled.
- Exact paper-order validation of the momentum-factor sandwiches in the spin
  couplings against the original HO-basis perturbation workflow.
- Unequal-mass antisymmetric spin–orbit and tensor off-diagonal mixing
  (perturbative in the text).
- Isoscalar annihilation and explicit $n\bar n$—$s\bar s$ large mixings
  (Table III). The current debug target is the pseudoscalar `^1S_0` block:
  `docs/residual_reports/pseudoscalar_annihilation_audit.md` shows that a
  positive rank-one radial/flavor annihilation update maps the current unmixed
  `[1n, 1s, 2n, 2s]` masses onto the GI isoscalar pseudoscalar masses, with
  Table-III/P1-like eigenvectors. This is a calibrated diagnostic/control, not
  the literal paper P1 implementation. Paper P1 must still implement Eq. (16)
  with the Eq. (18a) replacement, including the stated `A_eta p` exponential and
  perturbative `alpha_s` product. Paper P2 must still implement Eq. (18b) as a
  mass-dependent pole problem, with non-orthogonal poles expected by the paper.

## Formula-Audit Gate Before P1/P2 Coding

The next implementation pass should use the vision-OCR reference plus crop/PDF
audit of the main formulas. Treat
`paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` as the searchable checkpoint,
then update this map with one status per formula: exact, approximate, missing,
or suspect.

Minimum audit list:

- Hamiltonian pieces and sign/normalization conventions.
- Running `alpha_s` form and units.
- Smearing and Appendix A momentum-dependent `m/E` factors.
- Contact hyperfine ordering and any nonperturbative treatment of S waves.
- Eq. (16)/(17) annihilation normalization and wavefunction factor.
- Eq. (18a)/(18b) pseudoscalar replacements and parameter values.
- Table II/III parameters, units, and model labels.

Only after that audit should dispatch-controlled model choices be added for
`NoAnnihilation`, the existing calibrated diagnostic, literal paper P1, and
literal paper P2.
