# Conventions

> **Detail layer.** For paper-unit status and where each is realized, use the
> manifest/dashboard (`docs/paper_manifest/*.toml`; `cd docs && make dashboard`).
> This file is the standing **code-and-data conventions** reference (parameter
> loading, sector naming, quark ordering) that the whole codebase relies on.

This document will record the conventions used by the cleaned data and model
implementation.

## Current Code Conventions

### Parameter and mass loading

Solver switches (**`GIParameters`**) and flavor masses (**`QuarkMassTable`**) come from the same
TOML file but different tables: see **`load_parameters`**, **`load_quark_masses`**, and the
usual combined **`load_parameters_and_quark_masses`** in `src/quark_mass_table.jl`. Masses are
**not** fields on **`GIParameters`**. **`ReferenceState`** / **`ReferenceStateWithMasses`** are defined in
**`reference_state.jl`** (first IO include in **`GIModel.jl`**), including **`load_reference_spectrum`**. Reference CSV rows become **`ReferenceState`** via that loader; **`attach_constituent_masses`** (**`masses_from_content.jl`**) produces **`ReferenceStateWithMasses`** for **`compute_sector`** /
**`compare`**, calling **`resolve_constituent_masses`** (**`masses_from_content.jl`**) with sector/content strings.
Full call graph: **`docs/code_architecture.md`**.

- Energies and masses are in GeV internally.
- Distances are in GeV^-1.
- Reference spectrum CSV masses are read in GeV.
- Finite-difference bound-state eigenvectors are treated as reduced radial
  wavefunctions `u(r)` on a uniform mesh (so `ψ(r,Ω) = u(r)/r * Y_{LM}(Ω)`), with
  physical normalization `∫ |u(r)|^2 dr = 1` (discrete proxy `∑ |uᵢ|² h = 1`).
- `L` is stored as spectroscopic letters `S`, `P`, `D`, `F`, `G` and mapped to
  orbital angular momentum `0, 1, 2, 3, 4`.
- Multiplicity is `2S+1`; the current contact hyperfine term supports
  singlets (`1`) and triplets (`3`).
- For equal-mass quarkonia, the baseline solver predicts one radial level per
  `(n, L)` before spin shifts.
- The current residual reports compare individual Fig. 6 and Fig. 8 labels
  directly to the baseline prediction plus available spin shifts.
- Central potential “path” is parameterized and reported by
  `GIModel.central_potential_path(params)`:
  pointwise `V(r)`, optional `coulomb_1d_smear` (1D Gaussian renormalization of
  Coulomb $G(r)$ only), experimental `appendix_a_smearing` (3D (A7)–(A8)
  convolution of pointwise $G$ and $S$), `appendix_a_derivative_g` (first
  derivative-smearing proxy for $G$ with pointwise $S$),
  `appendix_a_closed_form` (closed-form (A12)-(A14) `G~+S~`), or
  `appendix_a_momentum_sandwich` (active central candidate
  `A(p)G~A(p)+S~`). The 1D, 3D, and derivative paths are comparison branches;
  the closed-form and momentum-sandwich paths are the current Appendix-A central
  reproduction candidates. See `docs/appendix_a_equation_audit.md`.

## Sector naming and reference CSVs

`data/reference_spectrum_*.csv` use a `sector` string per file, for example:
`charmonium`, `bottomonium`, `charmed`, `b_flavored`, `strange`, `isovector`,
`isoscalar`. These are labels for the reference rows and the residual reports.
**`attach_constituent_masses`** (**`masses_from_content.jl`**) wraps each **`ReferenceState`** using **`resolve_constituent_masses`**
(**`masses_from_content.jl`**) with `parse_quark_masses`, passing `String(state.sector)` and `String(state.quark_content)`.
Then **`compute_sector`** and **`compare`** (**`sector_comparison.jl`**; radial cache types in **`sector_solver.jl`**) map each annotated row to
Table II masses via the attached **`ConstituentMasses(m1, m2)`** (12-digit rounding) and orbital letter
`L` in a **`RadialChannelKey`**; identical keys share one cached radial
finite-difference solve before **`compare`** attaches level `n`, hyperfine, and
fine-structure shifts using **`RadialWaveOnUniformMesh`** for one column of the
cached eigenvectors.

- Equal-mass quarkonia: `c cbar`, `b bbar` in the `quark_content` column.
- Heavy-light: semicolon lists such as `c ubar; c dbar` for charmed, similarly
  for $b$-flavored rows.

## Required topics (stable vs future)

Addressed in this file or `docs/formula_map.md` for the diagnostic build:

- Spectroscopic notation, units, and radial $u(r)$ convention (see **Current
  Code Conventions** above). Tensor and spin–orbit operator forms are mapped in
  `docs/formula_map.md` to Eqs. (3)–(7) and Appendix A context.

Pending or only partially specified in code (do not over-interpret current
residuals):

- Full mixed-state mass-matrix labeling for $^1L_J$ / $^3L_J$ sectors,
- isoscalar ($I=0$) and hidden-flavor annihilation as in the full GI treatment.

## Quark order

`quark_content` strings follow the paper-style flavor listing; the parser in
`masses_from_content.jl` maps the first quark/anti-quark token pair to $(m_1, m_2)$
for the reduced Hamiltonian. If a row lists multiple $q\bar{q}$ pairs, the first
valid pair is used unless extended logic is added.
