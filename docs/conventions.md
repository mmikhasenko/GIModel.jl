# Conventions

This document will record the conventions used by the cleaned data and model
implementation.

## Current Code Conventions

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

## Sector naming and reference CSVs

`data/reference_spectrum_*.csv` use a `sector` string per file, for example:
`charmonium`, `bottomonium`, `charmed`, `b_flavored`, `strange`, `isovector`,
`isoscalar`. These are labels for the reference rows and the residual reports;
the solver’s `compare_sector` matches them to Table II quark flavors via
`quark_content` and `parse_quark_masses` in `src/GIModel/`.

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
