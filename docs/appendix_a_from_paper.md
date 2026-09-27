# Appendix A: structure from the 1985 paper (navigation)

This note orients the repository toward the **original** Godfrey–Isgur
implementation of relativistic smearing, following Appendix A of the published
article (PDF pages 35-38); see `docs/appendix_a_equation_audit.md` for the
equation-status ledger.

## Smearing setup

- **(A7)** – Definition of a **smearing function** in three dimensions.
- **(A8)** – The smeared form of a potential:
  `f̃(r) = ∫ d³r' ρ(|r - r'|) f(r')` (schematic: convolution with a normalized
  kernel; see the PDF for the exact definition and notation for `G` and `S`).
- **(A9)** – Universal **σ** and **s** (Table II). In this repo they feed
  `contact_smearing_sigma` and the (experimental) 3D blur
  `smear_3d_radial` / `smeared_central_values`.

## Pointwise And Smeared Central Potentials

- **(A10)** and **(A11)** – Pointwise Coulomb and confinement building blocks:
  `G(r)` and `S(r)`.
- The **Coulomb** and **confinement** building blocks in coordinate space
  (orientation around Eqs. (11)–(13) in the main text) match our `static_coulomb_G`
  and `static_confinement_S` (see `docs/formula_map.md`).
- **(A12)**, **(A13)**, and **(A14)** – Closed-form smeared `G~`, smeared `S~`,
  and `tau_k`. These are implemented by `smeared_coulomb_G_closed`,
  `smeared_confinement_S_closed`, and `appendix_a_closed_central_values`.
- The active central path adds the subsequent Coulomb momentum factor as
  `A(p) G~ A(p)`, assembled natively in either HO or FD representation.

## Effective operators and factorization

- **(A5)** and **(A6)** – Derivational context only. They motivate the later
  prescription, but the local markdown still has spin-label risk and should not
  be used as a direct coefficient source.
- **(A15)** – Contact, vector spin-orbit, and tensor operators are implemented
  literally. The contact kernel is the analytic Laplacian of `G~`; `11`, `22`,
  and `12` mass pairs remain separate; the Pauli-`S12` convention carries the
  required factor `1/12`.
- **(A16)** – Separate scalar-confinement `11` and `22` derivatives generate
  both diagonal and unequal-mass antisymmetric spin-orbit terms.
- **(A17)** – HO position and momentum matrix elements are assembled without a
  mesh, beta/basis are refined variationally, complete fixed sectors are
  diagonalized, and the later tensor/spin-orbit/annihilation blocks consume the
  resulting native waves.

**Repository consequence:** the flag `appendix_a_smearing` remains the older
experimental (A7)–(A8) 3D blur of pointwise `G` and `S`; keep it as a diagnostic
branch. The active central reproduction branch is
`appendix_a_momentum_sandwich`, which uses closed-form `G~`, `S~`, plus
`A(p)G~A(p)`. The original HO algorithm and its post-diagonalization mixing
layer are complete; FD remains an independent modern comparator.

## Quick PDF map

| PDF page (approx.) | Content |
| ---: | --- |
| 36 | (A7)–(A11) |
| 37 | (A12) start, (A9) text |
| 37–38 | (A12)–(A16) expanded potentials |

## Related code

- `smeared_central_values` / `smear_3d_radial` — `src/GIModel.jl`
- `smeared_coulomb_G_closed` / `smeared_confinement_S_closed` —
  `src/smearing_appendix_a.jl`
- `appendix_a_momentum_sandwich_matrix` — `src/hamiltonian.jl`
- `appendix_a_derivative_central_values` —
  `src/appendix_a_derivative_potential.jl`
- `src/appendix_a_status.jl` — explicit status (which path is active)
- `docs/appendix_a_equation_audit.md` — equation-source audit and next-stage
  ledger
- Appendix-A diagnostics now live in Julia tests and `docs/paper_gap_ledger.md`
  instead of tracked one-off report scripts.
