# Appendix A Equation Audit

This ledger records the current equation source status for Appendix A. The best
current OCR-derived text is `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md`;
the older split files under `paper/vision_ocr/pages/` are archived historical
OCR output and remain useful for provenance comparison only. Implementation
claims below should be updated only after checking against the PDF or page
images.

Historical page references used when this ledger was first written:

- `paper/vision_ocr/pages/page-036.md` for (A1)-(A9).
- `paper/vision_ocr/pages/page-037.md` for (A10)-(A15).
- `paper/vision_ocr/pages/page-038.md` for (A16)-(A17).

## Source Status

| Item | Local source status | Code status | Audit note |
| --- | --- | --- | --- |
| (A1)-(A4) scattering setup | clear enough for context | not implemented directly | These equations motivate the effective potential construction, but the code uses the later prescription rather than evaluating the scattering kernel. |
| (A5) vector/Coulomb effective kernel | context only | not implemented directly | The bracket structure is readable, but the spin labels remain risky as a code source: the text defines both `S_q` and `S_{\bar q}`, while the displayed expression uses `S_q` in every spin slot. Do not use this equation as a coefficient source without a dedicated manual check. |
| (A6) scalar/confinement effective kernel | context only | not implemented directly | Same caveat as (A5). The line-broken product is useful for deriving the paper's qualitative prescriptions, not for direct transcription. |
| (A7) Gaussian smearing kernel | clear | implemented for contact and diagnostic 3D smearing | The normalization matches `delta_sigma_3d` and `smear_3d_radial`; tests cover normalization/constant preservation. |
| (A8) smeared radial potential definition | clear | implemented as diagnostic convolution and closed-form central path | The direct convolution path remains diagnostic; the active central path uses the closed-form results below. |
| (A9) mass-dependent smearing width | clear | implemented | `contact_smearing_sigma` is the shared width for contact, central smearing, and smeared derivative kernels. |
| (A10)-(A11) pointwise `G(r)` and `S(r)` | clear | implemented | `static_coulomb_G`, `static_confinement_S`, and `alpha_s_r` follow these forms. |
| (A12)-(A14) closed-form `G~`, `S~`, and `tau_k` | clear | active FD analogue | `smeared_coulomb_G_closed`, `smeared_confinement_S_closed`, and `appendix_a_closed_central_values` implement these formulas; tests compare behavior and derivatives. |
| Coulomb momentum factor after (A14) | clear | active FD analogue | `appendix_a_momentum_sandwich_matrix` builds `A(p) G~ A(p)` on the FD `p^2` eigenbasis. |
| Spin-dependent `m/E` factor after (A14) | clear | active FD analogue | Contact, tensor, vector spin-orbit, and scalar spin-orbit use the two-sided `1/2 + epsilon_i` sandwich when the corresponding switches are enabled. |
| (A15) effective Coulomb-side spin operators | mostly clear | assigned FD analogue | Diagonal fine-structure kernels are present, active kernels use smeared `G~` derivatives, and partnered same-`J` tensor blocks are wired into sector comparison. |
| (A16) scalar/Thomas spin-orbit operator | mostly clear | partial | The equal-mass radial convention is implemented; unequal-mass antisymmetric spin-orbit mixing is now folded into open-flavor same-`J` physical assignment in the FD comparison path. |
| (A17) HO matrix-element factorization | clear | central basis comparison complete | `HarmonicOscillatorBasis` projects radial and momentum operators. The focused central comparison in `docs/residual_reports/appendix_a_ho_comparison.md` finds sub-MeV FD/HO agreement for the active Appendix-A central operator; post-diagonalization tensor/annihilation blocks and HO-order validation of the antisymmetric block remain separate stages. |

## Cleared By This Audit

- The active central path no longer depends on the older derivative proxy as the
  main Appendix-A candidate. The closed-form (A12)-(A14) source is clear enough
  to treat `appendix_a_closed_form` and `appendix_a_momentum_sandwich` as the
  primary spin-independent reproduction paths.
- A5 and A6 should remain derivational context, not direct coding targets.
  The practical implementation target starts at the paper's prescription:
  smearing (A7)-(A14), momentum factors, effective operators (A15)-(A16), and
  HO factorization (A17).

## Next Stages

1. Keep the FD active path as the headline reproduction candidate:
   `G_eff = A(p) G~ A(p) + S~`.
2. The central FD/HO Appendix-A comparison is complete:
   `docs/residual_reports/appendix_a_ho_comparison.md`.
3. Use the HO basis next for paper-order spin-dependent checks, rather than
   reopening A5/A6.
4. Tensor same-`J` mixing is now assigned through
   `assign_mixed_rows(::TensorMixing, ...)`; the unequal-mass antisymmetric
   spin-orbit block is also assigned for open-flavor partner rows.
5. Treat literal isoscalar annihilation/P1/P2 as a separate stage using the Eq. (16)-(18)
   ledger in `docs/formula_map.md`.
