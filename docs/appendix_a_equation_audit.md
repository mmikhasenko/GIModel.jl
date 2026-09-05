# Appendix A Equation Audit

> **Detail layer.** Per-equation status/code/page for A1-A17 is the manifest
> (`docs/paper_manifest/appendix_a.toml`, rendered in the dashboard). This file
> is the deeper **source-and-derivation audit** those units point back to. Its
> historical page references (page-036 for A1-A9, page-037 for A10-A15, page-038
> for A16-A17) are the OCR-verified anchors the manifest now uses.

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
| (A12)-(A14) closed-form `G~`, `S~`, and `tau_k` | clear | native HO and FD implementations | `smeared_coulomb_G_closed`, `smeared_confinement_S_closed`, and `appendix_a_closed_central_values` supply the radial functions. FD samples them on its native grid; HO integrates them directly in oscillator matrix elements. |
| Coulomb momentum factor after (A14) | clear | native HO and FD implementations | FD uses `appendix_a_momentum_sandwich_matrix` on its `p²` eigenbasis; HO uses `ho_momentum_sandwich_matrix` with the exact oscillator `p²` matrix. |
| Spin-dependent `m/E` factor after (A14) | clear | native HO and FD implementations | Contact, tensor, vector spin-orbit, and scalar spin-orbit use the same two-sided `1/2 + epsilon_i` sandwich through backend-native matrices. |
| (A15) effective Coulomb-side spin operators | mostly clear | implemented in both backends | `fine_structure_grid_matrices` and `ho_fine_structure_matrices` assemble the corresponding native fixed-sector blocks; cross-sector tensor elements dispatch on native waves. |
| (A16) scalar/Thomas spin-orbit operator | mostly clear | partial formula certification | The native fixed-sector symmetric term and unequal-mass antisymmetric same-`J` mixing are implemented for both representations. The fitted `k_spin_orbit` bridge still prevents a literal paper-strength certification. |
| (A17) HO matrix-element factorization | clear | complete fixed-sector implementation | Both sides of the factorization use exact oscillator matrix elements: `ho_p2_matrix` for momentum functions, `ho_operator_matrix` (generalized Gauss-Laguerre in Golub-Welsch form) for position functions, and spectral functions of exact `p²` for momentum sandwiches. `fixed_channel_solution` assembles central, contact, symmetric spin-orbit, and diagonal tensor matrices before diagonalization; later tensor, antisymmetric spin-orbit, and annihilation blocks consume the resulting signed `OscillatorWave`s directly. The adaptive controller refines beta and basis size and attaches an `OscillatorConvergence` certificate; no spatial mesh appears in this path. |

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

1. Keep FD as an independent implementation comparator, not the headline paper
   route or a hidden HO dependency.
2. Migrate the remaining flavor-sensitive observables to the final
   `physical_components` composition (PA-17).
3. Remove or explicitly isolate the fitted spin bridge factors and certify the
   end-to-end native-HO paper route (PA-18).
