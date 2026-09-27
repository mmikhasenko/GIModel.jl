# Research workbench

This directory holds exploratory studies that use GIModel but are not yet part
of the package API, reproduction layer, or polished report.

## Reference history

- **Closed investigation:** [What the GI model says about J=L mixing, 1985 → 2016](gi-later-papers/README.md)
  traces the singlet–triplet mixing angles through Godfrey's later papers and
  an independent reimplementation. From 1991 on, all of them give the small
  angles GIModel reproduces, not the 1985 captions. It ends with a compact
  paragraph for the report.

## Decays and resonances

- [Quark-model strong decays: from width tables to resonance poles](decays/quark_model_resonances.md)
  surveys calculations built from quark-model wave functions and strong-decay
  operators, identifies the limit of isolated-width calculations, and proposes
  a staged GIModel project from exact-wave amplitudes to coupled-channel poles.
- [Strong-decay and resonance TODO](decays/TODO.md) records the implementation
  audit and the staged work, beginning with the meaning and units of the
  existing `matrix_element` API.
- [Transition matrix-element API](decays/transition_matrix_element_api.md)
  proposes the public multiple-dispatch interface, algebra-generation boundary,
  compatibility path, and first vertical implementation slice.
- [Matrix-element API implementation plan](decays/matrix_element_api_implementation_plan.md)
  turns that proposal into reviewable phases, invariants, acceptance gates,
  tests, compatibility rules, alternatives, and a risk register.
- [Review of the matrix-element API plan](decays/matrix_element_api_review.md)
  checks that plan against the code, Appendix B/C, and the `^3P0` literature:
  it argues the helicity amplitude (Eq. C1, Table XI) is the primitive and the
  partial waves its projection, splits the algebra from the spatial integrals
  at the coefficient/orbital-label line, replaces `StateView` with a resolved
  `PhysicalState`, and closes eight of the plan's ten review questions.
