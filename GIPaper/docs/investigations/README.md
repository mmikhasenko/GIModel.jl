# GIPaper investigations

This directory contains reproducible, hypothesis-driven studies of mechanisms
inside the Godfrey–Isgur calculation.

It is deliberately separate from `../residual_reports/`:

- `residual_reports/` is the authoritative paper-reproduction ledger. It asks
  which published values differ and records every comparison row.
- `investigations/` asks why a localized discrepancy occurs. A study may
  separate operators, compare numerical realizations, or test a physical
  hypothesis, but it must link back to the canonical residual report rather
  than copy its paper data or implement a second comparison path.

## Standards for an investigation

Every numerical study here should:

1. name its executable generator and numerical settings;
2. use GIModel's public computation API wherever possible;
3. distinguish coordinate-kernel plots from complete momentum-sandwiched
   matrix elements;
4. include an independent numerical control when the conclusion depends on a
   solver;
5. state what was ruled out, what remains unresolved, and the next calculation
   that could discriminate between explanations;
6. keep generated figures under `figures/` and avoid embedding paper reference
   values already owned by a residual report.

## Studies

- [Report refresh after the contact correction](report_contact_refresh.md) —
  dependency map, regeneration commands, and validation for `report/main.tex`.

- **Closed; historical brief:** [Mixing matching targets](mixing_matching_targets.md) — distinguishes direct
  caption angles and tensor compositions from conditional energy-shift
  reconstructions, and checks the two-stage algorithm and its basis limits.
- **Closed:** [Why the mixing compositions failed](mixing_composition_investigation.md) —
  demonstrates the missing L>0 smeared contact term as the cause of the wrong
  singlet–triplet diagonal gaps (fixed in core), certifies both mixing stages,
  independently re-derives the antisymmetric element, and shows the remaining
  ground-state off-diagonal difference is in the 1985 captions: Godfrey's later
  GI-model papers agree with GIModel. Numbers: [results](mixing_composition_results.md).
- [Quantitative mixing-layer review](mixing_layer_review.md) (numbers predate
  the contact fix) — inventories the
  ten-panel spectrum, reports all mixing-only shifts and projected angles,
  checks production eigensystems, and distinguishes paper-angle agreement
  from conditional shift estimates.
- [Discrepancies beyond orbital P waves](non_pwave_discrepancies.md) — tests
  excited-eta composition and tensor two-photon flavor substitutions, checks
  printed P1-vector consistency, and distinguishes demonstrated explanations
  from remaining hypotheses across the other residual families.
- [Light and charmed P-wave fine structure](pwave_fine_structure.md) — separates
  vector spin-orbit, scalar/Thomas, and tensor kernels; compares central-wave
  first order with full fixed-sector diagonalization; and checks native HO
  against finite differences.
- [Differentiable state-energy readiness](differentiable_state_energy.md) —
  establishes that an isolated fixed-sector energy should be differentiated as
  an eigenvalue problem, verifies its Hellmann–Feynman pullback numerically, and
  identifies the adaptive solver and `Float64` result boundaries that must stay
  outside the reverse-mode kernel.

Regenerate the P-wave study from the repository root with:

```sh
julia GIPaper/scripts/investigate_pwave_fine_structure.jl
```

Reproduce the state-energy pullback check with:

```sh
julia GIPaper/scripts/investigate_state_energy_differentiability.jl
```

Regenerate the mixing-composition study (the `before` ledger needs a checkout of 88e9a99):

```sh
julia --project=GIPaper/scripts GIPaper/scripts/investigate_mixing_composition.jl --label after
julia --project=GIPaper/scripts GIPaper/scripts/investigate_mixing_composition.jl --report
```

Replay the non-P-wave composition diagnostics from the recorded native inputs:

```sh
julia --project=GIPaper/scripts GIPaper/scripts/investigate_excited_eta_moments.jl
julia --project=GIPaper/scripts GIPaper/scripts/investigate_tensor_photons.jl
```

These require the files produced by `trace_rate_inputs.jl` in `../input_traces/`.
They do not modify the canonical residual reports or refit model parameters.
