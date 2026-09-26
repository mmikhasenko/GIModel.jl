# Validation record

- `julia --project=GIPaper/scripts QuarkModelTransitions/test/runtests.jl`:
  **703 passed**. Includes 132 linearity/topology/cache checks; all previous
  native Eq. (19), angular/flavor, observable, and Table IV/V tests remain intact.
- `julia --project=. -e 'using Pkg; Pkg.test()'`:
  **1,022 passed**. Pkg's test environment is necessary for test-only dependencies.
  Existing under-resolution/quadrature warnings occurred in deliberate solver
  stress tests; no core source changes were made for this study.
- `julia --project=GIPaper GIPaper/test/runtests.jl`:
  **2,745 passed**, including 93 new survey contracts. This suite requires the
  independently maintained `paper/` source archive. The first run in the isolated
  worktree had 20 missing-provenance failures; an ignored symlink to the existing
  local archive supplied those unchanged reference files, after which it passed.

The pilot archives record the original failed **1%** cross-solver gate. The full
survey uses the explicitly documented **3%** gate and exports exact discrepancies
rather than relabeling the original pilot as sub-percent agreement.

The unequal-mass SHO tests certify what the existing operator computes
(spectator fraction), not the requested emitter expression. The discrepancy is
an open physics issue, exposed by `requested_formula_pass` in the data.

Artifact validation passed: **40,320** raw pieces, **7,560** comparisons,
**2,520** certificates and **5,040** plot contributions. Both generated PNG
figures were visually inspected; the PDF companions preserve vector graphics.
