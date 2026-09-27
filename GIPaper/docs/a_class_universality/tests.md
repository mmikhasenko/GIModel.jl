# Validation record

- `julia --project=GIPaper/scripts QuarkModelTransitions/test/runtests.jl`:
  **703 passed**. Includes 132 linearity/topology/cache checks; all previous
  native Eq. (19), angular/flavor, observable, and Table IV/V tests remain intact.
- `julia --project=. -e 'using Pkg; Pkg.test()'`:
  **1,022 passed**. Pkg's test environment is necessary for test-only dependencies.
  Existing under-resolution/quadrature warnings occurred in deliberate solver
  stress tests; no core source changes were made for this study.
- `julia --project=GIPaper GIPaper/test/runtests.jl`:
  **2,869 passed**, including 217 survey contracts. This suite requires the
  independently maintained `paper/` source archive. The first run in the isolated
  worktree had 20 missing-provenance failures; an ignored symlink to the existing
  local archive supplied those unchanged reference files, after which it passed.

The pilot archives record the original failed **1%** cross-solver gate. The full
survey uses the explicitly documented **3%** gate and exports exact discrepancies
rather than relabeling the original pilot as sub-percent agreement.

The unequal-mass SHO tests certify what the existing operator computes
(spectator fraction), not the requested emitter expression. The discrepancy is
an open physics issue, exposed by `requested_formula_pass` in the data.

The expanded artifact audit passed for **51,840** raw pieces, **9,720** comparisons,
**3,240** certificates and **7,290** plot contributions. Of those contributions,
810 belong to the new standalone threshold histogram. Exact publication-bin
edges are checked against every contribution.

The new A′/A″ tests verify the −1/4 and +1/8 SHO threshold limits for both emitter
lines, and compare finite-q analytic SHO overlaps to independently sampled mesh
waves. At q=0.30 GeV the new light pilot's refined-FD coefficient errors are
0.0054% (A′) and 0.0704% (A″). The original transition-package/core tests above
were run for the unchanged operator implementation in the preceding commit;
the complete GIPaper suite was rerun for this extension.

All **40,320** original raw operator-piece values were compared by trace id to
the preceding archive and are exactly unchanged (maximum relative difference
zero). The expanded threshold census has **798/810** numerical passes: all
180 new-family combinations pass, and only the original twelve flags remain.

The standalone threshold PNG and both validation figures were visually inspected.
All U_hg labels now use LaTeXStrings, including consistently sized h/g subscripts.
The final plotting run completed without warnings.
