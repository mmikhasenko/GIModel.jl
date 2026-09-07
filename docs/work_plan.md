# Remaining Work Plan

This is the repository's only work plan. Per-equation and per-table status lives
in `docs/paper_manifest/*.toml`; completed implementation history belongs in the
audits and Git history, not in parallel plans.

The native-HO spectrum algorithm, physical-state composition, observable wave
interface, and independent FD comparator are complete. One paper-reproduction
unit remains.

## Table VI — complete the photon-decay audit

The canonical transcription contains 79 rows in
`GIPaper/data/raw/digitized_tables/table_vi_photon_decays.csv`. The generated
report currently evaluates 43 transitions. Complete the remaining radial,
excited-state, mixing-induced, and nominally forbidden entries using the
existing Appendix-D overlap and physical-state interfaces.

The implementation should:

1. drive row coverage from the canonical CSV rather than extend another
   hand-maintained subset;
2. evaluate every applicable M1, E1, and M2 row, including signed mixed-state
   amplitudes through the shared physical-state composition;
3. classify rows that cannot be computed from information printed in the paper
   with an explicit reason instead of silently omitting them;
4. resolve or quantitatively characterize cancellation-sensitive residuals,
   especially `Upsilon'' -> eta_b gamma`, without fitting a new constant; and
5. regenerate the report, mark Table VI complete in the manifest, and pass the
   full project verification gate.

Acceptance is complete accounting: all 79 canonical rows must be computed or
explicitly classified, with formula, phase, mass, and mixing conventions stated.

## Completion boundary

Directly applying the Eq. (19) quark-emission operator to the calculated model
wavefunctions is a useful extension beyond the 1985 paper's numerical
single-beta SHO treatment. It is therefore listed in the README as a possible
future improvement, not as unfinished reproduction work.

When Table VI is closed, this file can be removed; ongoing unit status will
remain available from the manifest and generated dashboard.
