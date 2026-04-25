# Conventions

This document will record the conventions used by the cleaned data and model
implementation.

## Current Code Conventions

- Energies and masses are in GeV internally.
- Distances are in GeV^-1.
- Reference spectrum CSV masses are read in GeV.
- `L` is stored as spectroscopic letters `S`, `P`, `D`, `F`, `G` and mapped to
  orbital angular momentum `0, 1, 2, 3, 4`.
- Multiplicity is `2S+1`; the current contact hyperfine term supports
  singlets (`1`) and triplets (`3`).
- For equal-mass quarkonia, the baseline solver predicts one radial level per
  `(n, L)` before spin shifts.
- The current residual reports compare individual Fig. 6 and Fig. 8 labels
  directly to the baseline prediction plus available spin shifts.

Required topics:

- spectroscopic notation,
- sector naming,
- quark ordering,
- spin operator normalization,
- tensor operator convention,
- spin-orbit decomposition,
- units,
- basis normalization,
- radial wave-function convention,
- mixed-state labeling.
