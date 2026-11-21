# Weak-decay research exercise

This directory is a closed exercise on this branch. It computes weak meson
decays from native GI wavefunctions and writes them up in
`report/weak_decay_exploration.pdf` (LaTeX source alongside the PDF). Frozen
numerical results and provenance are in `results/` and `data/`. That report is
the deliverable; the module and scripts exist so the report can be reproduced.

Load the module with `include`. It is not a package, not part of a release,
and not an extension of the QuarkModelTransitions API. A PDF under `output/`
at the repository root is only a local export of this same report.

## Reproduce

From the repository root, instantiate its Julia environment if necessary, then:

```sh
JULIA=/path/to/julia PYTHON=/path/to/python bash WeakDecays/scripts/run_all.sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd WeakDecays/report/weak_decay_exploration.tex
```

Python requires NumPy, SciPy and Matplotlib. The Julia implementation uses only
root-project dependencies and standard libraries. Binary wave caches are excluded
from git and rebuilt by `solve_waves.jl`; CSV outputs are included. The benchmark
is deterministic. `data/runtime.json` records the environment used for this run.

The optional `import_benchmarks.py` and `compare_pdgnn.py` use the PDG.NN Python
environment to refresh snapshots and query its existing model, read-only. They
are **not** part of reproduction from the frozen data and do not retrain PDG.NN.

## Scope and units

- Active-quark S-wave impulse currents, P→P/V form factors and finite-lepton rates.
- Two recoil prescriptions; the scalar spectator boost is intentionally incomplete.
- Coherent S-wave state adapter; unsupported orbital components are rejected.
- Independent light-front P→P plus current at spacelike q² and zero transfer.
- Pole–z continuation experiment; no full BCL unitarity constraints or calibrated errors.
- Leptonic constants, radial excitations, factorized PP modes and finite-width folding.
- A two-pole elastic K-matrix example and a separately calibrated GI K* strong vertex.

All **new research widths are GeV**, unlike QMT's public widths in MeV. Physical
external masses and lifetimes are frozen separately from model wavefunctions.
There is no weak-data fit, CKM extraction, general four-quark calculation, full
Poincaré-covariant current, fitted coupled-channel model, or baryon calculation.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
include("WeakDecays/src/WeakDecayExploration.jl")
using .WeakDecayExploration
params, masses = load_parameters_and_quark_masses(default_parameters_path())
level = BasisState(1, "S", 1, 0)
solver = OscillatorSolver(nbasis=32, converge=false, nlevels_per_channel=1)
i = physical_state(compute_spectrum(params, Meson(masses, :c, :u);
    levels=[level], solver), level)
f = physical_state(compute_spectrum(params, Meson(masses, :s, :u);
    levels=[level], solver), level)
current = current_from_states(f, i, masses; active=:c=>:s)
# Explicit external masses: fixed solved waves, experimental kinematics.
result = integrated_rate(current, 1.86484, 0.493677, 0.000511, 0.9735)
result.width # GeV
```

The report contains both successes and failures. In particular the Galilean
heavy-limit slope violates a QCD sum-rule bound, and the unmodified GI pion/kaon
leptonic currents overshoot. Small numerical errors must not be mistaken for
small theory errors.
