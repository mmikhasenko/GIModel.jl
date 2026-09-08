# GIModel.jl report

`gi_reproduction.qmd` is the source of the project report,
**Recomputing Godfrey--Isgur: An Excursion into a 1985 Quark Model**.
It is a first-person account of the reproduction rather than a description of
a finished tool: what the paper asks you to compute, the FD-then-HO route the
project actually took, the four corrections (normalization and phase, a
converged-but-wrong approximation, the limits of convergence, and a scorecard
reporting the wrong run) that taught the physics, what working with coding
agents was like, what the reproduction currently rests on, and what it left
behind — the package, the learning track, and the notebooks.

The original paper and `LearningTrack/` supply the detailed theory; the report
does not repeat them. The development narrative is anchored to verifiable git
commits and can serve as the basis for a companion website story.

## Building

```bash
quarto render report/gi_reproduction.qmd
```

Requires Quarto and a LaTeX toolchain. The rendered `gi_reproduction.pdf` is
tracked so readers do not need the toolchain. The `.tex` file is generated.

## Evidence and scope

The report describes the current native HO spectrum route and independent FD
comparator. Full-depth paper reproduction is the release objective; the current
audit still lists complete Table VI photon-decay coverage and residual
characterization as unfinished. Keep that boundary synchronized with
`docs/reproduction_audit.md` and `docs/paper_manifest/`.

Quantitative examples identify their checked-in source audit and snapshot.
This editorial revision does not regenerate model results. Before changing
numerical claims, regenerate the relevant audit and inspect its solver settings,
input provenance, and treatment of calibrated controls. The full computation
gate is `bash scripts/verify_project.sh`.

Earlier spectrum and wavefunction images remain in `figures/` as historical
assets; the revised report does not use those older plots as current evidence.
The curated teaching demonstrations are under `examples/`.
