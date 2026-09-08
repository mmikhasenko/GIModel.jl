# GIModel.jl report

`gi_reproduction.qmd` is the source of the project report,
**Recomputing the Godfrey--Isgur Relativized Quark Model**. It follows the
standard shape: introduction and motivation (the model's status, and the
absence of any open constituent quark-model spectroscopy code), computation
and architecture (the Hamiltonian, Appendix-A relativization, the paper's
three-stage algorithm, Table II inputs, and the GIModel/GIPaper and
model/solver separations), implementation (the FD-then-HO route, the
corrections it required, agent-assisted development, provenance and current
status), usage (inspectable answers, the learning track and notebooks,
exploring the model), and conclusion and outlook.

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
