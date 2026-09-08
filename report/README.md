# GIModel.jl report

`gi_reproduction.qmd` is the source of the project report,
**GIModel.jl: An Inspectable Quark Model for Computation and Learning**.
The report explains the value of the software and its companion learning
material: inspectable numerical answers, language-model-assisted exploration,
Pluto and live documentation, flavor-sector visualizations, open-source reuse,
and future physics and differentiation work. The original paper and
`LearningTrack/` supply the detailed theory.

The development appendix follows the FD-to-native-HO route through verifiable
git commits, including normalization, phase, staging, and verification failures.
It can also serve as the narrative basis for a companion website story.

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
