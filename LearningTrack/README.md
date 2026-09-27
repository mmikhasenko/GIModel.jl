# Agentic Godfrey-Isgur

A learning track from elementary quantum mechanics to the Godfrey-Isgur model.
`part1/` contains nine pen-and-paper theory sheets; `part2/` begins the
computational track. Each sheet includes problems, worked solutions, and short
post-solution concept checks.

Problem statements specify the inputs, conventions, and requested result.
Necessary assumptions belong in the statement; short `Hint:` paragraphs supply
an identity or a starting method. Students should be able to begin from the
statement and its cited earlier results without reading the worked solution.
Conceptual design questions specify the expected output, such as a table,
controlled comparison, or brief pseudocode.

## Download the sheets

Each link downloads the PDF from the latest successful learning-sheet build on
`main`. The dedicated [Learning sheets release](https://github.com/mmikhasenko/GIModel.jl/releases/tag/learning-sheets)
records the source commit; it is separate from package releases.

| Sheet | PDF |
|---|---|
| Radial quantum mechanics | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/01_radial_quantum_mechanics.pdf) |
| Angular momentum and mesons | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/02_angular_momentum_and_mesons.pdf) |
| Variational and semirelativistic methods | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/03_variational_and_semirelativistic.pdf) |
| Color, confinement and running coupling | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/04_color_confinement_running_coupling.pdf) |
| Relativization and smearing | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/05_relativization_and_smearing.pdf) |
| Spin fine structure | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/06_spin_fine_structure.pdf) |
| State and flavor mixing | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/07_state_and_flavor_mixing.pdf) |
| Decay observables | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/08_decay_observables.pdf) |
| Full algorithm capstone | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/09_full_algorithm_capstone.pdf) |
| Computational discovery map | [Download](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/C1_computational_discovery.pdf) |

The [Learning sheets workflow](https://github.com/mmikhasenko/GIModel.jl/actions/workflows/LearningSheets.yml)
builds all sheets on relevant pushes and pull requests, with a cached TeX Live
installation. Pull-request builds retain a PDF bundle for 30 days; successful
`main` builds also refresh the individual downloads above and list them in the
workflow summary. A manual run on `main` can refresh the published PDFs.

## Compile

A LaTeX installation with `latexmk` is required. From this directory, run:

```sh
make
```

The generated PDFs are written to the gitignored `pdf/` directory. To remove
LaTeX intermediate files while keeping the PDFs, run:

```sh
make clean
```
