# Learning track

The repository includes a pen-and-paper course, in `LearningTrack/`, that
leads from elementary quantum mechanics to the complete Godfrey–Isgur
algorithm. It is written for physics students who want to understand the model
well enough to implement it.

## Contents

**Part 1: theory sheets**

1. [Radial quantum mechanics](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/01_radial_quantum_mechanics.pdf)
2. [Angular momentum and mesons](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/02_angular_momentum_and_mesons.pdf)
3. [Variational methods and the semirelativistic kinetic energy](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/03_variational_and_semirelativistic.pdf)
4. [Color, confinement and the running coupling](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/04_color_confinement_running_coupling.pdf)
5. [Relativization and smearing](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/05_relativization_and_smearing.pdf)
6. [Spin-dependent fine structure](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/06_spin_fine_structure.pdf)
7. [State and flavor mixing](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/07_state_and_flavor_mixing.pdf)
8. [Decay observables](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/08_decay_observables.pdf)
9. [The full algorithm: a capstone](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/09_full_algorithm_capstone.pdf)

**Part 2: computational track**

- C1: [computational discovery](https://github.com/mmikhasenko/GIModel.jl/releases/download/learning-sheets/C1_computational_discovery.pdf)

Each sheet contains problems, worked solutions and short concept checks after
the solutions. Problem statements give all inputs, conventions and the
requested result, so a student can start without reading the solution.

The links download the PDFs from the latest successful build on `main`, which
the Learning sheets workflow publishes to a dedicated
[release](https://github.com/mmikhasenko/GIModel.jl/releases/tag/learning-sheets).

## Building the PDFs

To build them yourself, a LaTeX installation with `latexmk` is required:

```bash
cd LearningTrack
make
```

The PDFs are written to `LearningTrack/pdf/`.

## Using it with the package

Each sheet corresponds to a part of this documentation:

| sheet | documentation |
|---|---|
| 1–3 | [The model](@ref), [Solvers and convergence](@ref) |
| 4–6 | [The model](@ref), [Computing a spectrum](@ref) |
| 7 | [Computing a spectrum](@ref), [Isoscalar flavor mixing](@ref) |
| 8 | [Transitions and decays](@ref) |
| 9 | [Charmonium, start to finish](@ref) |
