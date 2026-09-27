# Learning track

The repository includes a pen-and-paper course, in `LearningTrack/`, that
leads from elementary quantum mechanics to the complete Godfrey–Isgur
algorithm. It is written for physics students who want to understand the model
well enough to implement it.

## Contents

**Part 1: theory sheets**

1. Radial quantum mechanics
2. Angular momentum and mesons
3. Variational methods and the semirelativistic kinetic energy
4. Color, confinement and the running coupling
5. Relativization and smearing
6. Spin-dependent fine structure
7. State and flavor mixing
8. Decay observables
9. The full algorithm: a capstone

**Part 2: computational track**

- C1: computational discovery

Each sheet contains problems, worked solutions and short concept checks after
the solutions. Problem statements give all inputs, conventions and the
requested result, so a student can start without reading the solution.

## Building the PDFs

A LaTeX installation with `latexmk` is required:

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
