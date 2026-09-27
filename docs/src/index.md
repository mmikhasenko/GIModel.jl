# GIModel.jl

GIModel is a Julia implementation of the Godfrey–Isgur relativized quark
model of mesons: S. Godfrey and N. Isgur, *Phys. Rev. D* **32**, 189 (1985),
[doi:10.1103/PhysRevD.32.189](https://doi.org/10.1103/PhysRevD.32.189).

Given two quark flavors, it builds the relativized quark–antiquark Hamiltonian,
solves it, and returns the meson spectrum. Each state comes with its mass, the
contributions that make up that mass, its mixing with nearby states, and its
wavefunction. The submodule `GIModel.QuarkModelTransitions` uses those
wavefunctions to compute radiative, strong, leptonic and annihilation decays.

```@example index
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
spectrum = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(1))
```

Each row is one ``n\,{}^{2S+1}L_J`` level of charmonium. The columns split its
mass into the spin-independent part and the spin-dependent corrections; the
[Getting started](@ref) page explains the whole table.

## What the package does

- **Spectra** for any ``q_1\bar q_2`` pair of `u`, `d`, `s`, `c` and `b`,
  with all levels labeled by explicit quantum numbers.
- **The paper's algorithm**: a harmonic-oscillator expansion with a
  variationally chosen scale and automatic convergence, plus an independent
  finite-difference solver to cross-check it.
- **Mixing** between states, following the paper's stages: same-``J``
  spin-orbit and tensor mixing inside a meson, and annihilation mixing between
  ``n\bar n``, ``s\bar s``, ``c\bar c`` and ``b\bar b`` isoscalars.
- **Wavefunctions** in position and momentum space for every state, including
  the signed components of mixed states.
- **Transitions** computed from those wavefunctions: E1/M1/M2 photon emission,
  pseudoscalar emission, leptonic and two-photon widths, and gluonic widths.
- **A comparison with the 1985 paper** kept in a separate package, GIPaper,
  which records how each table is reproduced and where it differs.

## How to read these pages

The documentation is ordered so that each page builds on the previous ones.

1. [Getting started](@ref) installs the package and walks one calculation
   through from inputs to wavefunctions.
2. The **Manual** explains the concepts, one topic per page: the physics of
   the model, its inputs, the spectrum, the solvers, wavefunctions, flavor
   mixing, transitions and conventions.
3. The **Tutorials** are complete worked studies that combine those pieces:
   charmonium from start to finish, heavy-light mixing, the charm-to-bottom
   limit, η–η′ mixing, and strong decays.
4. **The 1985 paper** section describes how the package reproduces the
   original article, what agrees, and what differs.
5. The **API reference** lists every exported name with its docstring.

If you know quark models and want the code, read [Getting started](@ref) and
then [Computing a spectrum](@ref). If you are new to the Godfrey–Isgur model,
read [The model](@ref) first; the [Learning track](@ref) is a pen-and-paper
course on the same material.

## Citing

If you use GIModel in published work, cite the original article: S. Godfrey
and N. Isgur, *Phys. Rev. D* **32**, 189 (1985). Please also cite this
repository.
