```@meta
EditURL = "../../quarto/manual/solvers.qmd"
```



# Solvers and convergence

GIModel solves the same Hamiltonian in two independent ways. This page explains both, how to choose between them, and how to check that an answer is converged.

```julia
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
charmonium = Meson(mq, :c, :c)
levels = spectrum_levels(1; L_labels = ("S", "P"))
```

## Two methods, one answer

|  | [`OscillatorSolver`](@ref) | [`FiniteDifferenceSolver`](@ref) |
|----|----|----|
| method | expansion in harmonic-oscillator functions | matrices on a uniform radial mesh |
| used by the paper | yes | no |
| accuracy set by | basis size and oscillator scale $\beta$ (chosen automatically) | mesh spacing and box size (chosen by you) |
| convergence check | built in, with a certificate | your own refinement study |
| wave type | [`OscillatorWave`](@ref) | [`MeshWave`](@ref) |
| default in | [`compute_isoscalar_spectrum`](@ref) | [`compute_spectrum`](@ref) |

Both solvers produce the same kind of result, and everything downstream (mixing, wavefunctions, transitions) works with either. The two are separate implementations of the same operators, so their agreement is a strong check.

```julia
ho = compute_spectrum(params, charmonium; levels, solver = OscillatorSolver())
fd = compute_spectrum(params, charmonium; levels, solver = FiniteDifferenceSolver())
[(s.label, round(1000 * (f.mass_GeV - s.mass_GeV); digits = 3)) for (s, f) in zip(ho.states, fd.states)]
```

```
6-element Vector{Tuple{String, Float64}}:
 ("1^1S_0", -0.389)
 ("1^3S_1", -0.14)
 ("1^1P_1", -0.071)
 ("1^3P_0", -0.118)
 ("1^3P_1", -0.076)
 ("1^3P_2", -0.055)
```

The differences are in MeV. The default mesh is already well below 1 MeV.

## The harmonic-oscillator solver

The paper’s method expands each radial wavefunction as

```math
u(r) = \sum_{k=0}^{N-1} c_k\, R_{kL}(\beta; r),
```

where $R_{kL}$ are harmonic-oscillator radial functions with scale $\beta$. Kinetic and potential matrix elements are computed in closed form or by Gauss–Laguerre quadrature; no radial mesh is involved.

For each $(L, S, J)$ sector the solver:

1.  scans $\beta$ over `beta_grid` and picks the value that minimizes the highest requested level (the paper’s choice), then refines it;
2.  enlarges the basis from `nbasis` in steps of `basis_step` until every requested energy changes by less than `energy_tolerance_GeV` (0.1 MeV by default) on two successive refinements;
3.  fails loudly if the optimal $\beta$ lies at the edge of the grid or if `max_nbasis` is reached without convergence.

```julia
OscillatorSolver()
```

```
OscillatorSolver: nbasis = 24:8:80 adaptive, ΔE ≤ 0.1 MeV, beta in [0.25, 2.35] GeV (22 bracket points, refined to 0.002 GeV), 6 levels/channel
```

The outcome is recorded as an [`OscillatorConvergence`](@ref) certificate on each solved sector. The spectrum keeps its solved sectors in `spec.computation.channel_cache`:

```julia
for (key, solution) in ho.computation.channel_cache
    c = solution.convergence
    println(rpad("$(key.L_label) 2S+1=$(key.multiplicity) J=$(key.J)", 18),
        c.status, "  β = ", round(c.beta_GeV; digits = 3), " GeV  N = ", c.nbasis,
        "  ΔE = ", round(1e3 * c.energy_delta_GeV; sigdigits = 2), " MeV")
end
```

```
P 2S+1=3 J=0      converged  β = 1.094 GeV  N = 40  ΔE = 3.5e-5 MeV
P 2S+1=3 J=1      converged  β = 1.056 GeV  N = 40  ΔE = 2.9e-5 MeV
S 2S+1=1 J=0      converged  β = 1.353 GeV  N = 40  ΔE = 8.9e-5 MeV
S 2S+1=3 J=1      converged  β = 1.163 GeV  N = 40  ΔE = 2.9e-6 MeV
P 2S+1=1 J=1      converged  β = 1.102 GeV  N = 40  ΔE = 0.00013 MeV
P 2S+1=3 J=2      converged  β = 1.0 GeV  N = 40  ΔE = 3.9e-6 MeV
```

Each sector gets its own $\beta$, because each has a different size.

The certificate covers the **energies** that were requested. Other quantities, such as a wavefunction at the origin or a transition overlap, can converge more slowly. When such a quantity matters, check it directly, for example by tightening `energy_tolerance_GeV` or raising `nbasis`, and compare.

`OscillatorSolver(converge = false)` performs a single fixed-size solve without the check. Its certificate says `:unchecked`. This is meant for convergence studies, not for results.

## The finite-difference solver

The finite-difference solver represents the wavefunction on `ngrid` interior points $r_i = i\,r_\text{max}/(n_\text{grid}+1)$, with $u = 0$ at both ends. The kinetic operator $\sqrt{p^2 + m^2}$ is built as a matrix function of the discretized $p^2$.

```julia
FiniteDifferenceSolver()
```

```
FiniteDifferenceSolver: ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel
```

Its accuracy depends on two numbers you choose:

- **spacing** $h = r_\text{max}/(n_\text{grid}+1)$ must resolve the wave near the origin, which matters most for compact heavy quarkonia;
- **box size** $r_\text{max}$ (in $\mathrm{GeV}^{-1}$) must contain the wave, which matters most for light and radially excited states.

There is no automatic check, so refine and compare. The error falls roughly as $h^2$:

```julia
for ngrid in (450, 900, 1800)
    fine = compute_spectrum(params, charmonium; levels,
        solver = FiniteDifferenceSolver(ngrid = ngrid))
    shift = maximum(abs(f.mass_GeV - h.mass_GeV) for (f, h) in zip(fine.states, ho.states))
    println("ngrid = ", ngrid, ":  max |FD − HO| = ", round(1e3 * shift; digits = 3), " MeV")
end
```

```
ngrid = 450:  max |FD − HO| = 0.389 MeV
ngrid = 900:  max |FD − HO| = 0.097 MeV
ngrid = 1800:  max |FD − HO| = 0.024 MeV
```

To change one setting and keep the rest, pass an existing solver as the first argument:

```julia
FiniteDifferenceSolver(FiniteDifferenceSolver(); ngrid = 900)
```

```
FiniteDifferenceSolver: ngrid = 900, rmax = 24.0 GeV^-1 (h = 0.02664), relativistic, full, 6 levels/channel
```

The finite-difference solver has two extra options: `kinetic = :nonrelativistic` replaces $\sqrt{p^2 + m^2}$ by $p^2/2\mu$ for comparison, and `eigensolver = :krylov` uses an iterative eigensolver for large meshes.

## Which one to use

- For **results to quote**, use `OscillatorSolver()`. It is the paper’s method and certifies its own energies.
- For **quick exploration**, the default `FiniteDifferenceSolver()` is fast and accurate to a fraction of an MeV for most states.
- For a **cross-check**, compute with both. A difference larger than the tolerances points to a convergence problem, not to physics.
- For **wavefunction-sensitive quantities** (widths, radii, values at the origin), check convergence of that quantity itself with either solver.

## Solver settings are recorded

Reports and papers should state how numbers were computed. [`numerics_provenance`](@ref) turns solver settings into one line:

```julia
println(numerics_provenance(OscillatorSolver()))
println(numerics_provenance(FiniteDifferenceSolver()))
```

```
Numerics: `OscillatorSolver` — nbasis = 24:8:80 adaptive, ΔE ≤ 0.1 MeV, beta in [0.25, 2.35] GeV (22 bracket points, refined to 0.002 GeV), 6 levels/channel.
Numerics: `FiniteDifferenceSolver` — ngrid = 450, rmax = 24.0 GeV^-1 (h = 0.05322), relativistic, full, 6 levels/channel.
```

## Solving a single channel

The spectrum functions call [`channel_solution`](@ref) for each sector. You can call it directly to get the spin-independent radial eigenproblem of one orbital channel:

```julia
sol = channel_solution(params, charmonium.constituent_masses, 0;
    solver = OscillatorSolver())
sol.eigenvalues_GeV[1:3]
```

```
3-element Vector{Float64}:
 3.0645237923656308
 3.666192238050535
 4.0909574961787625
```

The result is a [`ChannelRadialSolution`](@ref): eigenvalues, one wave per level (`radial_wave(sol, n)`), and for the oscillator solver a convergence certificate. [`fixed_channel_solution`](@ref) does the same for a complete fixed $(L, S, J)$ Hamiltonian including spin terms.
