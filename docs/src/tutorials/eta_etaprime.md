```@meta
EditURL = "../../quarto/tutorials/eta_etaprime.qmd"
```



# η and η′

Without annihilation, the quark model would put an $n\bar n$ pseudoscalar degenerate with the pion and an $s\bar s$ pseudoscalar near 690 MeV. The physical $\eta(548)$ and $\eta'(958)$ look nothing like that. This tutorial follows the Godfrey–Isgur treatment of annihilation mixing in the pseudoscalar sector: how the masses move, what the states are made of, and why the answer depends on the radial basis.

```julia
using GIModel
using GIModel.QuarkModelTransitions
params, mq = load_parameters_and_quark_masses(default_parameters_path())
nn, ss = Meson(mq, :q, :q), Meson(mq, :s, :s)
pseudoscalars = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
solver = OscillatorSolver()
```

## Before annihilation

Solve the two flavor channels separately:

```julia
before = [spectrum_state(compute_spectrum(params, m; levels = pseudoscalars, solver), n, "S", 1, 0).mass_GeV
          for m in (nn, ss), n in 1:2]
```

```
2×2 Matrix{Float64}:
 0.148877  1.29131
 0.655116  1.57724
```

Rows are $n\bar n$ and $s\bar s$, columns are $1S$ and $2S$. The $n\bar n$ ground state is the pion, far below the $\eta$.

## Two models from the paper

The paper gives two forms of the annihilation amplitude for pseudoscalars, Eq. (18a) (`P1`) and Eq. (18b) (`P2`). Both add a nonperturbative term to the perturbative two-gluon amplitude; `P2` makes the nonperturbative term depend on the meson mass.

```julia
for model in (PaperP1Annihilation(), PaperP2Annihilation())
    iso = compute_isoscalar_spectrum(params, nn, ss; levels = pseudoscalars,
        pseudoscalar = model, solver)
    masses = sort([s.mass_GeV for s in iso.states])
    println(rpad(nameof(typeof(model)), 22), round.(masses; digits = 3))
end
```

```
PaperP1Annihilation   [0.531, 1.0, 1.496, 1.799]
PaperP2Annihilation   [0.49, 0.974, 1.262, 1.55]
```

The four eigenvalues correspond to $\eta$, $\eta'$ and their radial excitations. Annihilation lifts the pion-like $n\bar n$ state to the $\eta$ region and pushes the $s\bar s$-like state up toward the $\eta'$.

## What the states are made of

```julia
iso = compute_isoscalar_spectrum(params, nn, ss; levels = pseudoscalars,
    pseudoscalar = PaperP1Annihilation(), solver)
eta = spectrum_state(iso, BasisState(1, "S", 1, 0; flavors = (:q, :q)))
eta_prime = spectrum_state(iso, BasisState(1, "S", 1, 0; flavors = (:s, :s)))
composition(state) = [(join(c.basis.flavors), c.basis.n, round(c.coefficient; digits = 3))
                      for c in physical_components(iso, state)]
composition(eta)
```

```
4-element Vector{Tuple{String, Int64, Float64}}:
 ("qq", 1, 0.625)
 ("ss", 1, -0.772)
 ("qq", 2, 0.111)
 ("ss", 2, 0.041)
```

```julia
composition(eta_prime)
```

```
4-element Vector{Tuple{String, Int64, Float64}}:
 ("qq", 1, 0.568)
 ("ss", 1, 0.552)
 ("qq", 2, 0.591)
 ("ss", 2, 0.152)
```

In the flavor basis, SU(3) symmetry would give an octet $\eta_8 \propto n\bar n - \sqrt2\,s\bar s$ and a singlet $\eta_1 \propto \sqrt2\,n\bar n + s\bar s$. The model $\eta$ has opposite-sign $n\bar n$ and $s\bar s$ components and is octet-like; the $\eta'$ has same-sign components and is singlet-like. The $\eta'$ also contains a large $2S$ component: annihilation connects all states with the same quantum numbers, not only ground states.

## The radial basis matters

Repeat the calculation with only ground states in the block:

```julia
ground = [compute_spectrum(params, m; levels = pseudoscalars[1:1], solver) for m in (nn, ss)]
small = add_isoscalar_annihilation(params, ground;
    pseudoscalar = PaperP1Annihilation(),
    pseudoscalar_basis = [BasisState(1, "S", 1, 0; flavors = (f, f)) for f in (:q, :s)])
sort([s.mass_GeV for s in small.states])
```

```
2-element Vector{Float64}:
 0.5417374618268833
 1.2772003433559722
```

The $\eta$ barely moves, but the $\eta'$ rises by roughly 300 MeV. In the larger basis, mixing with the radially excited states pulls it down. The annihilation amplitude involves the wavefunction at the origin, which is large for $2S$ states, so these couplings are not small. The paper uses the four-state basis; [`compute_isoscalar_spectrum`](@ref) always builds that four-state pseudoscalar block, and the general [`add_isoscalar_annihilation`](@ref) lets you choose another.

## Two-photon widths

The two-photon amplitude of each flavor component is proportional to its squared charge. For the normalized $(u\bar u + d\bar d)/\sqrt2$ component the effective factor is $(e_u^2 + e_d^2)/\sqrt2$, and for $s\bar s$ it is $e_s^2$:

```julia
charges = [AnnihilationTerm((:q, :q), (4 / 9 + 1 / 9) / sqrt(2)),
           AnnihilationTerm((:s, :s), 1 / 9)]
two_photon = TwoPhotonAnnihilation(mq, charges)
for (name, state) in (("η", eta), ("η′", eta_prime))
    ps = physical_state(iso, state)
    println(name, ":  M = ", round(ps.mass_GeV; digits = 3), " GeV,  Γ(γγ) = ",
        round(1e6 * decay_width(TwoPhotonChannel(), two_photon, ps); digits = 2), " eV")
end
```

```
η:  M = 0.531 GeV,  Γ(γγ) = 173.85 eV
η′:  M = 1.0 GeV,  Γ(γγ) = 1182.51 eV
```

The flavor components interfere in the amplitude: the $\eta$ has opposite-sign components, so its $n\bar n$ and $s\bar s$ contributions partly cancel. The width also scales as $M^3$, so the model masses matter. [`mass_correction_factor`](@ref) shows the effect of evaluating at the measured mass, with the wavefunctions held fixed:

```julia
ps = physical_state(iso, eta)
factor = mass_correction_factor(TwoPhotonChannel(), two_photon, ps; target_mass = 0.548)
amplitude = matrix_element(TwoPhotonChannel(), two_photon, ps)
1e6 * 1000 * abs2(factor * amplitude.value)       # eV at M = 548 MeV
```

```
190.5679874027079
```

The paper’s own comparison of these widths is in the [Table VII report](../paper/reports/table_vii_annihilation_em.md).

## Vectors and tensors

For ${}^3S_1$ and ${}^3P_2$ states the paper uses the general annihilation amplitude of Eq. (16), with $A({}^3S_1) = 2.5$ and $A({}^3P_2) = -0.8$ from Table III. The mixing is much weaker because it requires three and two gluons respectively, and the $\omega$–$\phi$ and $f_2$–$f_2'$ pairs stay close to ideal mixing:

```julia
A = params.annihilation
vt = compute_isoscalar_spectrum(params, nn, ss;
    levels = [BasisState(1, "S", 3, 1), BasisState(1, "P", 3, 2)], solver,
    amplitudes = Dict(("S", 3, 1) => A.s1_A, ("P", 3, 2) => A.a_3p2))
```

```
MixedSpectrum: qq + ss, 4 levels — all values in GeV
  level        central    contact   fine str     mixing       mass
  1^3S_1 qq     0.6654     0.1035     0.0000     0.0099     0.7789
  1^3P_2 qq     1.2679     0.0086     0.0301    -0.0307     1.2758
  1^3S_1 ss     0.9482     0.0679     0.0000     0.0033     1.0194
  1^3P_2 ss     1.5040     0.0061     0.0295    -0.0078     1.5316
  (4 levels carry mixing; see `spec.states[i].mixings`)
```

```julia
omega = spectrum_state(vt, BasisState(1, "S", 3, 1; flavors = (:q, :q)))
[(join(c.basis.flavors), round(c.coefficient; digits = 3)) for c in physical_components(vt, omega)]
```

```
2-element Vector{Tuple{String, Float64}}:
 ("qq", 1.0)
 ("ss", -0.023)
```
