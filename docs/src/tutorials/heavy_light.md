```@meta
EditURL = "../../quarto/tutorials/heavy_light.qmd"
```



# Heavy-light mesons and mixing

When the quark and antiquark have different masses, the spin-orbit force has a part that does not conserve total quark spin. It mixes the ${}^1P_1$ and ${}^3P_1$ states into the two physical axial mesons, for example $K_1(1270)$/$K_1(1400)$ and $D_1(2420)$/$D_1(2430)$. This tutorial computes that mixing across flavors, takes the off-diagonal matrix element apart, and compares with the published angles.

```julia
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
P_levels = spectrum_levels(1; L_labels = ("P",))
```

## One system: the charmed-strange mesons

```julia
Ds = compute_spectrum(params, Meson(mq, :c, :s); levels = P_levels, solver = OscillatorSolver())
```

```
MixedSpectrum: cs, 4 levels — all values in GeV
  level     central    contact   fine str     mixing       mass
  1^1P_1     2.5650    -0.0147     0.0000    -0.0029     2.5474
  1^3P_0     2.5756     0.0062    -0.1011     0.0000     2.4807
  1^3P_1     2.5655     0.0050    -0.0189     0.0029     2.5545
  1^3P_2     2.5672     0.0037     0.0202     0.0000     2.5911
  (2 levels carry mixing; see `spec.states[i].mixings`)
```

The two $J = 1$ states move apart in the `mixing` column; ${}^3P_0$ and ${}^3P_2$ have no partner with the same $J$ and do not mix. The lower state is a superposition:

```julia
lower = spectrum_state(Ds, "1^1P_1")
[(c.basis.label, round(c.coefficient; digits = 3)) for c in physical_components(Ds, lower)]
```

```
2-element Vector{Tuple{String, Float64}}:
 ("1^1P_1", 0.771)
 ("1^3P_1", -0.637)
```

With the convention of [Conventions and units](@ref), the lower state is $\cos\theta\,|{}^1P_1\rangle + \sin\theta\,|{}^3P_1\rangle$. A small helper reads the angle and the singlet probability:

```julia
function lower_state_angle(spec)
    a, b = spectrum_state(spec, "1^1P_1"), spectrum_state(spec, "1^3P_1")
    low = a.mass_GeV < b.mass_GeV ? a : b
    coefficient(label) = sum((c.coefficient for c in physical_components(spec, low)
                              if c.basis.label == label); init = 0.0)
    s, t = coefficient("1^1P_1"), coefficient("1^3P_1")
    s < 0 && ((s, t) = (-s, -t))          # singlet component positive
    return (θ_deg = round(atan(t, s) * 180 / π; digits = 1), P_singlet = round(s^2; digits = 3))
end
lower_state_angle(Ds)
```

```
(θ_deg = -39.6, P_singlet = 0.594)
```

## All flavors at once

The quark comes first and the antiquark second, following the paper’s sector labels:

```julia
systems = [(:u, :s), (:c, :u), (:c, :s), (:b, :u), (:b, :s), (:b, :c), (:c, :c)]
results = map(systems) do (f1, f2)
    spec = compute_spectrum(params, Meson(mq, f1, f2); levels = P_levels, solver = OscillatorSolver())
    (f1, f2) => lower_state_angle(spec)
end
for ((f1, f2), r) in results
    println(rpad("$f1 $(f2)bar", 10), "θ = ", lpad(r.θ_deg, 6), "°   P(¹P₁) = ", r.P_singlet)
end
```

```
u sbar    θ =    4.3°   P(¹P₁) = 0.994
c ubar    θ =  -25.8°   P(¹P₁) = 0.811
c sbar    θ =  -39.6°   P(¹P₁) = 0.594
b ubar    θ =  -28.0°   P(¹P₁) = 0.78
b sbar    θ =  -41.1°   P(¹P₁) = 0.567
b cbar    θ =   69.2°   P(¹P₁) = 0.127
c cbar    θ =   90.0°   P(¹P₁) = 0.0
```

For $c\bar c$ nothing mixes: with equal masses the antisymmetric operator vanishes. The helper reports 90° only because the lower $J = 1$ state happens to be the pure ${}^3P_1$. The sign changes between $u\bar s$, where the antiquark is the heavier constituent, and $c\bar u$, where the quark is.

## Anatomy of the mixing element

The off-diagonal element has two pieces with opposite signs: one-gluon exchange (vector) and Thomas precession in the confining potential (scalar). [`spin_orbit_mixing_components`](@ref) evaluates them between the unmixed waves:

```julia
function mixing_anatomy(f1, f2)
    meson = Meson(mq, f1, f2)
    fixed = fixed_spectrum(params, meson; levels = P_levels, solver = OscillatorSolver())
    ws = radial_wave(fixed, spectrum_state(fixed, "1^1P_1"))
    wt = radial_wave(fixed, spectrum_state(fixed, "1^3P_1"))
    c = spin_orbit_mixing_components(params, meson.constituent_masses, "P", ws, wt)
    gap = spectrum_state(fixed, "1^3P_1").mass_GeV - spectrum_state(fixed, "1^1P_1").mass_GeV
    return (vector_MeV = round(1e3c.vector; digits = 1), thomas_MeV = round(1e3c.thomas; digits = 1),
            total_MeV = round(1e3c.total; digits = 1), diagonal_gap_MeV = round(1e3gap; digits = 1))
end
for (f1, f2) in systems[1:6]
    println(rpad("$f1 $(f2)bar", 10), mixing_anatomy(f1, f2))
end
```

```
u sbar    (vector_MeV = 7.9, thomas_MeV = -8.9, total_MeV = -1.1, diagonal_gap_MeV = 13.9)
c ubar    (vector_MeV = -25.3, thomas_MeV = 29.1, total_MeV = 3.7, diagonal_gap_MeV = 5.9)
c sbar    (vector_MeV = -18.0, thomas_MeV = 21.5, total_MeV = 3.5, diagonal_gap_MeV = 1.3)
b ubar    (vector_MeV = -31.7, thomas_MeV = 33.4, total_MeV = 1.7, diagonal_gap_MeV = 2.3)
b sbar    (vector_MeV = -24.9, thomas_MeV = 26.6, total_MeV = 1.6, diagonal_gap_MeV = 0.4)
b cbar    (vector_MeV = -9.4, thomas_MeV = 6.8, total_MeV = -2.6, diagonal_gap_MeV = -5.9)
```

The mixing angle is set by the ratio of the off-diagonal element to the gap between the unmixed ${}^1P_1$ and ${}^3P_1$ masses. Both are small differences: the element is a partial cancellation between the vector and Thomas terms, and the gap is a few MeV. The angle is therefore sensitive to details of the spin-dependent potentials, far more than the masses are.

## Comparison with the literature

The 1985 paper quotes the ground-state angles in the captions of Figs. 4, 7 and 9. The same authors’ model, with the same parameters, was used again in later papers; Godfrey and Kokoski (1991), Table I, lists angles for the same six systems. The table below puts the three side by side (degrees, in the convention above; the Godfrey–Kokoski $s\bar u$ value is converted to $u\bar s$):

```julia
published = Dict(
    (:u, :s) => (+34, +5), (:c, :u) => (-41, -26), (:c, :s) => (-44, -38),
    (:b, :u) => (-43, -31), (:b, :s) => (-45, -40), (:b, :c) => (-53, +68))
println(rpad("system", 10), lpad("GIModel", 9), lpad("GI 1985", 9), lpad("GK 1991", 9))
for ((f1, f2), r) in results[1:6]
    gi85, gk91 = published[(f1, f2)]
    println(rpad("$f1 $(f2)bar", 10), lpad(r.θ_deg, 9), lpad(gi85, 9), lpad(gk91, 9))
end
```

```
system      GIModel  GI 1985  GK 1991
u sbar          4.3       34        5
c ubar        -25.8      -41      -26
c sbar        -39.6      -44      -38
b ubar        -28.0      -43      -31
b sbar        -41.1      -45      -40
b cbar         69.2      -53       68
```

GIModel agrees with the later calculation to within a few degrees, while the 1985 captions differ, in one case with the opposite sign.

The angles above come from the oscillator solver; the finite-difference solver agrees within 0.1°. Because the off-diagonal element is a small difference of two larger terms, the angle is a sensitive test of both implementations, and the heavy test suite compares them. The spectrum masses agree with the 1985 figures throughout. [Mixing angles: a paper erratum](@ref) explains how this was traced to the 1985 caption values rather than to the implementation.

## What to take away

- Singlet–triplet mixing requires unequal masses and a state of each spin with the same $J$.
- Its size comes from a small off-diagonal element divided by a small diagonal gap, so it is much more sensitive than the masses.
- The sign of the angle depends on the quark order and the basis phase; always state the convention when quoting an angle.
