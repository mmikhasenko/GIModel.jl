# Computing a spectrum

[`compute_spectrum`](@ref) is the main entry point of the package. This page
explains what you pass to it, what it returns, and how to take the result
apart.

```@example spectra
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
nothing # hide
```

## Choosing levels

The `levels` keyword says which states to compute. It is a vector of
[`BasisState`](@ref)s, each naming one ``n\,{}^{2S+1}L_J`` multiplet by its
radial number `n`, orbital letter `L_label`, spin multiplicity ``2S+1`` and
total ``J``:

```@example spectra
BasisState(1, "P", 3, 2)        # 1³P₂, e.g. χ_c2
```

[`spectrum_levels`](@ref) builds the complete list for given radial numbers
and orbital letters. By default it includes S, P and D waves:

```@example spectra
[l.label for l in spectrum_levels(1)]
```

```@example spectra
[l.label for l in spectrum_levels(2; L_labels = ("S",))]
```

A hand-written vector works the same way. Computing only what you need is
faster, but it also changes the mixing blocks, as explained
[below](@ref "Mixing depends on the requested levels").

## The result: a `MixedSpectrum`

```@example spectra
ccbar = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(1))
```

The result is a [`MixedSpectrum`](@ref). Its `states` are
[`MixedState`](@ref)s, stored in the order you requested them. Retrieve one by
label, by quantum numbers, or by `BasisState`:

```@example spectra
spectrum_state(ccbar, "1^3P_2").mass_GeV
```

```@example spectra
spectrum_state(ccbar, 1, "P", 3, 2) === spectrum_state(ccbar, BasisState(1, "P", 3, 2))
```

A missing level throws an `ArgumentError` instead of returning `nothing`.

### What a state contains

```@example spectra
chi_c2 = spectrum_state(ccbar, "1^3P_2")
propertynames(chi_c2)
```

| property | meaning |
|---|---|
| `mass_GeV` | physical mass after all mixing |
| `n`, `L`, `multiplicity`, `J`, `label` | quantum numbers of the assigned level |
| `basis` | the same as a [`BasisState`](@ref), including `flavors` |
| `central_GeV` | expectation value of the spin-independent Hamiltonian |
| `contact_shift_GeV` | expectation value of the contact hyperfine term |
| `spin_orbit_vector_shift_GeV`, `spin_orbit_thomas_shift_GeV` | the two spin-orbit pieces |
| `spin_orbit_shift_GeV` | their sum |
| `tensor_shift_GeV` | diagonal tensor term |
| `fine_structure_shift_GeV` | spin-orbit plus tensor |
| `corrected` | the state before mixing (a [`CorrectedState`](@ref)) |
| `mixings` | records of every mixing step (empty if none) |

The spin-orbit and tensor pieces for ``\chi_{c2}``:

```@example spectra
(vector = chi_c2.spin_orbit_vector_shift_GeV,
 thomas = chi_c2.spin_orbit_thomas_shift_GeV,
 tensor = chi_c2.tensor_shift_GeV)
```

The vector (one-gluon) and Thomas (confinement) spin-orbit terms have opposite
signs and partly cancel, as described in [The model](@ref).

## The stages

`compute_spectrum` is two calls in sequence:

```julia
compute_spectrum(params, meson; kw...) ==
    add_intra_meson_mixing(fixed_spectrum(params, meson; kw...))
```

Each stage returns a [`Spectrum`](@ref) with its own state type, so you can
stop after either one:

| stage | function | returns | state type |
|---|---|---|---|
| 1: solve each ``(L, S, J)`` | [`fixed_spectrum`](@ref) | [`CorrectedSpectrum`](@ref) | [`CorrectedState`](@ref) |
| 2: tensor and spin-orbit mixing | [`add_intra_meson_mixing`](@ref) | [`MixedSpectrum`](@ref) | [`MixedState`](@ref) |

```@example spectra
D = Meson(mq, :c, :u)
fixed = fixed_spectrum(params, D; levels = spectrum_levels(1; L_labels = ("P",)))
```

```@example spectra
add_intra_meson_mixing(fixed)
```

The two ``J = 1`` levels move apart after mixing: the lower goes down and the
upper goes up, as in any two-level system.

### Why the central column depends on spin

In stage 1, each ``(L, S, J)`` sector is diagonalized with all its spin terms
included, so the four 1P wavefunctions differ slightly. The `central` column
is the spin-independent energy evaluated in *that* wavefunction, which is why
it is not identical across a multiplet. [`central_spectrum`](@ref) is a
separate diagnostic that solves the spin-independent problem alone:

```@example spectra
central_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(1; L_labels = ("P",)))
```

Here every 1P level has the same energy. `compute_spectrum` never calls this
function; it is useful for teaching and for isolating the spin-independent
binding.

## Mixing records

When a state mixes, its `mixings` field holds one [`StateMixing`](@ref) per
mixing step:

```@example spectra
mixed = add_intra_meson_mixing(fixed)
d1 = spectrum_state(mixed, "1^1P_1")
m = only(d1.mixings)
(m.mechanism, m.block_label, m.partner_labels)
```

```@example spectra
m.components            # signed amplitudes on the partner basis states
```

```@example spectra
(m.offdiag_GeV, m.partner_masses_GeV)
```

`components` is this state's eigenvector in the basis `partner_labels`. The
lower ``D_1`` state here is about 81% ``{}^1P_1`` and 19% ``{}^3P_1``. The
shared block with all its eigenvectors is in `m.result`, a
[`MixingResult`](@ref).

Mixed eigenstates are assigned to labels by order: the lowest mixed mass takes
the label of the lowest unmixed mass, and so on. A mixed state therefore keeps
the label of the level it came from, even when the mixing is large.

For a complete, flattened decomposition, including chains of mixing steps,
use [`physical_components`](@ref) (see [Wavefunctions](@ref)):

```@example spectra
[(c.basis.label, round(c.coefficient; digits = 3)) for c in physical_components(mixed, d1)]
```

### Mixing depends on the requested levels

A mixing block contains only the levels you requested. With
`spectrum_levels(2)`, the ``J/\psi`` mixes with ``2\,{}^3S_1``, ``1\,{}^3D_1``
and ``2\,{}^3D_1``; with only ``1\,{}^3S_1`` and ``1\,{}^3D_1``, it mixes with
one partner:

```@example spectra
small = compute_spectrum(params, Meson(mq, :c, :c);
    levels = [BasisState(1, "S", 3, 1), BasisState(1, "D", 3, 1)])
large = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(2))
[spectrum_state(s, "1^3S_1").mixings[1].partner_labels for s in (small, large)]
```

For low-lying states the effect on masses is small, but compositions of
nearly degenerate levels can depend on it. Request every level that can mix
with the state you study.

## Switching spin terms off

[`SpinTerms`](@ref) turns individual spin-dependent terms on or off. Unlike the
solver settings, changing it is *meant* to change the answer, so it is how you
measure what each term does:

```@example spectra
no_contact = compute_spectrum(params, Meson(mq, :c, :c);
    levels = spectrum_levels(1; L_labels = ("S",)),
    terms = SpinTerms(contact_hyperfine = false))
```

Without the contact term, ``\eta_c`` and ``J/\psi`` are nearly degenerate;
the small difference that remains comes from tensor mixing with D waves. The
four switches are:

| switch | term | stage |
|---|---|---|
| `contact_hyperfine` | spin–spin contact | 1 |
| `fine_structure` | diagonal spin-orbit and tensor | 1 |
| `same_j_spin_orbit` | ``{}^1L_J \leftrightarrow {}^3L_J`` mixing (unequal masses only) | 2 |
| `tensor` | ``L = J\pm1`` mixing | 2 |

Heavy-light mixing can be switched off the same way, for comparison:

```@example spectra
unmixed = compute_spectrum(params, D;
    levels = spectrum_levels(1; L_labels = ("P",)),
    terms = SpinTerms(same_j_spin_orbit = false))
[s.label => round(s.mass_GeV; digits = 4) for s in unmixed.states]
```

## Choosing the solver

The `solver` keyword chooses the numerical method. The default is a
[`FiniteDifferenceSolver`](@ref); `solver = OscillatorSolver()` uses the
paper's harmonic-oscillator method. Both solve the same Hamiltonian;
[Solvers and convergence](@ref) explains when the choice matters.

## Parameters travel with the spectrum

A spectrum remembers the parameters and solver that produced it:

```@example spectra
parameters(ccbar) === params
```

```@example spectra
ccbar.computation.solver
```

Later calculations, such as annihilation mixing and transitions, read these
instead of asking for them again, so they cannot be mixed up by accident.
