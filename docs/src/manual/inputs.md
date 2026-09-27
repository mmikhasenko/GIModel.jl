```@meta
EditURL = "../../quarto/manual/inputs.qmd"
```



# Inputs and parameters

Every GIModel calculation takes three kinds of input. The package keeps them in separate objects because they answer different questions:

| kind | question it answers | object | changing it should… |
|----|----|----|----|
| **model** | which Hamiltonian? | [`GIParameters`](@ref), [`QuarkMassTable`](@ref) | change the physics |
| **numerics** | how accurately is it solved? | [`OscillatorSolver`](@ref), [`FiniteDifferenceSolver`](@ref) | *not* change the answer beyond the tolerance |
| **physics switches** | which terms are included? | [`SpinTerms`](@ref) | change the answer, on purpose |

This page covers the model inputs. [Solvers and convergence](@ref) covers the numerics and [Computing a spectrum](@ref) covers the switches.

## The shipped parameter set

The package ships the Godfrey–Isgur parameters (Table II of the paper) as a TOML file. [`default_parameters_path`](@ref) returns its location, and [`load_parameters_and_quark_masses`](@ref) reads it:

```julia
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
```

An ordinary meson spectrum depends on 18 numbers:

| inputs | count | values | where they live |
|----|---:|----|----|
| constituent masses $m_{u,d}$, $m_s$, $m_c$, $m_b$ | 4 | 0.220, 0.419, 1.628, 4.977 GeV | `[masses]` (stored in MeV) |
| confinement slope $b$ and offset $c$ | 2 | 0.18 GeV², −0.253 GeV | `[potential]` |
| smearing $\sigma_0$ and $s$ | 2 | 1.80 GeV, 1.55 | `[relativistic_smearing]` |
| relativistic exponents $\epsilon_c$, $\epsilon_t$, $\epsilon_{so(v)}$, $\epsilon_{so(s)}$ | 4 | −0.168, 0.025, −0.035, 0.055 | `[relativistic_factors]` |
| running-coupling weights $\alpha_k$ | 3 | 0.25, 0.15, 0.20 | fixed constants in the code |
| running-coupling scales $\gamma_k$ | 3 | 0.5, $\sqrt{10}/2$, $\sqrt{1000}/2$ GeV | fixed constants in the code |

The six coupling constants define the paper’s parametrization of $\alpha_s$ and are not configurable. The other twelve are read from the file.

The file has one more section, `[annihilation]`, with six constants for isoscalar annihilation mixing (Table III). They only enter [Isoscalar flavor mixing](@ref).

## `GIParameters`: the interaction

[`GIParameters`](@ref) groups the interaction by the paper’s structure. Each field is a small immutable object:

```julia
propertynames(params)
```

```
(:potential, :central, :smearing, :factors, :fine_structure, :annihilation)
```

```julia
params.potential
```

```
ConfinementPotential{Float64}(0.18, -0.253)
```

```julia
params.smearing
```

```
RelativisticSmearing{Float64}(1.8, 1.55)
```

```julia
params.factors
```

```
RelativisticFactors{Float64}(-0.168, 0.025, -0.035, 0.055, true, true, true)
```

| field | type | contents |
|----|----|----|
| `potential` | [`ConfinementPotential`](@ref) | $b$ in GeV², $c$ in GeV |
| `central` | [`CentralPotentialMethod`](@ref) | how the spin-independent potential is built |
| `smearing` | [`RelativisticSmearing`](@ref) | $\sigma_0$ in GeV, $s$ |
| `factors` | [`RelativisticFactors`](@ref) | the four $\epsilon_i$ and three switches for how they are applied |
| `fine_structure` | [`FineStructure`](@ref) | master switch for spin-orbit and tensor terms |
| `annihilation` | [`AnnihilationAmplitudes`](@ref) | Table III constants |

Quark masses are deliberately **not** in `GIParameters`: the interaction is the same for every flavor.

### The central potential method

`params.central` selects how the spin-independent potential is constructed. The shipped file uses [`AppendixAMomentumSandwich`](@ref), which is the paper’s prescription: closed-form smeared $\tilde G$ and $\tilde S$ with the Coulomb term enclosed by $A(p)$ factors, as in [The model](@ref). The other methods ([`PointwiseCentral`](@ref), [`AppendixAClosedForm`](@ref), and the diagnostic smearing variants) exist to study how much each ingredient matters.

```julia
params.central
```

```
AppendixAMomentumSandwich()
```

## Changing a parameter

All parameter objects are immutable. To change one number, build a copy with that field overridden. Each type has a constructor that takes an existing object and keyword overrides:

```julia
steeper = GIParameters(params;
    potential = ConfinementPotential(params.potential; b = 0.20))
steeper.potential
```

```
ConfinementPotential{Float64}(0.2, -0.253)
```

Everything else in `steeper` is the same as in `params`. This is the recommended way to vary the model, because the result stays complete and consistent.

**Building parameters from scratch.** The keyword constructor `GIParameters(; potential, smearing, ...)` fills unspecified fields with simple defaults: a pointwise central potential and no relativistic momentum factors. That is useful for studying a stripped-down model, but it is *not* the Godfrey–Isgur model. Start from the shipped file and override fields unless you want those defaults.

## Writing your own parameter file

To keep a modified parameter set, copy the shipped file and edit it:

```julia
cp(default_parameters_path(), "my_parameters.toml")
# edit my_parameters.toml
params, mq = load_parameters_and_quark_masses("my_parameters.toml")
```

The loader is strict. Every section and field of the ordinary model must be present, unknown fields are rejected, and every number must be finite, so a typo fails loudly. The `[metadata]` section is descriptive and does not enter the Hamiltonian. The `[annihilation]` section is optional; if it is omitted, the defaults of [`AnnihilationAmplitudes`](@ref) are used.

The loaders also exist separately: [`load_parameters`](@ref) reads only the interaction, and [`load_quark_masses`](@ref) reads only the `[masses]` section.

## Quark masses and mesons

A [`QuarkMassTable`](@ref) is a dictionary from flavor name to mass in GeV:

```julia
mq["s"], mq["c"]
```

```
(0.419, 1.628)
```

The light quarks are degenerate in the paper. The table provides `"u"`, `"d"` and `"q"` (the light average) with the same mass; use `:q` when isospin does not matter and `:u`/`:d` when the flavor must be explicit, for example in electromagnetic transitions.

A [`Meson`](@ref) combines two flavors with their masses. The first argument is the quark, the second the antiquark:

```julia
Ds = Meson(mq, :c, :s)
```

```
Meson: c sbar   (m1 = 1.628, m2 = 0.419 GeV, reduced = 0.3332 GeV)
```

```julia
(reduced_mass(Ds), flavor_label(Ds), is_equal_flavor(Ds))
```

```
(0.33323497801660973, "cs", false)
```

The order matters for signs: the antisymmetric spin-orbit operator changes sign under quark–antiquark exchange, and so does the mixing angle. [Conventions and units](@ref) states the convention.

### Masses that are not in the table

A meson can also be built from explicit masses, which is how you vary a quark mass continuously:

```julia
Meson(:c, :c, ConstituentMasses(1.5, 1.5))
```

```
Meson: c cbar   (m1 = 1.5, m2 = 1.5 GeV, reduced = 0.75 GeV)
  self-conjugate: same-J antisymmetric spin-orbit mixing vanishes
```

or from quark objects, which carry a mass and a flavor symbol:

```julia
Q = HeavyQuark{:up}(2.5, :Q)       # an up-type heavy quark with a free mass
Meson(Q, LightQuark(mq["q"]))
```

```
Meson: Q qbar   (m1 = 2.5, m2 = 0.22 GeV, reduced = 0.2022 GeV)
```

[`LightQuark`](@ref), [`StrangeQuark`](@ref) and [`HeavyQuark`](@ref) are subtypes of [`AbstractQuark`](@ref). The type parameter of `HeavyQuark` is its weak-isospin class (`:up` for charge $+2/3$, `:down` for $-1/3$); only electromagnetic transitions depend on it. [From charm to bottom](@ref) uses this to move a quark mass from the charm to the bottom value.

## Inputs outside the spectrum

Transition calculations need further inputs, such as the strong-decay couplings $g$ and $h$ and the photon-vertex smearing exponents. They belong to the transition operators, not to `GIParameters`; see [Transitions and decays](@ref).
