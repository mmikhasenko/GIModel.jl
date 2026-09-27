```@meta
EditURL = "../../quarto/manual/isoscalar.qmd"
```



# Isoscalar flavor mixing

Self-conjugate isoscalar mesons, such as $\eta$, $\eta'$, $\omega$, $\phi$, $f_2$ and $f_2'$, are not pure flavor states. A $q\bar q$ pair can annihilate into gluons and reappear as a different flavor, which mixes $n\bar n = (u\bar u + d\bar d)/\sqrt2$ with $s\bar s$ (and, weakly, with $c\bar c$ and $b\bar b$). This is the third stage of the Godfrey–Isgur calculation (Sec. II D and Table III of the paper).

## Why it needs its own step

[`compute_spectrum`](@ref) solves one flavor channel. Annihilation connects two or more channels, so it needs their solved states together. The package therefore:

1.  solves each flavor channel through stages 1 and 2;
2.  builds the annihilation matrix between states with the same quantum numbers in different channels;
3.  diagonalizes it and returns one combined spectrum.

[`compute_isoscalar_spectrum`](@ref) does all three for the $n\bar n$–$s\bar s$ pair. [`add_isoscalar_annihilation`](@ref) does steps 2 and 3 for spectra you have already computed, including combinations with heavy flavors.

## What must be specified

Annihilation is *not* applied by default: the paper uses different prescriptions for different $J^{PC}$, so each one is requested explicitly.

- `pseudoscalar`: the model for ${}^1S_0$ states. Choose [`PaperP1Annihilation`](@ref) (Eq. (18a)) or [`PaperP2Annihilation`](@ref) (Eq. (18b)). The amplitude involves the smeared wavefunction at the origin (Eq. (17)), computed from the solved waves.
- `amplitudes`: a dictionary from `(L_label, multiplicity, J)` to the amplitude $A$ of the general Eq. (16), applied at every radial level present in both channels. Table III gives $A({}^3S_1) = 2.5$ and $A({}^3P_2) = -0.8$; both are stored in `params.annihilation`.

Channels that are not listed stay ideally mixed (pure $n\bar n$ or $s\bar s$), as the Table III caption prescribes.

## The η and η′

```julia
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
nn, ss = Meson(mq, :q, :q), Meson(mq, :s, :s)
A = params.annihilation
iso = compute_isoscalar_spectrum(params, nn, ss;
    levels = spectrum_levels(2; L_labels = ("S",)),
    pseudoscalar = PaperP1Annihilation(),
    amplitudes = Dict(("S", 3, 1) => A.s1_A))
```

```
MixedSpectrum: qq + ss, 8 levels — all values in GeV
  level        central    contact   fine str     mixing       mass
  1^1S_0 qq     0.7637    -0.6151     0.0000     0.3828     0.5313
  1^3S_1 qq     0.6671     0.1041     0.0000     0.0099     0.7810
  2^1S_0 qq     1.3966    -0.1050     0.0000     0.2040     1.4957
  2^3S_1 qq     1.4113     0.0443     0.0000     0.0013     1.4569
  1^1S_0 ss     0.9950    -0.3402     0.0000     0.3456     1.0004
  1^3S_1 ss     0.9491     0.0679     0.0000     0.0033     1.0203
  2^1S_0 ss     1.6531    -0.0761     0.0000     0.2229     1.7999
  2^3S_1 ss     1.6596     0.0293     0.0000     0.0006     1.6895
  (8 levels carry mixing; see `spec.states[i].mixings`)
```

The combined spectrum holds both channels. Every row is labeled with the channel it came from: each mixed eigenstate keeps the label of the unmixed state it is assigned to, by ascending mass. The `mixing` column shows how far annihilation moved each level. Pseudoscalars move by hundreds of MeV; vectors move by a few MeV.

Because labels repeat across channels, select states with a flavor-qualified [`BasisState`](@ref):

```julia
eta = spectrum_state(iso, BasisState(1, "S", 1, 0; flavors = (:q, :q)))
eta_prime = spectrum_state(iso, BasisState(1, "S", 1, 0; flavors = (:s, :s)))
(eta.mass_GeV, eta_prime.mass_GeV)
```

```
(0.5313070518403085, 1.0004316903836472)
```

### Flavor and radial composition

The annihilation block contains every ${}^1S_0$ level that was requested, so the physical $\eta$ mixes ground and radially excited states of both flavors:

```julia
[(join(c.basis.flavors), c.basis.label, round(c.coefficient; digits = 3))
 for c in physical_components(iso, eta)]
```

```
4-element Vector{Tuple{String, String, Float64}}:
 ("qq", "1^1S_0", 0.624)
 ("ss", "1^1S_0", -0.772)
 ("qq", "2^1S_0", 0.11)
 ("ss", "2^1S_0", 0.041)
```

The $\eta$ is octet-like, with $n\bar n$ and $s\bar s$ of opposite sign. The $\eta'$ has same-sign components:

```julia
[(join(c.basis.flavors), c.basis.label, round(c.coefficient; digits = 3))
 for c in physical_components(iso, eta_prime)]
```

```
4-element Vector{Tuple{String, String, Float64}}:
 ("qq", "1^1S_0", 0.568)
 ("ss", "1^1S_0", 0.552)
 ("qq", "2^1S_0", 0.592)
 ("ss", "2^1S_0", 0.153)
```

The vectors stay close to ideal mixing: the $\omega$ is nearly pure $n\bar n$ and the $\phi$ nearly pure $s\bar s$:

```julia
omega = spectrum_state(iso, BasisState(1, "S", 3, 1; flavors = (:q, :q)))
[(join(c.basis.flavors), c.basis.label, round(c.coefficient; digits = 3))
 for c in physical_components(iso, omega)]
```

```
2-element Vector{Tuple{String, String, Float64}}:
 ("qq", "1^3S_1", 1.0)
 ("ss", "1^3S_1", -0.024)
```

## Choosing the radial content

As with spectroscopic mixing, a block contains only the states it is given, and the result can depend strongly on that choice.

- For vectors, tensors and any other `amplitudes` entry, a block is built at every radial level present in both channels.
- For pseudoscalars, `compute_isoscalar_spectrum` always uses the paper’s four-state basis $[1n\bar n, 1s\bar s, 2n\bar n, 2s\bar s]$, so the requested levels must include $1\,{}^1S_0$ and $2\,{}^1S_0$.

To use another pseudoscalar basis, solve the channels yourself and pass the basis explicitly to [`add_isoscalar_annihilation`](@ref), as in the next section. The [η and η′](@ref) tutorial compares the two-state and four-state results: the $\eta'$ moves by about 300 MeV.

## Heavy flavors and combined spectra

For more than two flavors, solve each channel with the same parameters and solver, then combine them. The pseudoscalar block needs an explicit list of flavor-qualified basis states:

```julia
solver = OscillatorSolver()
ground = [BasisState(1, "S", 1, 0)]
spectra = [compute_spectrum(params, Meson(mq, f, f); levels = ground, solver)
           for f in (:q, :s, :c)]
basis = [BasisState(1, "S", 1, 0; flavors = (f, f)) for f in (:q, :s, :c)]
four = add_isoscalar_annihilation(params, spectra;
    pseudoscalar = PaperP1Annihilation(), pseudoscalar_basis = basis)
```

```
MixedSpectrum: qq + ss + cc, 3 levels — all values in GeV
  level        central    contact   fine str     mixing       mass
  1^1S_0 qq     0.7635    -0.6145     0.0000     0.3932     0.5423
  1^1S_0 ss     0.9950    -0.3398     0.0000     0.6226     1.2777
  1^1S_0 cc     3.0785    -0.1114     0.0000    -0.0023     2.9648
  (3 levels carry mixing; see `spec.states[i].mixings`)
```

The $\eta_c$ shifts by a few MeV and gains small light-quark components:

```julia
eta_c = spectrum_state(four, BasisState(1, "S", 1, 0; flavors = (:c, :c)))
[(join(c.basis.flavors), round(c.coefficient; digits = 4)) for c in physical_components(four, eta_c)]
```

```
3-element Vector{Tuple{String, Float64}}:
 ("qq", -0.0069)
 ("ss", -0.004)
 ("cc", 1.0)
```

## Using isoscalar states in transitions

The isoscalar spectrum records that its nonstrange channel means $(u\bar u+d\bar d)/\sqrt{2}$. Photon emission and electromagnetic leptonic currents resolve this channel into explicit `:u` and `:d` contributions, including their relative signs and mixing coefficients:

```julia
using GIModel.QuarkModelTransitions
phi = physical_state(iso, BasisState(1, "S", 3, 1; flavors = (:s, :s)))
decay_width(MasslessLeptonPair(), LeptonicCurrent(:electromagnetic, mq), phi)
```

```
0.002844628833069689
```

A standalone `compute_spectrum(..., Meson(mq, :q, :q))` remains ambiguous and is rejected by electromagnetic currents. Explicit `AnnihilationTerm` coefficients, including those for two-photon widths, still act on the original coarse flavor components, so effective charges are not counted twice.

## A calibrated alternative

[`CalibratedP1Annihilation`](@ref) fixes the pseudoscalar block so that its eigenvalues equal target masses you supply through `pseudoscalar_targets`. It is used to compare with the paper’s quoted masses when studying observables that depend on the mixing, and it is not a prediction.
