# Transitions and decays

`GIModel.QuarkModelTransitions` computes decay amplitudes and widths from the
states and wavefunctions of a GIModel spectrum. It is a submodule of GIModel,
so there is nothing extra to install:

```@example tr
using GIModel
using GIModel.QuarkModelTransitions
params, mq = load_parameters_and_quark_masses(default_parameters_path())
spec = compute_spectrum(params, Meson(mq, :c, :c);
    levels = spectrum_levels(2; L_labels = ("S", "P")), solver = OscillatorSolver())
nothing # hide
```

## The pattern

Every calculation follows the same four steps:

1. **Resolve** the states: [`physical_state`](@ref)`(spectrum, label)`.
2. **Choose** an operator: photon emission, pseudoscalar emission, a leptonic
   current, or an annihilation.
3. **Evaluate** the amplitude: [`matrix_element`](@ref)`(final, operator, initial)`.
4. **Convert** to a width in MeV: [`decay_width`](@ref).

```@example tr
psi = physical_state(spec, "1^3S_1")
eta_c = physical_state(spec, "1^1S_0")
amplitude = matrix_element(eta_c, PhotonEmission(mq), psi)
```

```@example tr
decay_width(amplitude)       # MeV
```

The argument order is always `(final, operator, initial)`, like
``\langle f|\,O\,|i\rangle``.

### Physical states

A [`PhysicalState`](@ref) is a self-contained copy of one meson: its label, mass,
``J^P``, and the signed components with their radial waves.

```@example tr
psi
```

Because it holds all components, a mixed parent or daughter contributes
coherently, with interference, and nothing needs to be re-assembled by hand.
Masses come from the spectrum, so decay momenta follow from the model
masses. [Comparing at other masses](@ref) shows how to use different masses.

## What can be computed

| operator | process | final state for `matrix_element` | final state for `decay_width` |
|---|---|---|---|
| [`PhotonEmission`](@ref) | ``A \to B\gamma`` (E1, M1, M2) | daughter `PhysicalState` | the same |
| [`PseudoscalarEmission`](@ref) | ``A \to B\,P`` | [`TwoMesonChannel`](@ref)`(B, P)` | the same |
| [`LeptonicCurrent`](@ref) | ``V \to e^+e^-``, ``P \to \ell\nu`` | [`Vacuum`](@ref)`()` | [`MasslessLeptonPair`](@ref)`()`, [`LeptonNeutrinoChannel`](@ref) |
| [`TwoPhotonAnnihilation`](@ref) | ``{}^1S_0, {}^3P_2 \to \gamma\gamma`` | [`TwoPhotonChannel`](@ref)`()` | the same |
| [`GluonicAnnihilation`](@ref) | ``\to gg``, ``\to ggg`` | not available | [`TwoGluonChannel`](@ref)`()`, [`ThreeGluonChannel`](@ref)`()` |

Widths are always in **MeV**. The amplitude units depend on the process:

| result | `.value` is |
|---|---|
| M1 photon amplitude | magnetic moment in nuclear magnetons |
| E1 and M2 photon amplitude | width amplitude in ``\mathrm{MeV}^{1/2}`` |
| pseudoscalar emission | dimensionless partial-wave amplitudes (see below) |
| leptonic current | dimensionless reduced current |
| two-photon | width amplitude in ``\mathrm{GeV}^{1/2}`` |

## Photon emission

[`PhotonEmission`](@ref) follows Appendix D of the paper. You do not name the
multipole: it is inferred from the two states. With `verbose = true` the
selection is logged:

```@example tr
chi_c1 = physical_state(spec, "1^3P_1")
e1 = matrix_element(psi, PhotonEmission(mq), chi_c1; verbose = true)
decay_width(e1)
```

The classes are `DirectM1` (same radial level), `HinderedM1` (different radial
levels), `AllowedE1` (spin-conserving S–P), `SpinFlipE1` and `SpinFlipM2`.
`PhotonEmission(mq; recoil_order = 2)` adds the ``(qr)^2`` correction where it
is implemented (currently M1).

The photon couples to the quark charges, so flavor must be explicit. States
built from `:u` and `:d` work; the averaged `:q` is rejected for photon
emission because it cannot distinguish ``e_u + e_d`` from ``e_u - e_d``.

## Pseudoscalar emission

[`PseudoscalarEmission`](@ref)`(g, h, masses)` is the elementary-emission
operator of Eq. (19),

```math
H_\text{int} \propto g\,\boldsymbol\sigma\cdot\boldsymbol q + h\,\boldsymbol\sigma\cdot\boldsymbol p',
```

where ``\boldsymbol q`` is the momentum of the emitted pseudoscalar and
``\boldsymbol p'`` the momentum of the quark that emits it. The couplings `g`
and `h` (in ``\mathrm{GeV}^{-1}``) are **not predicted** by the spectrum; you
supply them, typically fitted to one or two known widths.

The final channel is ordered: `TwoMesonChannel(surviving, emitted)`, and the
emitted meson must be a ``0^-`` state. The result holds one amplitude per
allowed partial wave ``(L, S)``:

```@example tr
# K*⁺ (u s̄) → K⁰ (d s̄) + π⁺ (u d̄), with unit direct coupling g and h = 0.
S1, S0 = BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
solver = OscillatorSolver()
us = compute_spectrum(params, Meson(mq, :u, :s); levels = [S1], solver)
ds = compute_spectrum(params, Meson(mq, :d, :s); levels = [S0], solver)
ud = compute_spectrum(params, Meson(mq, :u, :d); levels = [S0], solver)
K_star, K_zero, pion = physical_state(us, S1), physical_state(ds, S0), physical_state(ud, S0)
strong = matrix_element(TwoMesonChannel(K_zero, pion), PseudoscalarEmission(1.0, 0.0, mq), K_star)
```

```@example tr
[w => strong[w] for w in partial_waves(strong)]
```

A vector decaying to two pseudoscalars has one P-wave amplitude. The flavor
algebra is exact: a channel that violates charge or flavor, such as
``K^{*+} \to K^+\pi^+``, gives exactly zero. See
[Strong decays beyond the paper](@ref) for neutral isospin states, axial
mesons with two partial waves, and a calibration of `g` and `h`.

This operator uses the full solved wavefunctions of parent and daughter, with
all radial nodes and each state's own size. The paper's Table V instead used a
single oscillator scale ``\beta = 0.40`` GeV for every meson.
[The 1985 paper](@ref "Reproducing Godfrey–Isgur") section explains how that
table is reproduced separately.

## Leptonic decays

A vector meson annihilates into a lepton pair through the electromagnetic
current. `LeptonicCurrent(:electromagnetic, masses)` takes the charges from the
flavor components:

```@example tr
current = LeptonicCurrent(:electromagnetic, mq)
decay_width(MasslessLeptonPair(), current, psi)     # Γ(J/ψ → e⁺e⁻) in MeV
```

`matrix_element(Vacuum(), current, psi)` returns the dimensionless reduced
current instead of a width. Weak decays of pseudoscalars use a
[`LeptonNeutrinoChannel`](@ref) with an explicit CKM element.

## Annihilation into photons and gluons

[`TwoPhotonAnnihilation`](@ref) takes explicit charge factors as
[`AnnihilationTerm`](@ref)s: for ``c\bar c`` the effective squared charge is
``e_c^2 = 4/9``.

```@example tr
two_photon = TwoPhotonAnnihilation(mq, AnnihilationTerm((:c, :c), 4 / 9))
decay_width(TwoPhotonChannel(), two_photon, eta_c)  # Γ(η_c → γγ) in MeV
```

[`GluonicAnnihilation`](@ref) takes the strong coupling, either as a number or
as a function of the meson mass:

```@example tr
gluonic = GluonicAnnihilation(mq, M -> alpha_s_q(M))
(eta_c = decay_width(TwoGluonChannel(), gluonic, eta_c),
 psi = decay_width(ThreeGluonChannel(), gluonic, psi))
```

These annihilation widths depend on the smeared wavefunction at the origin
(Eq. (17)), which is sensitive to the short-distance part of the wave. Check
their convergence when you quote them.

## Comparing at other masses

`matrix_element` and `decay_width` use the model masses stored in the states.
To see how an amplitude would change at an experimental mass or momentum,
keeping the wavefunctions fixed, use [`mass_correction_factor`](@ref):

```@example tr
factor = mass_correction_factor(TwoPhotonChannel(), two_photon, eta_c; target_mass = 2.984)
```

The factor multiplies the *amplitude* (`.value`), not the width, except for
gluonic annihilation, where no amplitude exists. The required keyword depends
on the operator:

| operator | keyword |
|---|---|
| `LeptonicCurrent`, `TwoPhotonAnnihilation`, `GluonicAnnihilation` | `target_mass` (GeV) |
| `PhotonEmission` | `target_momentum` (GeV) |
| `PseudoscalarEmission` | `target_momentum` (GeV) and `partial_wave` |

## States you build yourself

A `PhysicalState` can also be constructed by hand from any radial waves, for
example a simple oscillator function for an analytic check, or a linear
combination of solved states for an isospin eigenstate:

```@example tr
wave = OscillatorWave(0, 0.5, [1.0])     # one oscillator function, β = 0.5 GeV
toy = PhysicalState("toy ψ", 3.10, [(
    basis = BasisState(1, "S", 3, 1; flavors = (:c, :c)), coefficient = 1.0, wave = wave)])
toy.mass_GeV, toy.J
```

The package does not check such states: their masses and waves are your
assumptions.

## Charge radii

[`charge_radius_squared`](@ref) computes the mean-square charge radius of a
``{}^1S_0`` meson from its radial wave (Table VII(d) of the paper), in
``\mathrm{GeV}^{-2}``:

```@example tr
pion_wave = radial_wave(ud, "1^1S_0")
r2 = charge_radius_squared(pion_wave, mq["u"], 2 / 3, mq["d"], 1 / 3)
r2 * 0.1973^2            # fm²
```

## What is not covered

The operators above are the Godfrey–Isgur ones. The module does not implement
quark-pair-creation (``{}^3P_0``) decay models, emission of mesons other than
pseudoscalars, weak hadronic or semileptonic transitions, or coupled-channel
effects. See [Scope and limitations](@ref).
