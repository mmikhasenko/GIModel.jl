```@raw html
---
layout: home

hero:
  name: "GIModel.jl"
  text: "Mesons from the Godfrey–Isgur quark model"
  tagline: Masses, wavefunctions and decays of quark–antiquark mesons, from the relativized Hamiltonian of Phys. Rev. D 32, 189 (1985), in Julia.
  image:
    src: /home/wavefunction_cloud.png
    alt: Charmonium 1¹D₂ wavefunction drawn as a cloud of points
  actions:
    - theme: brand
      text: Get started
      link: /getting_started
    - theme: alt
      text: Tutorials
      link: /tutorials/charmonium
    - theme: alt
      text: View on GitHub
      link: https://github.com/mmikhasenko/GIModel.jl

features:
  - icon: M
    title: Masses
    details: Every n ²ˢ⁺¹L_J level of any u, d, s, c, b pair, with its mass split into central, contact, spin-orbit and tensor parts, and the mixing between levels.
    link: /manual/spectra
  - icon: ψ
    title: Wavefunctions
    details: Position- and momentum-space wavefunctions for every state, including the signed flavor and spin components of mixed states.
    link: /manual/wavefunctions
  - icon: γ
    title: Transitions
    details: Radiative E1, M1 and M2, strong pseudoscalar emission, γγ, gg and ggg annihilation, and leptonic widths, all from the solved wavefunctions.
    link: /manual/transitions
  - icon: "85"
    title: Checked against the paper
    details: Spectra of all seven flavor sectors within a few MeV, and every table of the 1985 paper reproduced, with each difference documented.
    link: /paper/results
---
```

````@raw html
<div class="gi-showcase">

<h2>Masses</h2>
<div class="gi-row">
<div>

One call solves the relativized Hamiltonian for a quark–antiquark pair and
returns every requested level, with the paper's harmonic-oscillator method and
an automatic convergence check. Here, the spectrum of the B<sub>c</sub> mesons:

```julia
using GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
Bc = compute_spectrum(params, Meson(mq, :b, :c);
    levels = spectrum_levels(2), solver = OscillatorSolver())
spectrum_state(Bc, "1^3S_1").mass_GeV     # 6.333
```

The two J<sup>P</sup> = 1<sup>+</sup> states are superpositions of the spin
singlet and triplet, because the quark masses differ.
[Computing a spectrum →](/manual/spectra)

</div>
<img src="/home/bc_spectrum.png" alt="The B_c meson spectrum computed with GIModel">
</div>

<h2>Wavefunctions</h2>
<div class="gi-row">
<div>

The cloud at the top of this page is the charmonium 1¹D₂ state: points drawn
from |ψ(<b>r</b>)|², coloured by the sign of ψ. Every state carries its radial
wavefunction in position and momentum space, and mixed states expose each
signed component.

[Wavefunctions →](/manual/wavefunctions)

</div>
<div>

```julia
spec = compute_spectrum(params, Meson(mq, :c, :c);
    levels = spectrum_levels(2), solver = OscillatorSolver())
wave = radial_wave(spec, "1^1D_2")
sample_wave(wave, range(0, 12; length = 241))   # u(r) on a grid
physical_components(spec, "1^3S_1")             # J/ψ: 1S with a small D admixture
```

</div>
</div>

<h2>Transitions</h2>
<div class="gi-row">
<div>

```julia
using GIModel.QuarkModelTransitions
psi, eta_c = physical_state(spec, "1^3S_1"), physical_state(spec, "1^1S_0")
decay_width(matrix_element(eta_c, PhotonEmission(mq), psi))       # M1, MeV
decay_width(MasslessLeptonPair(), LeptonicCurrent(:electromagnetic, mq), psi)
```

</div>
<div>

Decay amplitudes use the same solved states, so mixing and radial structure
enter coherently: photon emission (E1, M1, M2), pseudoscalar emission with all
partial waves, annihilation into γγ, gg and ggg, and leptonic decays.
[Transitions and decays →](/manual/transitions)

</div>
</div>

<h2>Get it</h2>

```julia
using Pkg
Pkg.add(url = "https://github.com/mmikhasenko/GIModel.jl")
```

Then follow [Getting started](/getting_started). If you use GIModel in
published work, please cite S. Godfrey and N. Isgur, *Phys. Rev. D* **32**, 189
(1985), and this repository.

</div>
````
