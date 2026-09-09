# Public API documentation graph

Explore 24 public entries along the input → spectrum → wavefunction workflow.
Arrows lead from a help entry to a related entry linked in its docstring.
Start with the suggested reading paths, or expand the complete graph.


## Suggested reading paths

A selected subset of actual links makes the main workflow easier to read.

```mermaid
flowchart LR
    QuarkMassTable --> load_quark_masses
    load_quark_masses --> load_parameters_and_quark_masses
    load_parameters --> load_parameters_and_quark_masses
    load_parameters_and_quark_masses --> Meson
    ConstituentMasses --> Meson
    Meson --> reduced_mass
    Meson --> flavor_label
    Meson --> is_equal_flavor
    Meson --> spectrum_levels
    spectrum_levels --> BasisState
    Meson --> compute_spectrum
    compute_spectrum --> MixedSpectrum
    MixedSpectrum --> MixedState
    MixedSpectrum --> spectrum_state
    spectrum_state --> physical_components
    MixedState --> StateMixing
    physical_components --> RadialWave
    physical_components --> physical_state_amplitude
    spectrum_state --> radial_wave
    radial_wave --> sample_wave
    radial_wave --> wave_norm
    BasisState --> compute_isoscalar_spectrum
```

## Complete graph

<details>
<summary>Show all core documentation links</summary>

```mermaid
flowchart LR
    load_parameters["load_parameters"]
    load_quark_masses["load_quark_masses"]
    load_parameters_and_quark_masses["load_parameters_and_quark_masses"]
    GIParameters["GIParameters"]
    QuarkMassTable["QuarkMassTable"]
    ConstituentMasses["ConstituentMasses"]
    Meson["Meson"]
    reduced_mass["reduced_mass"]
    flavor_label["flavor_label"]
    is_equal_flavor["is_equal_flavor"]
    BasisState["BasisState"]
    spectrum_levels["spectrum_levels"]
    compute_spectrum["compute_spectrum"]
    MixedSpectrum["MixedSpectrum"]
    MixedState["MixedState"]
    StateMixing["StateMixing"]
    spectrum_state["spectrum_state"]
    radial_wave["radial_wave"]
    physical_components["physical_components"]
    RadialWave["RadialWave"]
    sample_wave["sample_wave"]
    wave_norm["wave_norm"]
    physical_state_amplitude["physical_state_amplitude"]
    compute_isoscalar_spectrum["compute_isoscalar_spectrum"]
    BasisState --> compute_isoscalar_spectrum
    BasisState --> compute_spectrum
    BasisState --> physical_components
    BasisState --> spectrum_levels
    BasisState --> spectrum_state
    ConstituentMasses --> GIParameters
    ConstituentMasses --> Meson
    ConstituentMasses --> QuarkMassTable
    ConstituentMasses --> reduced_mass
    GIParameters --> Meson
    GIParameters --> MixedSpectrum
    GIParameters --> compute_spectrum
    GIParameters --> load_parameters
    GIParameters --> load_parameters_and_quark_masses
    GIParameters --> spectrum_levels
    GIParameters --> spectrum_state
    Meson --> ConstituentMasses
    Meson --> QuarkMassTable
    Meson --> compute_spectrum
    Meson --> flavor_label
    Meson --> is_equal_flavor
    Meson --> load_parameters
    Meson --> load_quark_masses
    Meson --> physical_components
    Meson --> reduced_mass
    Meson --> spectrum_levels
    Meson --> spectrum_state
    MixedSpectrum --> Meson
    MixedSpectrum --> MixedState
    MixedSpectrum --> StateMixing
    MixedSpectrum --> compute_isoscalar_spectrum
    MixedSpectrum --> compute_spectrum
    MixedSpectrum --> physical_components
    MixedSpectrum --> radial_wave
    MixedSpectrum --> spectrum_state
    MixedState --> MixedSpectrum
    MixedState --> StateMixing
    MixedState --> compute_spectrum
    MixedState --> physical_components
    MixedState --> radial_wave
    MixedState --> spectrum_state
    QuarkMassTable --> ConstituentMasses
    QuarkMassTable --> Meson
    QuarkMassTable --> load_parameters_and_quark_masses
    QuarkMassTable --> load_quark_masses
    RadialWave --> physical_components
    RadialWave --> radial_wave
    RadialWave --> sample_wave
    RadialWave --> wave_norm
    StateMixing --> MixedState
    compute_isoscalar_spectrum --> BasisState
    compute_isoscalar_spectrum --> Meson
    compute_isoscalar_spectrum --> MixedSpectrum
    compute_isoscalar_spectrum --> MixedState
    compute_isoscalar_spectrum --> compute_spectrum
    compute_isoscalar_spectrum --> physical_components
    compute_isoscalar_spectrum --> spectrum_state
    compute_spectrum --> Meson
    compute_spectrum --> MixedSpectrum
    compute_spectrum --> MixedState
    compute_spectrum --> load_parameters_and_quark_masses
    compute_spectrum --> physical_components
    compute_spectrum --> spectrum_levels
    compute_spectrum --> spectrum_state
    flavor_label --> Meson
    flavor_label --> compute_spectrum
    flavor_label --> is_equal_flavor
    is_equal_flavor --> Meson
    is_equal_flavor --> compute_spectrum
    is_equal_flavor --> flavor_label
    load_parameters --> GIParameters
    load_parameters --> Meson
    load_parameters --> compute_spectrum
    load_parameters --> load_parameters_and_quark_masses
    load_parameters_and_quark_masses --> GIParameters
    load_parameters_and_quark_masses --> Meson
    load_parameters_and_quark_masses --> QuarkMassTable
    load_parameters_and_quark_masses --> compute_spectrum
    load_parameters_and_quark_masses --> load_parameters
    load_parameters_and_quark_masses --> load_quark_masses
    load_quark_masses --> Meson
    load_quark_masses --> QuarkMassTable
    load_quark_masses --> load_parameters_and_quark_masses
    physical_components --> MixedState
    physical_components --> RadialWave
    physical_components --> StateMixing
    physical_components --> compute_spectrum
    physical_components --> physical_state_amplitude
    physical_components --> radial_wave
    physical_components --> sample_wave
    physical_components --> spectrum_state
    physical_state_amplitude --> physical_components
    radial_wave --> MixedState
    radial_wave --> RadialWave
    radial_wave --> compute_spectrum
    radial_wave --> physical_components
    radial_wave --> sample_wave
    radial_wave --> spectrum_state
    radial_wave --> wave_norm
    reduced_mass --> Meson
    reduced_mass --> compute_spectrum
    reduced_mass --> spectrum_state
    sample_wave --> RadialWave
    sample_wave --> physical_components
    sample_wave --> radial_wave
    sample_wave --> wave_norm
    spectrum_levels --> BasisState
    spectrum_levels --> Meson
    spectrum_levels --> compute_spectrum
    spectrum_state --> BasisState
    spectrum_state --> MixedSpectrum
    spectrum_state --> MixedState
    spectrum_state --> compute_spectrum
    spectrum_state --> physical_components
    spectrum_state --> radial_wave
    wave_norm --> RadialWave
    wave_norm --> physical_components
    wave_norm --> radial_wave
```

</details>

## Link coverage

The selected entries contain **118** links. From `Meson`,
**24 of 24** entries are reachable (including itself).

| Entry | Incoming links | Outgoing links |
|---|---:|---:|
| `load_parameters` | 3 | 4 |
| `load_quark_masses` | 3 | 3 |
| `load_parameters_and_quark_masses` | 5 | 6 |
| `GIParameters` | 3 | 7 |
| `QuarkMassTable` | 4 | 4 |
| `ConstituentMasses` | 2 | 4 |
| `Meson` | 13 | 11 |
| `reduced_mass` | 2 | 3 |
| `flavor_label` | 2 | 3 |
| `is_equal_flavor` | 2 | 3 |
| `BasisState` | 3 | 5 |
| `spectrum_levels` | 4 | 3 |
| `compute_spectrum` | 15 | 7 |
| `MixedSpectrum` | 5 | 8 |
| `MixedState` | 7 | 6 |
| `StateMixing` | 3 | 1 |
| `spectrum_state` | 10 | 6 |
| `radial_wave` | 7 | 7 |
| `physical_components` | 12 | 8 |
| `RadialWave` | 4 | 4 |
| `sample_wave` | 3 | 4 |
| `wave_norm` | 3 | 3 |
| `physical_state_amplitude` | 1 | 1 |
| `compute_isoscalar_spectrum` | 2 | 7 |

Counts combine constructor and method docstrings and exclude self-links.
This graph covers explicit docstring links among the selected entries;
it does not represent function calls or every exported method.

See [Discoverability](discoverability.md) for help conventions and regeneration.
