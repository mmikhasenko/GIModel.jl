# Public documentation graph

Generated from all 21 supported API names (all exported).
All 21 public docstrings have independently executable examples (21 blocks).
Every entry reaches every other entry through explicit help links. GIModel links are validated separately.

Regenerate and execute the examples from the repository root:

```sh
julia --project=. QuarkModelTransitions/scripts/audit_documentation.jl
```

The graph describes help navigation, not function calls or numerical dependencies.

| Entry | Access | Incoming | Outgoing |
|---|---|---:|---:|
| `AnnihilationTerm` | exported | 2 | 2 |
| `GluonicAnnihilation` | exported | 3 | 4 |
| `LeptonNeutrinoChannel` | exported | 2 | 2 |
| `LeptonicCurrent` | exported | 4 | 7 |
| `MasslessLeptonPair` | exported | 2 | 2 |
| `PartialWave` | exported | 3 | 1 |
| `PhotonEmission` | exported | 4 | 2 |
| `PhysicalState` | exported | 4 | 4 |
| `PseudoscalarEmission` | exported | 4 | 5 |
| `ThreeGluonChannel` | exported | 2 | 2 |
| `TwoGluonChannel` | exported | 2 | 2 |
| `TwoMesonChannel` | exported | 4 | 2 |
| `TwoPhotonAnnihilation` | exported | 4 | 5 |
| `TwoPhotonChannel` | exported | 3 | 2 |
| `Vacuum` | exported | 2 | 2 |
| `charge_radius_squared` | exported | 1 | 1 |
| `decay_width` | exported | 12 | 16 |
| `mass_correction_factor` | exported | 7 | 9 |
| `matrix_element` | exported | 18 | 15 |
| `partial_waves` | exported | 3 | 3 |
| `physical_state` | exported | 4 | 2 |

<details>
<summary>Complete graph (90 links)</summary>

```mermaid
flowchart LR
    AnnihilationTerm --> LeptonicCurrent
    AnnihilationTerm --> TwoPhotonAnnihilation
    GluonicAnnihilation --> ThreeGluonChannel
    GluonicAnnihilation --> TwoGluonChannel
    GluonicAnnihilation --> decay_width
    GluonicAnnihilation --> mass_correction_factor
    LeptonNeutrinoChannel --> decay_width
    LeptonNeutrinoChannel --> matrix_element
    LeptonicCurrent --> AnnihilationTerm
    LeptonicCurrent --> LeptonNeutrinoChannel
    LeptonicCurrent --> MasslessLeptonPair
    LeptonicCurrent --> Vacuum
    LeptonicCurrent --> decay_width
    LeptonicCurrent --> mass_correction_factor
    LeptonicCurrent --> matrix_element
    MasslessLeptonPair --> decay_width
    MasslessLeptonPair --> matrix_element
    PartialWave --> matrix_element
    PhotonEmission --> mass_correction_factor
    PhotonEmission --> matrix_element
    PhysicalState --> PhotonEmission
    PhysicalState --> TwoMesonChannel
    PhysicalState --> matrix_element
    PhysicalState --> physical_state
    PseudoscalarEmission --> TwoMesonChannel
    PseudoscalarEmission --> mass_correction_factor
    PseudoscalarEmission --> matrix_element
    PseudoscalarEmission --> partial_waves
    PseudoscalarEmission --> physical_state
    ThreeGluonChannel --> decay_width
    ThreeGluonChannel --> matrix_element
    TwoGluonChannel --> decay_width
    TwoGluonChannel --> matrix_element
    TwoMesonChannel --> PseudoscalarEmission
    TwoMesonChannel --> matrix_element
    TwoPhotonAnnihilation --> AnnihilationTerm
    TwoPhotonAnnihilation --> TwoPhotonChannel
    TwoPhotonAnnihilation --> decay_width
    TwoPhotonAnnihilation --> mass_correction_factor
    TwoPhotonAnnihilation --> matrix_element
    TwoPhotonChannel --> decay_width
    TwoPhotonChannel --> matrix_element
    Vacuum --> decay_width
    Vacuum --> matrix_element
    charge_radius_squared --> matrix_element
    decay_width --> GluonicAnnihilation
    decay_width --> LeptonNeutrinoChannel
    decay_width --> LeptonicCurrent
    decay_width --> MasslessLeptonPair
    decay_width --> PhotonEmission
    decay_width --> PhysicalState
    decay_width --> PseudoscalarEmission
    decay_width --> ThreeGluonChannel
    decay_width --> TwoGluonChannel
    decay_width --> TwoMesonChannel
    decay_width --> TwoPhotonAnnihilation
    decay_width --> TwoPhotonChannel
    decay_width --> mass_correction_factor
    decay_width --> matrix_element
    decay_width --> partial_waves
    decay_width --> physical_state
    mass_correction_factor --> GluonicAnnihilation
    mass_correction_factor --> LeptonicCurrent
    mass_correction_factor --> PartialWave
    mass_correction_factor --> PhotonEmission
    mass_correction_factor --> PhysicalState
    mass_correction_factor --> PseudoscalarEmission
    mass_correction_factor --> TwoPhotonAnnihilation
    mass_correction_factor --> decay_width
    mass_correction_factor --> matrix_element
    matrix_element --> GluonicAnnihilation
    matrix_element --> LeptonicCurrent
    matrix_element --> PartialWave
    matrix_element --> PhotonEmission
    matrix_element --> PhysicalState
    matrix_element --> PseudoscalarEmission
    matrix_element --> TwoMesonChannel
    matrix_element --> TwoPhotonAnnihilation
    matrix_element --> TwoPhotonChannel
    matrix_element --> Vacuum
    matrix_element --> charge_radius_squared
    matrix_element --> decay_width
    matrix_element --> mass_correction_factor
    matrix_element --> partial_waves
    matrix_element --> physical_state
    partial_waves --> PartialWave
    partial_waves --> decay_width
    partial_waves --> matrix_element
    physical_state --> PhysicalState
    physical_state --> matrix_element
```

</details>
