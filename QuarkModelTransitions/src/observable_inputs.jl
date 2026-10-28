"""
Published phenomenological defaults used by electromagnetic observables.

These are additional transition-model inputs, not spectrum parameters or
numerical tolerances. Individual functions accept keyword overrides.

`electric_exponent` is the power p of the smearing |m/sqrt(<E><E>)|^p on the
electric radial moments E_n. GI 1985 fitted p = 0.5 to A2 -> pi gamma, with
spin-flip denominators that were sqrt(2) too small. Measured E1 widths prefer
p = 0, the choice of Godfrey's later GI-model papers, and no single p fits the
light mesons (GIPaper report `radiative_exponent.md`), so the default is 0.
Pass `electric_exponent = 0.5` to reproduce the paper. `magnetic_exponent`
(0.7, fitted to rho -> pi gamma) is the paper's value.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
operator = PhotonEmission(QuarkMassTable("c" => 1.628))
@assert operator.magnetic_exponent == QMT.ELECTROMAGNETIC_DEFAULTS.magnetic_exponent
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
const ELECTROMAGNETIC_DEFAULTS = (
    magnetic_exponent = 0.7,
    electric_exponent = 0.0,
    charge_radius_exponent = 0.2,
    recoil_beta_GeV = 0.40,
)

"""
Calibration defaults for the frozen Table IV/V model (amplitudes in MeV^1/2).

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert QMT.STRONG_DECAY_DEFAULTS.beta_GeV == 0.4
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
const STRONG_DECAY_DEFAULTS = (
    rho_amplitude = 12.4,
    B_amplitude = -11.0,
    beta_GeV = 0.40,
)

"""
Number of MeV in one GeV; used by the generic width interface.

## Example

```julia
import GIModel.QuarkModelTransitions as QMT
@assert 0.001 * QMT.GEV_TO_MEV == 1.0
```

## Related

[`decay_width`](@ref).
"""
const GEV_TO_MEV = 1000.0
