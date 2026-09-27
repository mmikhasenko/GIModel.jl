"""
    mass_correction_factor(final, operator, initial; target_mass, kwargs...)
    mass_correction_factor(final, operator, initial; target_momentum, kwargs...)

Compare another external mass/momentum with the default input-state kinematics,
holding wavefunctions, mixing, constituent masses and operator parameters fixed.
This is not a new spectrum calculation. The operator's docstring specifies the
required keyword; unsupported keywords are errors, not silently ignored.

| Operator | Required keyword | Multiplies |
|:--|:--|:--|
| [`LeptonicCurrent`](@ref), [`TwoPhotonAnnihilation`](@ref) | `target_mass` (GeV) | Numerical current/amplitude |
| [`PhotonEmission`](@ref) | `target_momentum` (GeV) | `result.value` |
| [`PseudoscalarEmission`](@ref) | `target_momentum` (GeV), `partial_wave` | `result[partial_wave]` |
| [`GluonicAnnihilation`](@ref) | `target_mass` (GeV) | Width, since no matrix element is defined |

Emission corrections may reevaluate momentum-dependent integrals. A ratio is
undefined when the reference amplitude vanishes: those cases throw `DomainError`.
Closed reference emission channels also throw; a factor cannot open a channel
whose reference amplitude is unavailable. Factors can be signed or complex.
Do not square a current/amplitude factor and assume it is a width correction:
additional phase-space and mass factors may enter the width separately.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = TwoPhotonAnnihilation(masses, AnnihilationTerm((:c, :c), 4/9))
amplitude = matrix_element(TwoPhotonChannel(), operator, final)
factor = mass_correction_factor(TwoPhotonChannel(), operator, final; target_mass=2.98)
comparison_amplitude = factor * amplitude.value
@assert factor ≈ (2.98/final.mass_GeV)^1.5
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref), [`PhysicalState`](@ref),
[`PartialWave`](@ref), `OnShell`, `CMKinematics`.
"""
function mass_correction_factor end

function _correction_mass(m::Real)
    isfinite(m) && m > 0 || throw(ArgumentError("target_mass must be finite and positive (GeV)"))
    return float(m)
end

function _correction_ratio(target::Number, reference::Number)
    isfinite(target) && isfinite(reference) ||
        throw(DomainError((target, reference), "correction requires finite amplitudes"))
    iszero(reference) && throw(DomainError(reference,
        "a multiplicative correction is undefined for a zero reference amplitude"))
    ratio = target / reference
    isfinite(ratio) || throw(DomainError(ratio, "correction ratio overflowed"))
    return ratio
end
