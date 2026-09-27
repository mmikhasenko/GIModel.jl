"""
`TwoPhotonAnnihilation(masses, terms)`: 1S0 or 3P2 annihilation.
Coefficients are effective squared quark charges. No helicity decomposition
is inferred from the published integrated formula.

Normalization: signed GI width amplitude `a` in `sqrt(GeV)`, with
`Gamma[MeV] = GEV_TO_MEV * abs2(a)`. Integrated phase space, spin averaging,
and identical-photon symmetry are already included; do not multiply by them again.

`mass_correction_factor(TwoPhotonChannel(), operator, initial; target_mass)`
requires a positive mass in GeV and multiplies `result.value` at fixed waves.
Its absolute square therefore also corrects the integrated two-photon width.

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
result = matrix_element(TwoPhotonChannel(), operator, final)
@assert decay_width(result) ≈ 1000abs2(result.value)
```

## Related

[`AnnihilationTerm`](@ref), [`TwoPhotonChannel`](@ref), `AnnihilationAmplitude`, [`mass_correction_factor`](@ref), [`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct TwoPhotonAnnihilation{T} <: AnnihilationOperator
    masses::QuarkMassTable
    terms::T
    function TwoPhotonAnnihilation(masses::QuarkMassTable, terms)
        ms, ts = _annihilation_inputs(masses, terms)
        new{typeof(ts)}(ms, ts)
    end
end

matrix_element(final::TwoPhotonChannel, op::TwoPhotonAnnihilation, initial::PhysicalState; npoints::Integer=900) =
    _ann_compose(final, op, initial, npoints)


function mass_correction_factor(final::TwoPhotonChannel, op::TwoPhotonAnnihilation,
                                initial::PhysicalState; target_mass::Real,
                                npoints::Integer=900)
    target = PhysicalState(initial.label, _correction_mass(target_mass), initial.components)
    return _correction_ratio(matrix_element(final, op, target; npoints).value,
                             matrix_element(final, op, initial; npoints).value)
end


decay_width(a::AnnihilationAmplitude{<:TwoPhotonAnnihilation}) = GEV_TO_MEV * abs2(a.value)
decay_width(final::TwoPhotonChannel, op::TwoPhotonAnnihilation, initial::PhysicalState; kwargs...) =
    decay_width(matrix_element(final, op, initial; kwargs...))

# --- Flavor coefficient and radial kernel ---

_ann_flavor_coefficient(op::TwoPhotonAnnihilation, c) = _explicit_ann_flavor_coefficient(op, c)

function _ann_kernel(::TwoPhotonChannel, op::TwoPhotonAnnihilation, c, M, npoints)
    kind = _ann_quantum(c) == (0,1,0) ? :P : _ann_quantum(c) == (1,3,2) ? :P2 :
        throw(ArgumentError("two-photon kernel supports only 1S0 and 3P2"))
    a, b = c.basis.flavors
    a == b || throw(ArgumentError("two-photon annihilation requires flavor-diagonal components"))
    return two_photon_amplitude(kind, c.wave, _ann_mass(op,a), M, 1.0; npoints)
end

