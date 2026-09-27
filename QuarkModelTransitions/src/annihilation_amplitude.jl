"""
Reduced annihilation result with explicit normalization and component provenance.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = PhotonEmission(masses)
result = matrix_element(Vacuum(), LeptonicCurrent(:electromagnetic, masses), initial)
@assert result isa QuarkModelTransitions.AnnihilationAmplitude
@assert result.provenance.amplitude_units == :dimensionless
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct AnnihilationAmplitude{O,F,I,T,P}
    operator::O
    final::F
    initial::I
    value::ComplexF64
    terms::T
    provenance::P
end

function _ann_mass(op, f)
    key = String(f == :n ? :q : f)
    haskey(op.masses, key) || throw(ArgumentError("missing constituent mass for $f"))
    return op.masses[key]
end
_ann_quantum(c) = (orbital_angular_momentum(c.basis.L_label), c.basis.multiplicity, c.basis.J)

function _explicit_ann_flavor_coefficient(op, c)
    index = findfirst(t -> t.flavors == c.basis.flavors, op.terms)
    isnothing(index) && throw(ArgumentError("missing flavor coefficient for $(c.basis.flavors)"))
    return op.terms[index].coefficient
end

_ann_components(op, initial) = initial.components

function _ann_compose(final, op, initial, npoints)
    npoints > 0 || throw(ArgumentError("npoints must be positive"))
    terms = Tuple(begin
        flavor = _ann_flavor_coefficient(op, c)
        kernel = _ann_kernel(final, op, c, initial.mass_GeV, npoints)
        (basis = c.basis, mixing = c.coefficient, flavor_coefficient = flavor,
         prefactor = c.coefficient * flavor,
         kernel = kernel, contribution = c.coefficient * flavor * kernel)
    end for c in _ann_components(op, initial))
    photon = final isa TwoPhotonChannel
    return AnnihilationAmplitude(op, final, initial,
        ComplexF64(sum(t.contribution for t in terms)), terms,
        (normalization = photon ? :GI_width_amplitude : :GI_reduced_current,
         amplitude_units = photon ? :GeV_sqrt : :dimensionless,
         width_units = :MeV, npoints = npoints,
         source = :mock_meson_annihilation, helicity_available = false))
end
