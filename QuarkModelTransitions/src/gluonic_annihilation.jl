"""
`GluonicAnnihilation(masses, alpha_s)`: integrated lowest-order gluonic width.
`alpha_s` is either a fixed number or a callable `alpha_s(mass_GeV)`, evaluated
at the initial state's mass. `mass_correction_factor` requires `target_mass`
in GeV and multiplies the **width**, reevaluating a callable coupling there.
A fixed coupling stays fixed. Coherent mixed gluonic
annihilation is not provided by this prescription; use a single pure component.

Normalization: `decay_width` returns MeV from the published integrated
lowest-order rate, including phase space and spin/color/symmetry factors.
No differential or helicity matrix element is exposed.

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
operator = GluonicAnnihilation(masses, M -> GIModel.alpha_s_q(M))
@assert decay_width(TwoGluonChannel(), operator, final) > 0
```

## Related

[`TwoGluonChannel`](@ref), [`ThreeGluonChannel`](@ref), [`mass_correction_factor`](@ref), [`decay_width`](@ref).
"""
struct GluonicAnnihilation{A} <: AnnihilationOperator
    masses::QuarkMassTable
    alpha_s::A
    function GluonicAnnihilation(masses::QuarkMassTable, alpha_s)
        if alpha_s isa Real
            isfinite(alpha_s) && alpha_s >= 0 || throw(ArgumentError("invalid alpha_s"))
            alpha_s = Float64(alpha_s)
        else
            applicable(alpha_s, 1.0) || throw(ArgumentError("alpha_s must be a number or callable at a mass in GeV"))
        end
        all(isfinite(m) && m > 0 for m in values(masses)) || throw(ArgumentError("invalid masses"))
        new{typeof(alpha_s)}(copy(masses), alpha_s)
    end
end

function _gluonic_alpha_s(op::GluonicAnnihilation, mass::Real)
    a = op.alpha_s isa Real ? op.alpha_s : op.alpha_s(mass)
    a isa Real && isfinite(a) && a >= 0 || throw(ArgumentError("alpha_s must return a finite nonnegative real coupling"))
    return a
end

function mass_correction_factor(final::Union{TwoGluonChannel,ThreeGluonChannel},
                                op::GluonicAnnihilation, initial::PhysicalState;
                                target_mass::Real, npoints::Integer=900)
    target = PhysicalState(initial.label, _correction_mass(target_mass), initial.components)
    return _correction_ratio(decay_width(final, op, target; npoints),
                             decay_width(final, op, initial; npoints))
end

function decay_width(final::Union{TwoGluonChannel,ThreeGluonChannel}, op::GluonicAnnihilation,
                     initial::PhysicalState; npoints::Integer=900)
    length(initial.components) == 1 || throw(ArgumentError("mixed gluonic annihilation is not implemented"))
    c = only(initial.components)
    a, b = c.basis.flavors
    a == b || throw(ArgumentError("gluonic annihilation requires flavor-diagonal constituents"))
    quantum = _ann_quantum(c)
    channel = final isa ThreeGluonChannel ?
        (quantum == (0,3,1) ? :S1_3g : nothing) :
        get(Dict((0,1,0)=>:S0_2g, (1,3,0)=>:P0_2g, (1,3,2)=>:P2_2g), quantum, nothing)
    isnothing(channel) && throw(ArgumentError("unsupported gluon channel for $(c.basis.label)"))
    m = _ann_mass(op,a)
    S = wavefunction_origin_smearing(c.wave, m; L=quantum[1], npoints)
    return GEV_TO_MEV * abs2(c.coefficient) * gluonic_annihilation_width(
        channel, S, _gluonic_alpha_s(op, initial.mass_GeV), m)
end

