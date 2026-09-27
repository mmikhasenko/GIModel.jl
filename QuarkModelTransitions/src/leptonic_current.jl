"""
`LeptonicCurrent(kind, masses, terms)`: GI reduced current, including flavor factors.
Kinds are :P_P, :V_V, :Vp_V, :Pp_A1. Returns a dimensionless current factor,
not a physical width; in particular the pseudoscalar convention is f_P/M_P.

`LeptonicCurrent(:electromagnetic, masses)` infers V or V′ from each component's
³S₁ or ³D₁ assignment and calculates its coefficient with
`vector_current_prefactor`. Supply explicit uū/dd̄ components for rho/omega;
the mass-sector aliases q/n do not specify a charge or an isospin state.
States obtained from an isoscalar spectrum carry that information explicitly;
the electromagnetic current resolves their nonstrange channel as (uū+dd̄)/√2.
Result terms expose `prefactor` (including state mixing) and `kernel` (V or V′),
whose product is `contribution`. Different components can have different kernels;
in general there is no single factorizable V for a mixed state.

Normalization: flavor-weighted dimensionless GI current factors of Appendix D.
They contain mock-meson mass and relativistic weights, but no leptonic phase
space. `decay_width` uses D7 for `LeptonNeutrinoChannel` (explicit CKM magnitude)
or D8 for `MasslessLeptonPair`, returning MeV. Lepton masses are neglected in D8.

`mass_correction_factor(Vacuum(), operator, initial; target_mass)` requires a
positive mass in GeV and multiplies `result.value` at fixed wavefunctions.
It excludes the width's separate mass and leptonic phase-space factors.

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
operator = LeptonicCurrent(:electromagnetic, masses)
@assert decay_width(MasslessLeptonPair(), operator, initial) > 0
```

## Related

[`Vacuum`](@ref), [`LeptonNeutrinoChannel`](@ref), [`MasslessLeptonPair`](@ref), [`AnnihilationTerm`](@ref), `vector_current_prefactor`, [`mass_correction_factor`](@ref), [`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct LeptonicCurrent{T} <: AnnihilationOperator
    kind::Symbol
    masses::QuarkMassTable
    terms::T
    function LeptonicCurrent(kind::Symbol, masses::QuarkMassTable, terms)
        if kind == :electromagnetic
            (terms isa Tuple || terms isa AbstractVector) && isempty(terms) ||
                throw(ArgumentError("electromagnetic coefficients are inferred from the state"))
            all(isfinite(m) && m > 0 for m in values(masses)) || throw(ArgumentError("invalid masses"))
            return new{Tuple{}}(kind, copy(masses), ())
        end
        haskey(LEPTONIC_FACTOR_KINDS, kind) || throw(ArgumentError("unsupported current $kind"))
        ms, ts = _annihilation_inputs(masses, terms)
        new{typeof(ts)}(kind, ms, ts)
    end
end

LeptonicCurrent(kind::Symbol, masses::QuarkMassTable) = LeptonicCurrent(kind, masses, ())

matrix_element(final::Vacuum, op::LeptonicCurrent, initial::PhysicalState; npoints::Integer=900) =
    _ann_compose(final, op, initial, npoints)


function mass_correction_factor(final::Vacuum, op::LeptonicCurrent,
                                initial::PhysicalState; target_mass::Real,
                                npoints::Integer=900)
    target = PhysicalState(initial.label, _correction_mass(target_mass), initial.components;
                           provenance = initial.provenance)
    return _correction_ratio(matrix_element(final, op, target; npoints).value,
                             matrix_element(final, op, initial; npoints).value)
end

function decay_width(final::LeptonNeutrinoChannel, op::LeptonicCurrent,
                     initial::PhysicalState; G::Real=G_FERMI_GEV, npoints::Integer=900)
    op.kind == :P_P || throw(ArgumentError("lepton-neutrino decay requires a pseudoscalar current"))
    a = matrix_element(Vacuum(), op, initial; npoints)
    return GEV_TO_MEV * final.ckm^2 * leptonic_pseudoscalar_width(
        abs(a.value), initial.mass_GeV, final.lepton_mass_GeV; G)
end
function decay_width(::MasslessLeptonPair, op::LeptonicCurrent,
                     initial::PhysicalState; alpha::Real=ALPHA_EM, npoints::Integer=900)
    op.kind in (:V_V, :Vp_V, :electromagnetic) || throw(ArgumentError("dilepton decay requires a vector current"))
    a = matrix_element(Vacuum(), op, initial; npoints)
    return GEV_TO_MEV * dilepton_vector_width(abs(a.value), initial.mass_GeV; alpha)
end

# --- Flavor coefficients and radial kernel ---

"""
vector_current_prefactor(basis::BasisState)

Coefficient multiplying GI's V (³S₁) or V′ (³D₁) for one neutral flavor component.
It is `2sqrt(3) * Q_f * a_L`, where `a_0=1` and `a_2=sqrt(2)/3`.
The orbital factors are the normalized spin/angular reduction of the vector
current in GI Appendix D, with its distinct definitions of V and V′; they are
not fitted or radial-level dependent. This is not a new symbolic derivation.
State mixing coefficients multiply this result separately. Flavor aliases q/n
are ambiguous and rejected. Charges are in units of e (no extra antiquark charge).

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
basis = BasisState(1, "S", 3, 1; flavors=(:c, :c))
@assert QMT.vector_current_prefactor(basis) ≈ sqrt(16/3)
```

## Related

[`LeptonicCurrent`](@ref).
"""
function vector_current_prefactor(basis::BasisState)
    L = orbital_angular_momentum(basis.L_label)
    basis.multiplicity == 3 && basis.J == 1 && L in (0,2) ||
        throw(ArgumentError("electromagnetic vector current requires ³S₁ or ³D₁"))
    f, antif = basis.flavors
    f == antif || throw(ArgumentError("electromagnetic annihilation requires flavor-diagonal constituents"))
    charge = f in (:u,:c,:t) ? 2/3 : f in (:d,:s,:b) ? -1/3 :
        throw(ArgumentError("specify explicit u, d, s, c, b or t flavor, not $f"))
    return 2sqrt(3) * charge * (L == 0 ? 1 : sqrt(2)/3)
end

_ann_flavor_coefficient(op::LeptonicCurrent, c) =
    op.kind == :electromagnetic ? vector_current_prefactor(c.basis) :
    _explicit_ann_flavor_coefficient(op, c)

function _ann_kernel(::Vacuum, op::LeptonicCurrent, c, M, npoints)
    kind = op.kind == :electromagnetic ? (_ann_quantum(c)[1] == 0 ? :V_V : :Vp_V) : op.kind
    expected = Dict(:P_P => (0,1,0), :V_V => (0,3,1), :Vp_V => (2,3,1), :Pp_A1 => (1,3,1))[kind]
    _ann_quantum(c) == expected || throw(ArgumentError("$(op.kind) does not support $(c.basis.label)"))
    a, b = c.basis.flavors
    return leptonic_decay_factor(kind, c.wave, _ann_mass(op,a), _ann_mass(op,b), M; npoints)
end


_ann_components(op::LeptonicCurrent, initial::PhysicalState) =
    op.kind == :electromagnetic ? _electromagnetic_components(initial) : initial.components
