# =============================================================================
# Transition-amplitude domain model
# =============================================================================

"""
Abstract supertype of operators evaluated by `matrix_element`. Implementations determine the final-state type and normalization.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert PhotonEmission <: QMT.TransitionOperator
```

## Related

- [`superpose`](@ref) — construct coherent combinations.
- [`matrix_element`](@ref) — evaluate an operator between states.
"""
abstract type TransitionOperator end
abstract type AnnihilationOperator <: TransitionOperator end
"""
Abstract transition operator for an ordered two-meson final channel.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
@assert PseudoscalarEmission <: QMT.StrongDecayOperator
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
"""
abstract type StrongDecayOperator <: TransitionOperator end
abstract type TransitionState end

"""One pure basis contribution to a resolved physical state."""
struct StateComponent{B,C<:Number,W<:RadialWave}
    basis::B
    coefficient::C
    wave::W
end

"""
    PhysicalState(label, mass_GeV, components; provenance=(source=:manual,))
    physical_state(spectrum, state_or_label)

An eagerly resolved external meson state. Its physical mass, common `J^P`,
pure-basis components, mixing coefficients, native radial waves, and provenance
are independent of the [`Spectrum`](@ref) that produced it. Prefer
[`physical_state`](@ref) to obtain these inputs from a GIModel calculation.
The manual constructor also accepts assumed radial waves and supplied masses;
constructing this object does not solve a Hamiltonian or certify those inputs.

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
@assert initial.J == 1
@assert length(initial.components) == 1
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
- [`physical_state`](@ref) — adapt a solved spectrum state.
- [`superpose`](@ref) — construct coherent combinations.
- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
- [`TwoMesonChannel`](@ref) — ordered surviving and emitted daughters.
"""
struct PhysicalState{C<:Tuple,P} <: TransitionState
    label::String
    mass_GeV::Float64
    J::Int
    parity::Int
    components::C
    provenance::P
end

_meson_parity(basis::BasisState) = iseven(orbital_angular_momentum(basis.L_label)) ? -1 : 1

function _state_component(component)
    basis = component.basis
    basis isa BasisState || throw(ArgumentError(
        "physical-state components must carry a BasisState, got $(typeof(basis))",
    ))
    coefficient = component.coefficient
    coefficient isa Number || throw(ArgumentError(
        "physical-state coefficients must be numeric, got $(typeof(coefficient))",
    ))
    isfinite(coefficient) || throw(ArgumentError(
        "physical-state coefficients must be finite, got $coefficient",
    ))
    wave = component.wave
    wave isa RadialWave || throw(ArgumentError(
        "physical-state components must carry a RadialWave, got $(typeof(wave))",
    ))
    return StateComponent(basis, coefficient, wave)
end

function PhysicalState(
    label::AbstractString,
    mass_GeV::Real,
    components;
    provenance = (source = :manual,),
)
    mass = Float64(mass_GeV)
    isfinite(mass) && mass > 0 || throw(ArgumentError(
        "physical-state mass must be finite and positive, got $mass_GeV",
    ))
    resolved = Tuple(_state_component(component) for component in components)
    isempty(resolved) && throw(ArgumentError("a PhysicalState needs at least one component"))
    J = first(resolved).basis.J
    parity = _meson_parity(first(resolved).basis)
    all(component.basis.J == J for component in resolved) || throw(ArgumentError(
        "all physical-state components must have the same total J",
    ))
    all(_meson_parity(component.basis) == parity for component in resolved) ||
        throw(ArgumentError("all physical-state components must have the same parity"))
    sum(abs2(component.coefficient) for component in resolved) > 0 ||
        throw(ArgumentError("physical-state component coefficients cannot all vanish"))
    return PhysicalState{typeof(resolved),typeof(provenance)}(
        String(label), mass, J, parity, resolved, provenance,
    )
end

# Resolve only a channel explicitly identified by the isoscalar solver. Keep
# the original state's coarse components for operators whose AnnihilationTerm
# already supplies effective flavor coefficients.
function _electromagnetic_components(state::PhysicalState)
    (hasproperty(state.provenance, :nonstrange_isoscalar) &&
     state.provenance.nonstrange_isoscalar) || return state.components
    return Tuple(Iterators.flatten(map(state.components) do c
        basis = c.basis
        if basis.flavors in ((:q, :q), (:u, :u), (:d, :d))
            return Tuple(StateComponent(
                BasisState(basis.n, basis.L_label, basis.multiplicity, basis.J;
                           label = basis.label, flavors = (flavor, flavor)),
                c.coefficient / sqrt(2), c.wave,
            ) for flavor in (:u, :d))
        end
        return (c,)
    end))
end

function _components_norm_squared(components)
    value = 0.0 + 0.0im
    for a in components, b in components
        ba, bb = a.basis, b.basis
        (ba.flavors, ba.L_label, ba.multiplicity, ba.J) ==
            (bb.flavors, bb.L_label, bb.multiplicity, bb.J) || continue
        value += conj(a.coefficient) * b.coefficient * radial_overlap(a.wave, b.wave, _ -> 1.0)
    end
    result = real(value)
    isfinite(result) && result > 0 || throw(ArgumentError("state combination has zero or invalid norm"))
    return result
end

"""
    superpose(states, coefficients; mass_GeV, label="superposition", normalize=true)

Form a coherent linear combination of [`PhysicalState`](@ref)s. Supply the
physical mass explicitly: a superposition of states with different masses is
not automatically a mass eigenstate. Normalization includes radial overlaps
and interference, not just the sum of coefficient squares. Known isoscalar
light channels are resolved into explicit uū/dd̄ components before combining.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
wave = OscillatorWave(0, 0.5, [1.0])
up = PhysicalState("uubar", 0.14, [(basis=BasisState(1,"S",1,0; flavors=(:u,:u)), coefficient=1.0, wave)])
down = PhysicalState("ddbar", 0.14, [(basis=BasisState(1,"S",1,0; flavors=(:d,:d)), coefficient=1.0, wave)])
pi0 = superpose([up, down], [1, -1]; label="pi0", mass_GeV=0.135)
@assert pi0.components[2].coefficient ≈ -1/sqrt(2)
```

## Related

[`PhysicalState`](@ref), [`physical_state`](@ref).
"""
function superpose(states, coefficients; mass_GeV::Real,
                   label::AbstractString = "superposition", normalize::Bool = true)
    length(states) == length(coefficients) && !isempty(states) ||
        throw(ArgumentError("supply one coefficient per state and at least one state"))
    all(state -> state isa PhysicalState, states) || throw(ArgumentError("superpose requires PhysicalStates"))
    all(c -> c isa Number && isfinite(c), coefficients) || throw(ArgumentError("coefficients must be finite numbers"))
    components = Tuple(StateComponent(c.basis, weight * c.coefficient, c.wave)
                       for (state, weight) in zip(states, coefficients)
                       for c in _electromagnetic_components(state))
    # Validate quantum numbers before evaluating overlaps.
    candidate = PhysicalState(label, mass_GeV, components)
    norm_squared = _components_norm_squared(candidate.components)
    scale = normalize ? sqrt(norm_squared) : 1.0
    return PhysicalState(label, mass_GeV,
        (StateComponent(c.basis, c.coefficient / scale, c.wave) for c in components);
        provenance = (source = :superposition, input_norm_squared = norm_squared,
                      normalized = normalize))
end


_state_mass(state::CentralState) = state.central_GeV
_state_mass(state::Union{CorrectedState,MixedState}) = state.mass_GeV

function _physical_state_provenance(spec::Spectrum, state)
    channels = Tuple((meson.flavor1, meson.flavor2) for meson in spec.channels)
    return (
        source = :spectrum,
        spectrum_stage = nameof(typeof(state)),
        channels = channels,
        state_label = state.label,
        nonstrange_isoscalar = spec.nonstrange_isoscalar,
    )
end

"""
    physical_state(spectrum, state_or_label)

Resolve a solver state, `BasisState`, or label into a `PhysicalState`, retaining its mass in GeV, coherent components, radial waves, and provenance.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
basis = BasisState(1, "S", 3, 1; flavors=(:c, :c))
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels=[basis], solver=FiniteDifferenceSolver(ngrid=90, rmax=12.0, nlevels_per_channel=1))
state = physical_state(spectrum, basis)
@assert state.mass_GeV > 0
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
- [`PhysicalState`](@ref) — resolved meson components and mass.
"""
function physical_state(
    spec::Spectrum,
    state::Union{CentralState,CorrectedState,MixedState},
)
    components = physical_components(spec, state)
    return PhysicalState(
        state.label,
        _state_mass(state),
        components;
        provenance = _physical_state_provenance(spec, state),
    )
end

physical_state(spec::Spectrum, label::AbstractString) =
    physical_state(spec, spectrum_state(spec, label))

physical_state(spec::Spectrum, basis::BasisState) =
    physical_state(spec, spectrum_state(spec, basis))

"""
    ReferenceState(label, mass_GeV; J=nothing, parity=nothing, provenance=...)

A mass-and-label external state for a reference backend. It deliberately has no
radial wave. Supply both `J` and parity (`+1` or `-1`) when angular selection
rules are required; omit both for legacy rows whose quantum numbers are not yet
encoded. Solver-native operators must reject this state kind.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
initial = QMT.ReferenceState("rho", 0.77; J=1, parity=-1)
final = TwoMesonChannel(QMT.ReferenceState("pi+", 0.14; J=0, parity=-1), QMT.ReferenceState("pi-", 0.14; J=0, parity=-1))
@assert initial.mass_GeV == 0.77
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
struct ReferenceState{P} <: TransitionState
    label::String
    mass_GeV::Float64
    J::Union{Nothing,Int}
    parity::Union{Nothing,Int}
    provenance::P
end

function ReferenceState(
    label::AbstractString,
    mass_GeV::Real;
    J::Union{Nothing,Integer} = nothing,
    parity::Union{Nothing,Integer} = nothing,
    provenance = (source = :reference,),
)
    mass = Float64(mass_GeV)
    isfinite(mass) && mass > 0 || throw(ArgumentError(
        "reference-state mass must be finite and positive, got $mass_GeV",
    ))
    isnothing(J) == isnothing(parity) || throw(ArgumentError(
        "ReferenceState requires both J and parity, or neither",
    ))
    j = isnothing(J) ? nothing : Int(J)
    p = isnothing(parity) ? nothing : Int(parity)
    isnothing(j) || j >= 0 || throw(ArgumentError("J must be non-negative"))
    isnothing(p) || p in (-1, 1) || throw(ArgumentError("parity must be +1 or -1"))
    return ReferenceState{typeof(provenance)}(String(label), mass, j, p, provenance)
end

_basis_identity(basis::BasisState) = (
    basis.n,
    basis.L_label,
    basis.multiplicity,
    basis.J,
    basis.flavors,
)

function _physical_ray(state::PhysicalState)
    ordered = sort(
        collect(state.components);
        by = component -> repr(_basis_identity(component.basis)),
    )
    anchor = findfirst(component -> !iszero(component.coefficient), ordered)
    isnothing(anchor) && error("PhysicalState invariant violated: all coefficients vanish")
    anchor_coefficient = ordered[anchor].coefficient
    return Tuple(
        (_basis_identity(component.basis), component.coefficient / anchor_coefficient) for
        component in ordered
    )
end

function _state_identity(state::PhysicalState)
    return (state.label, state.mass_GeV, state.J, state.parity, _physical_ray(state))
end

_state_identity(state::ReferenceState) =
    (state.label, state.mass_GeV, state.J, state.parity, :reference)

_same_external_state(a::TransitionState, b::TransitionState) =
    _state_identity(a) == _state_identity(b)

"""
    TwoMesonChannel(first, second)

An ordered pair of daughter states. The constructor preserves the supplied
order; it never sorts the daughters.

For [`PseudoscalarEmission`](@ref), use
`TwoMesonChannel(surviving, emitted)`: `second` must be the emitted `0^-`
state. Reverse the arguments to choose the other assignment in a
two-pseudoscalar channel. Identical-particle rules use state identity, not
sorting.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
a1, vector, pion_basis = BasisState(1, "P", 3, 1), BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
charged = compute_spectrum(params, Meson(masses, :u, :d);
    levels=[a1, vector, pion_basis], solver)
initial = physical_state(charged, a1)
rho = physical_state(charged, vector)
# Resolve pi0 = (u ubar - d dbar)/sqrt(2) in the isospin limit.
up = compute_spectrum(params, Meson(masses, :u, :u); levels=[pion_basis], solver)
down = compute_spectrum(params, Meson(masses, :d, :d); levels=[pion_basis], solver)
u, d = physical_state(up, pion_basis), physical_state(down, pion_basis)
pion = PhysicalState("pi0", (u.mass_GeV + d.mass_GeV)/2,
    [(basis=c.basis, coefficient=sign*c.coefficient/sqrt(2), wave=c.wave)
     for (state, sign) in ((u, 1), (d, -1)) for c in state.components];
    provenance=(source=:isospin_combination,))
final = TwoMesonChannel(rho, pion)
g, h = 0.7, 0.3 # Illustrative GeV^-1 inputs, not fitted GI predictions.
operator = PseudoscalarEmission(g, h, masses)
@assert final.second.J == 0
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
"""
struct TwoMesonChannel{A<:TransitionState,B<:TransitionState}
    first::A
    second::B
    function TwoMesonChannel(first::A, second::B) where {
        A<:TransitionState,
        B<:TransitionState,
    }
        return new{A,B}(first, second)
    end
end

"""
Relative orbital angular momentum `L` and coupled daughter spin `S`.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
wave = PartialWave(1, 0)
@assert wave.relative_L == 1
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
"""
struct PartialWave
    relative_L::Int
    channel_spin::Int
    function PartialWave(relative_L::Integer, channel_spin::Integer)
        relative_L >= 0 || throw(ArgumentError("relative L must be non-negative"))
        channel_spin >= 0 || throw(ArgumentError("channel spin must be non-negative"))
        new(Int(relative_L), Int(channel_spin))
    end
end

"""Internal exchange rule derived from the resolved daughter states and `(L,S)`."""
struct _DaughterSymmetry{T<:Real}
    identical::Bool
    exchange_phase::Int
    normalization::T
end

function _daughter_symmetry(final::TwoMesonChannel, wave::PartialWave)
    identical = _same_external_state(final.first, final.second)
    exchange_phase = isodd(wave.relative_L + wave.channel_spin) ? -1 : 1
    normalization = identical ? inv(sqrt(2.0)) : 1.0
    return _DaughterSymmetry(identical, exchange_phase, normalization)
end

Base.isless(a::PartialWave, b::PartialWave) =
    (a.relative_L, a.channel_spin) < (b.relative_L, b.channel_spin)

function _state_quantum_numbers(state::TransitionState)
    isnothing(state.J) && throw(ArgumentError(
        "state `$(state.label)` has no J^P metadata; supply it before deriving partial waves",
    ))
    isnothing(state.parity) && throw(ArgumentError(
        "state `$(state.label)` has no J^P metadata; supply it before deriving partial waves",
    ))
    return state.J, state.parity
end

"""
    allowed_partial_waves(final, initial)

Derive every `(L,S)` allowed by angular momentum, parity, and identical-boson
exchange symmetry. The external states must carry `J^P` metadata.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
initial = QMT.ReferenceState("rho", 0.77; J=1, parity=-1)
final = TwoMesonChannel(QMT.ReferenceState("pi+", 0.14; J=0, parity=-1), QMT.ReferenceState("pi-", 0.14; J=0, parity=-1))
@assert QuarkModelTransitions.allowed_partial_waves(final, initial) == [PartialWave(1, 0)]
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
"""
function allowed_partial_waves(final::TwoMesonChannel, initial::TransitionState)
    J, parent_parity = _state_quantum_numbers(initial)
    j1, parity1 = _state_quantum_numbers(final.first)
    j2, parity2 = _state_quantum_numbers(final.second)
    identical = _same_external_state(final.first, final.second)
    waves = PartialWave[]
    for S in abs(j1 - j2):(j1 + j2)
        for L in abs(J - S):(J + S)
            parent_parity == parity1 * parity2 * (isodd(L) ? -1 : 1) || continue
            identical && isodd(L + S) && continue
            push!(waves, PartialWave(L, S))
        end
    end
    sort!(unique!(waves))
    return waves
end

"""Doubled daughter helicities, allowing integer and half-integer labels."""
struct HelicityLabel
    twice_lambda1::Int
    twice_lambda2::Int
end

"""A concrete linear map from an ordered helicity basis to partial waves."""
struct PartialWaveProjection{H<:Tuple,W<:Tuple,M<:AbstractMatrix}
    helicities::H
    partial_waves::W
    matrix::M
end

"""
    partial_wave_projection(final, initial)

Return the helicity-to-partial-wave map. The columns use the GI Appendix-C
reduced helicities `h_0 = H_0`, `h_m = sqrt(2) H_m` for `m > 0`. Spin-zero
daughters have the unique one-column projection; a vector plus pseudoscalar is
projected with conventional Clebsch--Gordan coefficients and reproduces Table
XI. Higher-spin compression remains gated on its own phase audit.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
initial = QMT.ReferenceState("rho", 0.77; J=1, parity=-1)
final = TwoMesonChannel(QMT.ReferenceState("pi+", 0.14; J=0, parity=-1), QMT.ReferenceState("pi-", 0.14; J=0, parity=-1))
projection = QMT.partial_wave_projection(final, initial)
@assert projection.matrix == ones(1, 1)
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
"""
function partial_wave_projection(final::TwoMesonChannel, initial::TransitionState)
    J, _ = _state_quantum_numbers(initial)
    j1, _ = _state_quantum_numbers(final.first)
    j2, _ = _state_quantum_numbers(final.second)
    waves = Tuple(allowed_partial_waves(final, initial))
    if j1 == 0 && j2 == 0
        helicities = (HelicityLabel(0, 0),)
        return PartialWaveProjection(helicities, waves, ones(Float64, length(waves), 1))
    end
    if (j1, j2) in ((1, 0), (0, 1))
        helicities = j1 == 1 ?
            (HelicityLabel(0, 0), HelicityLabel(2, 0)) :
            (HelicityLabel(0, 0), HelicityLabel(0, 2))
        matrix = Matrix{Float64}(undef, length(waves), 2)
        for (row, wave) in pairs(waves)
            L = wave.relative_L
            scale = sqrt((2L + 1) / (2J + 1))
            matrix[row, 1] = scale * CG(L, 0, 1, 0, J, 0)
            matrix[row, 2] = sqrt(2) * scale * CG(L, 0, 1, 1, J, 1)
        end
        return PartialWaveProjection(helicities, waves, matrix)
    end
    throw(ArgumentError(
        "helicity-to-partial-wave projection is implemented for spin-zero pairs " *
        "and vector-pseudoscalar channels; higher daughter spins require a " *
        "separate parity/helicity-compression audit",
    ))
end

# =============================================================================
# Coherent physical-state composition
# =============================================================================

"""One retained pure-component triple in a physical decay amplitude."""
struct _ComponentTransitionTerm{A,B,C,W,N,K,V}
    initial_component::A
    first_component::B
    second_component::C
    mixing_coefficient::W
    identical_normalization::N
    kernel_value::K
    contribution::V
end

"""Internal result of composing one physical ket with two physical bras."""
struct _CoherentComposition{V,T<:Tuple,S,M}
    value::V
    terms::T
    symmetry::S
    external_masses_GeV::M
end

"""
    _compose_physical_decay(kernel, final, initial, wave)

Compose a pure-basis decay kernel over one physical parent and two physical
daughters. `kernel(first_component, second_component, initial_component)` must
exclude mixing coefficients and identical-particle normalization. Every
component triple is retained; final-state coefficients are conjugated.
"""
function _compose_physical_decay(
    kernel,
    final::TwoMesonChannel{<:PhysicalState,<:PhysicalState},
    initial::PhysicalState,
    wave::PartialWave,
)
    wave in allowed_partial_waves(final, initial) || throw(ArgumentError(
        "partial wave (L,S)=($(wave.relative_L),$(wave.channel_spin)) is forbidden " *
        "for $(initial.label) -> $(final.first.label) + $(final.second.label)",
    ))
    symmetry = _daughter_symmetry(final, wave)
    terms = Tuple(
        begin
            mixing = parent.coefficient * conj(first.coefficient) * conj(second.coefficient)
            kernel_value = kernel(first, second, parent)
            kernel_value isa Number || throw(ArgumentError(
                "a pure-component decay kernel must return a Number, got $(typeof(kernel_value))",
            ))
            contribution = symmetry.normalization * mixing * kernel_value
            _ComponentTransitionTerm(
                parent,
                first,
                second,
                mixing,
                symmetry.normalization,
                kernel_value,
                contribution,
            )
        end for parent in initial.components for first in final.first.components for
        second in final.second.components
    )
    value = sum(term.contribution for term in terms)
    masses = (initial.mass_GeV, final.first.mass_GeV, final.second.mass_GeV)
    return _CoherentComposition(value, terms, symmetry, masses)
end

function _validated_momentum(momentum_GeV::Number)
    momentum = float(momentum_GeV)
    isfinite(momentum) || throw(ArgumentError("momentum must be finite, got $momentum_GeV"))
    momentum isa Real && momentum < 0 && throw(ArgumentError(
        "real momentum must be non-negative, got $momentum_GeV",
    ))
    return momentum
end

struct ClosedChannelError <: Exception
    parent::String
    parent_mass_GeV::Float64
    daughters::Tuple{String,String}
    threshold_GeV::Float64
end

function Base.showerror(io::IO, error::ClosedChannelError)
    print(
        io,
        "on-shell channel ", error.parent, " -> ", error.daughters[1], " + ",
        error.daughters[2], " is closed: M=", error.parent_mass_GeV,
        " GeV, threshold=", error.threshold_GeV,
        " GeV; no on-shell amplitude exists for these input masses",
    )
end

function _on_shell_momentum(final::TwoMesonChannel, initial::TransitionState)
    M = initial.mass_GeV
    m1 = final.first.mass_GeV
    m2 = final.second.mass_GeV
    M > m1 + m2 || throw(ClosedChannelError(
        initial.label,
        M,
        (final.first.label, final.second.label),
        m1 + m2,
    ))
    return sqrt((M^2 - (m1 + m2)^2) * (M^2 - (m1 - m2)^2)) / (2M)
end

abstract type AmplitudeNormalization end
"""
Native Eq. (C2) normalization. For strong amplitudes, the partial width in MeV is `1000q/(2π(2J+1))` times the sum of squared partial-wave values, with `q` in GeV.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
a1, vector, pion_basis = BasisState(1, "P", 3, 1), BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
charged = compute_spectrum(params, Meson(masses, :u, :d);
    levels=[a1, vector, pion_basis], solver)
initial = physical_state(charged, a1)
rho = physical_state(charged, vector)
# Resolve pi0 = (u ubar - d dbar)/sqrt(2) in the isospin limit.
up = compute_spectrum(params, Meson(masses, :u, :u); levels=[pion_basis], solver)
down = compute_spectrum(params, Meson(masses, :d, :d); levels=[pion_basis], solver)
u, d = physical_state(up, pion_basis), physical_state(down, pion_basis)
pion = PhysicalState("pi0", (u.mass_GeV + d.mass_GeV)/2,
    [(basis=c.basis, coefficient=sign*c.coefficient/sqrt(2), wave=c.wave)
     for (state, sign) in ((u, 1), (d, -1)) for c in state.components];
    provenance=(source=:isospin_combination,))
final = TwoMesonChannel(rho, pion)
g, h = 0.7, 0.3 # Illustrative GeV^-1 inputs, not fitted GI predictions.
operator = PseudoscalarEmission(g, h, masses)
amplitude = matrix_element(final, operator, initial)
@assert amplitude.normalization isa QMT.RelativisticTwoBodyNormalization
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
"""
struct RelativisticTwoBodyNormalization <: AmplitudeNormalization end
"""
Frozen Table V normalization: each amplitude has units MeV^(1/2), and its squared modulus is a partial width in MeV.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
initial = QMT.ReferenceState("rho", 0.77; J=1, parity=-1)
final = TwoMesonChannel(QMT.ReferenceState("pi+", 0.14; J=0, parity=-1), QMT.ReferenceState("pi-", 0.14; J=0, parity=-1))
operator = QMT.TableVReference(model, channel, PartialWave(1, 0))
amplitude = matrix_element(final, operator, initial)
@assert amplitude.normalization isa QMT.GITableVNormalization
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
struct GITableVNormalization <: AmplitudeNormalization end

"""One inspectable contribution to a transition amplitude."""
struct TransitionTerm{T<:Number,P}
    label::String
    coefficient::T
    reduced::T
    spatial::T
    contribution::T
    provenance::P
end

"""
    TransitionAmplitude

A complete transition result with primitive helicity values, every available
partial-wave projection, typed term decomposition, normalization, and
provenance. The computed momentum is stored directly as `momentum_GeV`.
Collections are tuples so no field erases its element type.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
a1, vector, pion_basis = BasisState(1, "P", 3, 1), BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
charged = compute_spectrum(params, Meson(masses, :u, :d);
    levels=[a1, vector, pion_basis], solver)
initial = physical_state(charged, a1)
rho = physical_state(charged, vector)
# Resolve pi0 = (u ubar - d dbar)/sqrt(2) in the isospin limit.
up = compute_spectrum(params, Meson(masses, :u, :u); levels=[pion_basis], solver)
down = compute_spectrum(params, Meson(masses, :d, :d); levels=[pion_basis], solver)
u, d = physical_state(up, pion_basis), physical_state(down, pion_basis)
pion = PhysicalState("pi0", (u.mass_GeV + d.mass_GeV)/2,
    [(basis=c.basis, coefficient=sign*c.coefficient/sqrt(2), wave=c.wave)
     for (state, sign) in ((u, 1), (d, -1)) for c in state.components];
    provenance=(source=:isospin_combination,))
final = TwoMesonChannel(rho, pion)
g, h = 0.7, 0.3 # Illustrative GeV^-1 inputs, not fitted GI predictions.
operator = PseudoscalarEmission(g, h, masses)
amplitude = matrix_element(final, operator, initial)
@assert amplitude isa QuarkModelTransitions.TransitionAmplitude
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
- [`partial_waves`](@ref) — inspect the available angular components.
- [`decay_width`](@ref) — convert a transition result to MeV.
"""
struct TransitionAmplitude{O,I,F,K,N,H<:Tuple,W<:Tuple,R<:Tuple,P}
    operator::O
    initial::I
    final::F
    momentum_GeV::K
    normalization::N
    helicity::H
    partial_wave_amplitudes::W
    terms::R
    provenance::P
end

"""
    partial_waves(amplitude)

Return the available `PartialWave` keys of a strong transition. Read each value with `amplitude[wave]`; photon results instead record a transition class and multipole.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
a1, vector, pion_basis = BasisState(1, "P", 3, 1), BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
charged = compute_spectrum(params, Meson(masses, :u, :d);
    levels=[a1, vector, pion_basis], solver)
initial = physical_state(charged, a1)
rho = physical_state(charged, vector)
# Resolve pi0 = (u ubar - d dbar)/sqrt(2) in the isospin limit.
up = compute_spectrum(params, Meson(masses, :u, :u); levels=[pion_basis], solver)
down = compute_spectrum(params, Meson(masses, :d, :d); levels=[pion_basis], solver)
u, d = physical_state(up, pion_basis), physical_state(down, pion_basis)
pion = PhysicalState("pi0", (u.mass_GeV + d.mass_GeV)/2,
    [(basis=c.basis, coefficient=sign*c.coefficient/sqrt(2), wave=c.wave)
     for (state, sign) in ((u, 1), (d, -1)) for c in state.components];
    provenance=(source=:isospin_combination,))
final = TwoMesonChannel(rho, pion)
g, h = 0.7, 0.3 # Illustrative GeV^-1 inputs, not fitted GI predictions.
operator = PseudoscalarEmission(g, h, masses)
amplitude = matrix_element(final, operator, initial)
waves = partial_waves(amplitude)
@assert waves == [PartialWave(0, 1), PartialWave(2, 1)]
@assert isfinite(amplitude[first(waves)])
```

## Related

- [`matrix_element`](@ref) — evaluate an operator between states.
- [`PartialWave`](@ref) — orbital and coupled-spin labels.
- `TransitionAmplitude` — strong result and its decomposition.
- [`decay_width`](@ref) — convert a transition result to MeV.
"""
partial_waves(amplitude::TransitionAmplitude) =
    [first(item) for item in amplitude.partial_wave_amplitudes]

function Base.getindex(amplitude::TransitionAmplitude, wave::PartialWave)
    index = findfirst(item -> first(item) == wave, amplitude.partial_wave_amplitudes)
    isnothing(index) && throw(KeyError(wave))
    return last(amplitude.partial_wave_amplitudes[index])
end

"""
    matrix_element(final, operator, initial; kwargs...)

Evaluate a transition using the masses, wavefunctions, and coherent components
of the supplied states. Obtain GIModel solutions with [`physical_state`](@ref).
The returned object records its normalization and contributions; use
[`decay_width`](@ref) for the supported physical-width conversion in MeV.

## Initial states

[`PhysicalState`](@ref) for native operators, `ReferenceState` for the
qualified `TableVReference` compatibility backend.

## Final states

[`PhysicalState`](@ref) for the daughter in photon emission,
[`TwoMesonChannel`](@ref) for pseudoscalar emission,
[`Vacuum`](@ref) for a meson-current matrix element,
[`TwoPhotonChannel`](@ref) for two-photon annihilation.
For pseudoscalar emission, order the daughters as `(surviving, emitted)`;
the second must be a resolved elementary `0^-` state. Reference channels must
match the labels of their supplied row.

## Operators

[`PhotonEmission`](@ref), [`PseudoscalarEmission`](@ref),
[`LeptonicCurrent`](@ref), [`TwoPhotonAnnihilation`](@ref),
`TableVReference` (internal reference backend).
[`GluonicAnnihilation`](@ref) supports only `decay_width`, not `matrix_element`.

`PhotonEmission` returns a `RadiativeAmplitude`: M1 values are moments
in nuclear magnetons, E1/M2 values are integrated amplitudes in MeV^(1/2).
The default constructor infers the multipole and flavor current from the states;
`verbose=true` reports the selection.

`PseudoscalarEmission` returns a `TransitionAmplitude` with dimensionless
Eq. (19) partial waves, compressed helicities, and coherent component terms.
Its `g,h` couplings are additional phenomenological inputs supplied by the caller.
`TableVReference` returns a `TransitionAmplitude` for one supplied paper-row
partial wave in MeV^(1/2).

`LeptonicCurrent` and `TwoPhotonAnnihilation` return an
`AnnihilationAmplitude`. Its `.value` is respectively a dimensionless GI
reduced current or a signed two-photon width amplitude in GeV^(1/2). There is no
universal invariant-amplitude normalization shared by these results.

## Masses and comparisons

Emission momenta are calculated from the input-state masses. There is no public
momentum override. Use [`mass_correction_factor`](@ref) separately to compare
another mass or momentum while keeping waves, mixing, and operator parameters
fixed. A closed reference emission channel has no amplitude correction ratio.
Explicit flavor components are required where a current or isospin combination
must be resolved; an averaged `:q` label alone does not specify that information.

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
amplitude = matrix_element(final, operator, initial)
@assert decay_width(amplitude) > 0
```

## Related

Results: [`decay_width`](@ref), [`mass_correction_factor`](@ref),
[`partial_waves`](@ref), `allowed_partial_waves`, [`PartialWave`](@ref),
`partial_wave_projection`.

Internal conventions: `RelativisticTwoBodyNormalization`,
`GITableVNormalization`, `TransitionOperator`,
`StrongDecayOperator`.

Primitive observables: `leptonic_decay_factor`,
`gluonic_annihilation_width`, `two_photon_amplitude`,
[`charge_radius_squared`](@ref).
"""
function matrix_element end
function partial_width end

_validate_transition(final, operator, initial) = nothing

function matrix_element(
    final::TwoMesonChannel,
    operator::StrongDecayOperator,
    initial::TransitionState,
)
    states = (initial, final.first, final.second)
    reference = findfirst(state -> state isa ReferenceState, states)
    if !isnothing(reference)
        state = states[reference]
        throw(ArgumentError(
            "$(nameof(typeof(operator))) requires resolved PhysicalState inputs; " *
            "`$(state.label)` is a wave-free ReferenceState",
        ))
    end
    throw(ArgumentError(
        "matrix_element is not implemented for operator $(nameof(typeof(operator)))",
    ))
end

"""
    decay_width(final, operator, initial; kwargs...)
    decay_width(amplitude)

Return a numerical partial width in **MeV**. The conversion includes the
operator's prescribed phase space, spin sums/averages, and symmetry factors.
Do not apply those factors a second time. Masses and emission momenta come from
the input states; a closed native emission channel returns zero.

## Initial states

[`PhysicalState`](@ref) from [`physical_state`](@ref), `ReferenceState`
for qualified `TableVReference` compatibility calculations.

## Final states

[`PhysicalState`](@ref) for photon emission, [`TwoMesonChannel`](@ref) for
pseudoscalar emission or reference rows, [`LeptonNeutrinoChannel`](@ref) for
pseudoscalar leptonic decay, [`MasslessLeptonPair`](@ref) for vector dilepton
decay, [`TwoPhotonChannel`](@ref) for two photons, [`TwoGluonChannel`](@ref),
[`ThreeGluonChannel`](@ref) for integrated gluonic rates.

## Operators

[`PhotonEmission`](@ref), [`PseudoscalarEmission`](@ref),
[`LeptonicCurrent`](@ref), [`TwoPhotonAnnihilation`](@ref),
[`GluonicAnnihilation`](@ref), `TableVReference` (internal reference backend).

The final-state type must match the operator. `LeptonicCurrent(:P_P, ...)`
requires a lepton-neutrino channel with an explicit CKM magnitude; vector
currents use the massless-pair approximation. Axial tau decay remains available
through the primitive `axial_tau_width`, in GeV.

## Amplitude inputs

`TransitionAmplitude`, `RadiativeAmplitude`,
`AnnihilationAmplitude` produced by `TwoPhotonAnnihilation`,
`StrongDecayAmplitude` from the reference backend.
For a leptonic current use the three-argument form to specify a physical final
channel: its `Vacuum()` matrix element is a current, not a decay width.

Native strong widths use `RelativisticTwoBodyNormalization`; reference
rows use `GITableVNormalization`. M1 moments are converted with photon
phase space; E1/M2 amplitudes are squared. Two-photon amplitudes are squared
and converted from GeV to MeV. Gluonic formulas expose integrated rates only.

## Mass corrections

[`mass_correction_factor`](@ref) gives a fixed-wave comparison at another mass
or momentum. Most factors multiply amplitudes or currents, so their squared
modulus does not generally include the change in the width's remaining phase
space. Gluonic correction factors multiply the width directly.

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
current = matrix_element(Vacuum(), operator, initial)
width_MeV = decay_width(MasslessLeptonPair(), operator, initial)
@assert width_MeV > 0
@assert isfinite(current.value)
```

## Related

[`matrix_element`](@ref), [`partial_waves`](@ref), `m1_radiative_width`,
`leptonic_pseudoscalar_width`, `dilepton_vector_width`,
`gluonic_annihilation_width`, `two_photon_amplitude`,
`GEV_TO_MEV`.
"""
function decay_width end

decay_width(amplitude::TransitionAmplitude) =
    partial_width(amplitude.normalization, amplitude)

function decay_width(
    final::TwoMesonChannel,
    operator::StrongDecayOperator,
    initial::TransitionState,
)
    _validate_transition(final, operator, initial)
    initial.mass_GeV <= final.first.mass_GeV + final.second.mass_GeV && return 0.0
    return decay_width(matrix_element(final, operator, initial))
end

function Base.show(io::IO, state::PhysicalState)
    print(io, "PhysicalState(\"", state.label, "\", ", state.mass_GeV, " GeV, ",
        length(state.components), " component", length(state.components) == 1 ? "" : "s", ")")
end

function Base.show(io::IO, state::ReferenceState)
    print(io, "ReferenceState(\"", state.label, "\", ", state.mass_GeV, " GeV)")
end

function Base.show(io::IO, channel::TwoMesonChannel)
    print(io, "TwoMesonChannel(", channel.first.label, ", ", channel.second.label, ")")
end

function Base.show(io::IO, ::MIME"text/plain", amplitude::TransitionAmplitude)
    println(io, "TransitionAmplitude")
    println(io, "  transition     ", amplitude.initial.label, " -> ",
        amplitude.final.first.label, " + ", amplitude.final.second.label)
    println(io, "  operator       ", nameof(typeof(amplitude.operator)))
    println(io, "  momentum       ", amplitude.momentum_GeV, " GeV")
    println(io, "  normalization  ", nameof(typeof(amplitude.normalization)))
    println(io, "  helicities     ", length(amplitude.helicity))
    println(io, "  partial waves  ", length(amplitude.partial_wave_amplitudes))
    for (wave, value) in amplitude.partial_wave_amplitudes
        println(io, "    (L,S)=(", wave.relative_L, ",", wave.channel_spin, ")  ", value)
    end
    print(io, "  terms          ", length(amplitude.terms))
    return nothing
end

Base.show(io::IO, amplitude::TransitionAmplitude) = print(
    io,
    "TransitionAmplitude(", amplitude.initial.label, " -> ",
    amplitude.final.first.label, " + ", amplitude.final.second.label, ", ",
    length(amplitude.partial_wave_amplitudes), " partial waves)",
)
