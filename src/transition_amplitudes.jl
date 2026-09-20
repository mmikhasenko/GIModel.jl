# =============================================================================
# Transition-amplitude domain model
# =============================================================================

abstract type TransitionOperator end
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
are independent of the [`Spectrum`](@ref) that produced it.
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

_state_mass(state::CentralState) = state.central_GeV
_state_mass(state::Union{CorrectedState,MixedState}) = state.mass_GeV

function _physical_state_provenance(spec::Spectrum, state)
    channels = Tuple((meson.flavor1, meson.flavor2) for meson in spec.channels)
    return (
        source = :spectrum,
        spectrum_stage = nameof(typeof(state)),
        channels = channels,
        state_label = state.label,
    )
end

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

_state_sort_key(state::TransitionState) = repr(_state_identity(state))
_same_external_state(a::TransitionState, b::TransitionState) =
    _state_identity(a) == _state_identity(b)

"""A canonically ordered pair of daughter states, without a user-selected partial wave."""
struct TwoMesonChannel{A<:TransitionState,B<:TransitionState}
    first::A
    second::B
    function TwoMesonChannel(first::A, second::B) where {
        A<:TransitionState,
        B<:TransitionState,
    }
        if _state_sort_key(second) < _state_sort_key(first)
            return new{B,A}(second, first)
        end
        return new{A,B}(first, second)
    end
end

"""Relative orbital angular momentum `L` and coupled daughter spin `S`."""
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

Return the helicity-to-partial-wave map. Phase 1 implements the unique
spin-zero-daughter projection. General Jacob--Wick projection is intentionally
gated on the Phase 3 convention ledger.
"""
function partial_wave_projection(final::TwoMesonChannel, initial::TransitionState)
    j1, _ = _state_quantum_numbers(final.first)
    j2, _ = _state_quantum_numbers(final.second)
    waves = Tuple(allowed_partial_waves(final, initial))
    if j1 == 0 && j2 == 0
        helicities = (HelicityLabel(0, 0),)
        return PartialWaveProjection(helicities, waves, ones(Float64, length(waves), 1))
    end
    throw(ArgumentError(
        "general helicity-to-partial-wave projection is not implemented yet; " *
        "the Phase 3 Jacob-Wick convention gate must land first",
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

abstract type TransitionKinematics end
struct OnShell <: TransitionKinematics end

struct CMKinematics{T<:Number} <: TransitionKinematics
    momentum_GeV::T
    function CMKinematics(momentum_GeV::Number)
        momentum = float(momentum_GeV)
        isfinite(real(momentum)) && isfinite(imag(momentum)) || throw(ArgumentError(
            "center-of-mass momentum must be finite, got $momentum_GeV",
        ))
        momentum isa Real && momentum < 0 && throw(ArgumentError(
            "real center-of-mass momentum must be non-negative, got $momentum_GeV",
        ))
        return new{typeof(momentum)}(momentum)
    end
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
        " GeV; use explicit CMKinematics(k) for an off-shell vertex",
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

_resolve_kinematics(final, initial, kinematics::CMKinematics) = kinematics
_resolve_kinematics(final, initial, ::OnShell) =
    CMKinematics(_on_shell_momentum(final, initial))

abstract type AmplitudeNormalization end
struct RelativisticTwoBodyNormalization <: AmplitudeNormalization end
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
provenance. Collections are tuples so no field erases its element type.
"""
struct TransitionAmplitude{O,I,F,K,N,H<:Tuple,W<:Tuple,R<:Tuple,P}
    operator::O
    initial::I
    final::F
    kinematics::K
    normalization::N
    helicity::H
    partial_wave_amplitudes::W
    terms::R
    provenance::P
end

partial_waves(amplitude::TransitionAmplitude) =
    [first(item) for item in amplitude.partial_wave_amplitudes]

function Base.getindex(amplitude::TransitionAmplitude, wave::PartialWave)
    index = findfirst(item -> first(item) == wave, amplitude.partial_wave_amplitudes)
    isnothing(index) && throw(KeyError(wave))
    return last(amplitude.partial_wave_amplitudes[index])
end

function matrix_element end
function partial_width end

_validate_transition(final, operator, initial) = nothing

function matrix_element(
    final::TwoMesonChannel,
    operator::StrongDecayOperator,
    initial::TransitionState;
    kinematics::TransitionKinematics = OnShell(),
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

decay_width(amplitude::TransitionAmplitude) =
    partial_width(amplitude.normalization, amplitude)

function decay_width(
    final::TwoMesonChannel,
    operator::StrongDecayOperator,
    initial::TransitionState,
)
    _validate_transition(final, operator, initial)
    initial.mass_GeV <= final.first.mass_GeV + final.second.mass_GeV && return 0.0
    return decay_width(matrix_element(final, operator, initial; kinematics = OnShell()))
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
    momentum = hasproperty(amplitude.kinematics, :momentum_GeV) ?
        string(amplitude.kinematics.momentum_GeV, " GeV") :
        string(nameof(typeof(amplitude.kinematics)))
    println(io, "  momentum       ", momentum)
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
