# =============================================================================
# Sparse flavor states and elementary-emission topologies (GI Appendix B)
# =============================================================================

abstract type _EmissionTopology end
struct _QuarkEmission <: _EmissionTopology end
struct _AntiquarkEmission <: _EmissionTopology end

"""Normalized sparse `q qbar` flavor ket in the GI Appendix-B phase convention."""
struct _FlavorState{T<:Number}
    label::String
    components::Dict{Tuple{Symbol,Symbol},T}
end

function _FlavorState(label::AbstractString, terms)
    collected = collect(terms)
    isempty(collected) && throw(ArgumentError("a flavor state needs at least one component"))
    T = promote_type((typeof(last(term)) for term in collected)...)
    components = Dict{Tuple{Symbol,Symbol},T}()
    for term in collected
        key = first(term)
        key isa Tuple{Symbol,Symbol} || throw(ArgumentError(
            "flavor keys must be `(quark, antiquark)` symbol pairs, got $key",
        ))
        coefficient = convert(T, last(term))
        isfinite(coefficient) || throw(ArgumentError(
            "flavor coefficients must be finite, got $coefficient",
        ))
        components[key] = get(components, key, zero(T)) + coefficient
    end
    filter!(pair -> !iszero(last(pair)), components)
    norm2 = sum(abs2, values(components))
    isapprox(norm2, one(norm2); atol = 64eps(float(real(one(T)))), rtol = 0) ||
        throw(ArgumentError("flavor state `$label` is not normalized: norm²=$norm2"))
    return _FlavorState{T}(String(label), components)
end

_flavor_coefficient(state::_FlavorState, quark::Symbol, antiquark::Symbol) =
    get(state.components, (quark, antiquark), zero(eltype(values(state.components))))

"""Named normalized flavor states in the conventions of GI Eqs. (B1)--(B15)."""
function _appendix_b_flavor(name::Symbol)
    r2 = inv(sqrt(2.0))
    r3 = inv(sqrt(3.0))
    r6 = inv(sqrt(6.0))
    name === :pi_plus && return _FlavorState("pi+", [(:u, :d) => -1.0])
    name === :pi_zero && return _FlavorState("pi0", [(:u, :u) => r2, (:d, :d) => -r2])
    name === :pi_minus && return _FlavorState("pi-", [(:d, :u) => 1.0])
    name === :K_plus && return _FlavorState("K+", [(:u, :s) => -1.0])
    name === :K_zero && return _FlavorState("K0", [(:d, :s) => -1.0])
    name === :Kbar_zero && return _FlavorState("Kbar0", [(:s, :d) => -1.0])
    name === :K_minus && return _FlavorState("K-", [(:s, :u) => 1.0])
    name === :eta8 && return _FlavorState(
        "eta8", [(:u, :u) => r6, (:d, :d) => r6, (:s, :s) => -2r6],
    )
    name === :eta1 && return _FlavorState(
        "eta1", [(:u, :u) => r3, (:d, :d) => r3, (:s, :s) => r3],
    )
    name === :M_ns && return _FlavorState(
        "M_ns", [(:u, :u) => r2, (:d, :d) => r2],
    )
    name === :M_s && return _FlavorState("M_s", [(:s, :s) => 1.0])
    name === :eta && return _FlavorState(
        "eta", [(:u, :u) => 0.5, (:d, :d) => 0.5, (:s, :s) => -r2],
    )
    name === :eta_prime && return _FlavorState(
        "eta'", [(:u, :u) => 0.5, (:d, :d) => 0.5, (:s, :s) => r2],
    )
    throw(ArgumentError("unknown GI Appendix-B flavor state `$name`"))
end

function _heavy_flavor_state(heavy::Symbol, light_antiquark::Symbol)
    coefficient = light_antiquark === :d ? -1.0 : 1.0
    return _FlavorState("$heavy$(light_antiquark)bar", [
        (heavy, light_antiquark) => coefficient,
    ])
end

"""`X_q^P(final, initial)` derived from the emitted pseudoscalar flavor ket."""
function _flavor_transfer(
    ::_QuarkEmission,
    emitted::_FlavorState,
    final_flavor::Symbol,
    initial_flavor::Symbol,
)
    return sqrt(2.0) * conj(_flavor_coefficient(
        emitted, initial_flavor, final_flavor,
    ))
end

"""`X_qbar^P(final, initial)` with the anti-fundamental minus sign of Eq. (B25)."""
function _flavor_transfer(
    ::_AntiquarkEmission,
    emitted::_FlavorState,
    final_flavor::Symbol,
    initial_flavor::Symbol,
)
    return -sqrt(2.0) * conj(_flavor_coefficient(
        emitted, final_flavor, initial_flavor,
    ))
end

"""
    _flavor_factor(topology, daughter, emitted, parent)

Evaluate the flavor part of one elementary-emission topology. The spectator
flavor is enforced exactly; unsupported transitions return an exact zero.
"""
function _flavor_factor(
    topology::_QuarkEmission,
    daughter::_FlavorState,
    emitted::_FlavorState,
    parent::_FlavorState,
)
    T = promote_type(
        eltype(values(daughter.components)),
        eltype(values(emitted.components)),
        eltype(values(parent.components)),
        Float64,
    )
    result = zero(T)
    for ((initial_quark, spectator), parent_coefficient) in parent.components
        for ((final_quark, daughter_spectator), daughter_coefficient) in daughter.components
            spectator == daughter_spectator || continue
            result += conj(daughter_coefficient) *
                      _flavor_transfer(
                          topology, emitted, final_quark, initial_quark,
                      ) * parent_coefficient
        end
    end
    return result
end

function _flavor_factor(
    topology::_AntiquarkEmission,
    daughter::_FlavorState,
    emitted::_FlavorState,
    parent::_FlavorState,
)
    T = promote_type(
        eltype(values(daughter.components)),
        eltype(values(emitted.components)),
        eltype(values(parent.components)),
        Float64,
    )
    result = zero(T)
    for ((spectator, initial_antiquark), parent_coefficient) in parent.components
        for ((daughter_spectator, final_antiquark), daughter_coefficient) in
            daughter.components
            spectator == daughter_spectator || continue
            result += conj(daughter_coefficient) *
                      _flavor_transfer(
                          topology, emitted, final_antiquark, initial_antiquark,
                      ) * parent_coefficient
        end
    end
    return result
end
