# =============================================================================
# Two-spin-1/2 algebra in the GI Appendix-B / Condon--Shortley convention
# =============================================================================

"""Sparse coupled quark-antiquark spin ket; projections are stored doubled."""
struct _SpinState{T<:Number}
    total_spin::Int
    twice_projection::Int
    components::Dict{Tuple{Int,Int},T}
end

function _coupled_spin(total_spin::Integer, projection::Integer)
    S = Int(total_spin)
    m = Int(projection)
    S in (0, 1) || throw(ArgumentError("two spin-1/2 constituents have S=0 or 1"))
    abs(m) <= S || throw(ArgumentError("spin projection m=$m is outside S=$S"))
    r2 = inv(sqrt(2.0))
    components = if S == 0
        m == 0 || throw(ArgumentError("a singlet has only m=0"))
        Dict((1, -1) => r2, (-1, 1) => -r2)
    elseif m == 1
        Dict((1, 1) => 1.0)
    elseif m == 0
        Dict((1, -1) => r2, (-1, 1) => r2)
    else
        Dict((-1, -1) => 1.0)
    end
    return _SpinState(S, 2m, components)
end

"""Action of the Condon--Shortley spherical Pauli component `mu=-1,0,+1`."""
function _pauli_spherical(twice_spin::Int, mu::Int)
    twice_spin in (-1, 1) || throw(ArgumentError("spin projection must be +/-1 doubled"))
    mu in -1:1 || throw(ArgumentError("spherical Pauli component must be -1, 0, or +1"))
    mu == 0 && return twice_spin, Float64(twice_spin)
    mu == 1 && return twice_spin == -1 ? (1, -sqrt(2.0)) : (twice_spin, 0.0)
    return twice_spin == 1 ? (-1, sqrt(2.0)) : (twice_spin, 0.0)
end

function _spin_factor(
    ::_QuarkEmission,
    final::_SpinState,
    initial::_SpinState,
    mu::Integer,
)
    result = 0.0
    for ((mq, mqbar), initial_coefficient) in initial.components
        final_mq, pauli = _pauli_spherical(mq, Int(mu))
        iszero(pauli) && continue
        result += conj(get(final.components, (final_mq, mqbar), 0.0)) *
                  pauli * initial_coefficient
    end
    return result
end

function _spin_factor(
    ::_AntiquarkEmission,
    final::_SpinState,
    initial::_SpinState,
    mu::Integer,
)
    result = 0.0
    for ((mq, mqbar), initial_coefficient) in initial.components
        final_mqbar, pauli = _pauli_spherical(mqbar, Int(mu))
        iszero(pauli) && continue
        result += conj(get(final.components, (mq, final_mqbar), 0.0)) *
                  pauli * initial_coefficient
    end
    return result
end
