"""
AnnihilationTerm(flavors, coefficient)

Charge/flavor coefficient for an ordered constituent sector. Mixing coefficients
belong to the PhysicalState and must not be included here a second time.

## Example

```julia
using QuarkModelTransitions
term = AnnihilationTerm((:c, :c), 4/9)
@assert term.coefficient == 4/9
```

## Related

[`LeptonicCurrent`](@ref), [`TwoPhotonAnnihilation`](@ref).
"""
struct AnnihilationTerm
    flavors::Tuple{Symbol,Symbol}
    coefficient::Float64
    function AnnihilationTerm(flavors::Tuple{Symbol,Symbol}, coefficient::Real)
        isfinite(coefficient) || throw(ArgumentError("coefficient must be finite"))
        new(flavors, Float64(coefficient))
    end
end

function _annihilation_inputs(masses, terms)
    all(isfinite(m) && m > 0 for m in values(masses)) ||
        throw(ArgumentError("constituent masses must be finite and positive"))
    ts = terms isa AnnihilationTerm ? (terms,) : Tuple(terms)
    !isempty(ts) && all(t isa AnnihilationTerm for t in ts) ||
        throw(ArgumentError("supply AnnihilationTerm coefficients"))
    length(unique(t.flavors for t in ts)) == length(ts) ||
        throw(ArgumentError("duplicate flavor terms"))
    return copy(masses), ts
end
