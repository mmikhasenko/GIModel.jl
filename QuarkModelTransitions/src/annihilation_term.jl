"""
AnnihilationTerm(flavors, coefficient)

Charge/flavor coefficient for an ordered constituent sector. Mixing coefficients
belong to the PhysicalState and must not be included here a second time.

For `LeptonicCurrent(:P_P, ...)`, the coefficient `2sqrt(3)` converts the
Appendix-D pseudoscalar overlap P into `f_P/M_P = 2sqrt(3) P`. The `sqrt(3)`
is the color-singlet current factor (three colors contracted with a normalized
color state); the factor 2 belongs to this current/overlap normalization.
It is not an electric charge or a CKM element. CKM belongs to
[`LeptonNeutrinoChannel`](@ref), and state flavor/mixing weights belong to the
state. For electromagnetic vector annihilation the corresponding coefficient
is `2sqrt(3) Q_f` for S waves and additionally `sqrt(2)/3` for D waves;
[`LeptonicCurrent`](@ref)`(:electromagnetic, masses)` supplies these automatically.
Other operator kinds use their own documented normalization; `2sqrt(3)` is
not a universal coefficient for arbitrary currents.

## Example

```julia
using GIModel.QuarkModelTransitions
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
