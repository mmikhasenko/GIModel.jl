# Public API (exported from GIModel.jl):
#   orbital_angular_momentum, orbital_label

const ALPHA_COEFFS = (0.25, 0.15, 0.20)
const ALPHA_GAMMAS = (0.5, sqrt(10.0) / 2, sqrt(1000.0) / 2)
const L_SYMBOLS = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)
const L_LABELS = Dict(v => k for (k, v) in L_SYMBOLS)

"""Return the orbital angular momentum encoded by a spectroscopic letter."""
function orbital_angular_momentum(label::AbstractString)
    key = String(label)
    haskey(L_SYMBOLS, key) || throw(ArgumentError("unsupported orbital label `$label`"))
    return L_SYMBOLS[key]
end

"""Return the spectroscopic orbital letter for angular momentum `L`."""
function orbital_label(L::Integer)
    value = Int(L)
    haskey(L_LABELS, value) || throw(ArgumentError("unsupported orbital angular momentum `$L`"))
    return L_LABELS[value]
end
