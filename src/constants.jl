# Public API (exported from GIModel.jl): (none — internal constants only)

const ALPHA_COEFFS = (0.25, 0.15, 0.20)
const ALPHA_GAMMAS = (0.5, sqrt(10.0) / 2, sqrt(1000.0) / 2)
const L_SYMBOLS = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)
const L_LABELS = Dict(v => k for (k, v) in L_SYMBOLS)
