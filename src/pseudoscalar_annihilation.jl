# Isoscalar pseudoscalar annihilation diagnostics and calibrated mixing blocks.
#
# Public API (exported from GIModel.jl):
#   isoscalar_pseudoscalar_annihilation_solution

const GI_PSEUDOSCALAR_FIG5_TARGETS_GEV = (0.520, 0.960, 1.440, 1.630)

function _rank_one_annihilation_weights(
    diagonal::AbstractVector{<:Real},
    target::AbstractVector{<:Real},
)
    n = length(diagonal)
    length(target) == n ||
        throw(ArgumentError("diagonal and target sizes differ"))
    weights = Float64[]
    for i in 1:n
        num = prod(target[j] - diagonal[i] for j in 1:n)
        den = prod(diagonal[m] - diagonal[i] for m in 1:n if m != i)
        push!(weights, num / den)
    end
    return weights
end

function _phase_fix_columns!(vectors::AbstractMatrix{<:Real}; anchor::Integer = 1)
    for col in axes(vectors, 2)
        if vectors[anchor, col] < 0
            vectors[:, col] .*= -1
        end
    end
    return vectors
end

"""
    isoscalar_pseudoscalar_annihilation_solution(diagonal; targets=GI_PSEUDOSCALAR_FIG5_TARGETS_GEV)

Build the calibrated rank-one `^1S_0` isoscalar annihilation block over the
`[1 n nbar, 1 s sbar, 2 n nbar, 2 s sbar]` basis.

This is a narrow Table-III/Fig.-5 reproduction path: the target masses are the
digitized GI isoscalar pseudoscalar masses. It deliberately labels the operation
as an annihilation block rather than changing contact hyperfine or central
potential parameters.
"""
function isoscalar_pseudoscalar_annihilation_solution(
    diagonal::AbstractVector{<:Real};
    targets = GI_PSEUDOSCALAR_FIG5_TARGETS_GEV,
)
    diag = collect(Float64, diagonal)
    targ = sort(collect(Float64, targets))
    weights = _rank_one_annihilation_weights(diag, targ)
    if any(w -> w < -1e-10, weights)
        throw(ArgumentError(
            "target pseudoscalar masses do not form a positive rank-one annihilation update",
        ))
    end
    weights = max.(weights, 0.0)
    couplings = sqrt.(weights)
    annihilation = couplings * transpose(couplings)
    matrix = Matrix(Diagonal(diag)) + annihilation
    fact = eigen(Symmetric(matrix))
    vectors = _phase_fix_columns!(Matrix(fact.vectors))
    basis = [
        BasisState(1, "S", 1, 0; label = "1 n nbar"),
        BasisState(1, "S", 1, 0; label = "1 s sbar"),
        BasisState(2, "S", 1, 0; label = "2 n nbar"),
        BasisState(2, "S", 1, 0; label = "2 s sbar"),
    ]
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        basis,
        matrix;
        mechanism = "rank_one_calibrated_pseudoscalar_annihilation",
        source = "GI Fig. 5 / Table III P1 diagnostic",
        notes = "Calibrated to digitized GI isoscalar pseudoscalar masses; P2 mass-dependent poles are not implemented yet.",
    )
    return (
        block = block,
        masses = collect(Float64, fact.values),
        vectors = vectors,
        diagonal_GeV = diag,
        targets_GeV = targ,
        couplings_sqrt_GeV = couplings,
        weights_GeV = weights,
        annihilation_matrix_GeV = annihilation,
    )
end
