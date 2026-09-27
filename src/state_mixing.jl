# Generic state-mixing bookkeeping and diagonalization.
#
# Public API (exported from GIModel.jl):
#   MixingMechanism, AntisymmetricSpinOrbit, TensorMixing, IsoscalarAnnihilation
#   BasisState, MixingBlock, MixingResult, diagonalize_mixing_block

"""
    MixingMechanism

Marker supertype for post-fixed-sector mixing mechanisms. These objects name
which mass-matrix contribution is being applied to model states; they are not
radial solver paths. The same bookkeeping covers intra-channel spin mixing and
cross-channel flavor annihilation.
"""
abstract type MixingMechanism end

"""Open-flavor unequal-mass `^1L_J`/`^3L_J` antisymmetric spin-orbit mixing."""
struct AntisymmetricSpinOrbit <: MixingMechanism end

"""Triplet same-`J`, different-`L` tensor mixing, e.g. `^3S_1`/`^3D_1`."""
struct TensorMixing <: MixingMechanism end

"""Self-conjugate isoscalar flavor/radial annihilation mixing."""
struct IsoscalarAnnihilation <: MixingMechanism end

"""
    BasisState(n, L_label, multiplicity, J; label="", flavors=nothing)

Quantum labels for a requested level or mixing component: radial `n ≥ 1`,
orbital `L_label` ("S", "P", …), spin `multiplicity = 2S + 1`, and total `J`.
Multiplicity is `1` for singlets or `3` for triplets; it selects spin interactions
and mixing partners.

[`spectrum_levels`](@ref) leaves `flavors = nothing` for reuse across mesons.
Solvers attach `(flavor1, flavor2)` to a copy. This distinguishes nonstrange and
strange components in [`compute_isoscalar_spectrum`](@ref).

## Example

```julia
using GIModel
singlet = BasisState(1, "P", 1, 1)  # 1^1P_1, S = 0
triplet = BasisState(1, "P", 3, 1)  # 1^3P_1, S = 1
strange = BasisState(1, "P", 3, 1; flavors=(:s, :s))
```

## Related

[`compute_spectrum`](@ref), [`spectrum_state`](@ref), [`MixingBlock`](@ref),
[`physical_components`](@ref).
"""
struct BasisState
    n::Int
    L_label::String
    multiplicity::Int
    J::Int
    label::String
    flavors::Union{Nothing,Tuple{Symbol,Symbol}}
    function BasisState(
        n::Integer,
        L_label::AbstractString,
        multiplicity::Integer,
        J::Integer;
        label::AbstractString = "",
        flavors = nothing,
    )
        nf = Int(n)
        mf = Int(multiplicity)
        jf = Int(J)
        lf = String(L_label)
        nf >= 1 || throw(ArgumentError("BasisState: n must be positive"))
        mf in (1, 3) || throw(ArgumentError(
            "BasisState: multiplicity must be 1 or 3",
        ))
        jf >= 0 || throw(ArgumentError("BasisState: J must be non-negative"))
        haskey(L_SYMBOLS, lf) || throw(ArgumentError(
            "BasisState: unknown orbital label `$lf`; expected one of S, P, D, F, G",
        ))
        # J must be reachable by coupling L and S = (multiplicity - 1)/2; otherwise
        # the fixed-sector Hamiltonian would silently use meaningless angular factors.
        lnum, snum = L_SYMBOLS[lf], (mf - 1) ÷ 2
        abs(lnum - snum) <= jf <= lnum + snum || throw(ArgumentError(
            "BasisState: J = $jf cannot be formed from L = $lnum and S = $snum " *
            "(allowed: $(abs(lnum - snum)):$(lnum + snum))",
        ))
        label_s = isempty(label) ? @sprintf("%d^%d%s_%d", nf, mf, lf, jf) : String(label)
        flavor_pair = if isnothing(flavors)
            nothing
        else
            length(flavors) == 2 || throw(ArgumentError(
                "BasisState: `flavors` must contain exactly two symbols",
            ))
            (_canonical_flavor(Symbol(flavors[1])), _canonical_flavor(Symbol(flavors[2])))
        end
        return new(nf, lf, mf, jf, label_s, flavor_pair)
    end
end

_with_flavors(state::BasisState, meson::Meson; label::AbstractString = state.label) =
    BasisState(
        state.n, state.L_label, state.multiplicity, state.J;
        label = label,
        flavors = (meson.flavor1, meson.flavor2),
    )

"""
    MixingBlock(name, basis, matrix; mechanism="", source="", notes="")

Hermitian mass matrix over pure [`BasisState`](@ref) entries. The matrix is in
GeV and may be any size; small `2x2` same-`J` spin-orbit blocks and larger
tensor/flavor/radial mixing blocks use the same representation.
"""
struct MixingBlock
    name::String
    basis::Vector{BasisState}
    matrix::Matrix{Float64}
    mechanism::String
    source::String
    notes::String
    function MixingBlock(
        name::AbstractString,
        basis::AbstractVector{BasisState},
        matrix::AbstractMatrix{<:Real};
        mechanism::AbstractString = "",
        source::AbstractString = "",
        notes::AbstractString = "",
    )
        n = length(basis)
        n > 0 || throw(ArgumentError("MixingBlock $name: basis must not be empty"))
        identities = [
            (state.n, state.L_label, state.multiplicity, state.J, state.flavors) for
            state in basis
        ]
        length(unique(identities)) == n || throw(ArgumentError(
            "MixingBlock `$name`: basis contains duplicate spectroscopic/flavor identities",
        ))
        size(matrix) == (n, n) ||
            throw(ArgumentError(
                "MixingBlock `$name`: matrix size $(size(matrix)) does not match basis length $n",
            ))
        mat = Matrix{Float64}(matrix)
        isapprox(mat, mat'; rtol = 1e-12, atol = 1e-12) ||
            throw(ArgumentError("MixingBlock `$name`: matrix must be symmetric/Hermitian"))
        return new(
            String(name),
            collect(basis),
            mat,
            String(mechanism),
            String(source),
            String(notes),
        )
    end
end

"""
    MixingResult(block, masses, vectors[, pole_matrices])

Solution of a [`MixingBlock`](@ref). Columns of `vectors` are components in
`block.basis`, ordered by ascending `masses`. `pole_matrices` is empty for an
ordinary Hermitian eigensystem. It contains one effective matrix per pole only
for the mass-dependent Eq. (18b) fixed-point problem.
"""
struct MixingResult
    block::MixingBlock
    masses::Vector{Float64}
    vectors::Matrix{Float64}
    pole_matrices::Vector{Matrix{Float64}}
    function MixingResult(
        block::MixingBlock,
        masses::AbstractVector{<:Real},
        vectors::AbstractMatrix{<:Real},
        pole_matrices::AbstractVector{<:AbstractMatrix} = Matrix{Float64}[],
    )
        n = length(block.basis)
        length(masses) == n || throw(ArgumentError(
            "MixingResult: mass count does not match block basis",
        ))
        values = collect(Float64, masses)
        all(isfinite, values) || throw(ArgumentError(
            "MixingResult: masses must be finite",
        ))
        issorted(values) || throw(ArgumentError(
            "MixingResult: masses must be ordered ascending with the vector columns",
        ))
        size(vectors) == (n, n) || throw(ArgumentError(
            "MixingResult: vector matrix size does not match block basis",
        ))
        all(isfinite, vectors) || throw(ArgumentError(
            "MixingResult: vectors must be finite",
        ))
        all(
            isapprox(sum(abs2, view(vectors, :, column)), 1.0; rtol = 1e-8, atol = 1e-10)
            for column in axes(vectors, 2)
        ) || throw(ArgumentError(
            "MixingResult: every vector column must have unit norm",
        ))
        isempty(pole_matrices) || length(pole_matrices) == n || throw(ArgumentError(
            "MixingResult: pole_matrices must be empty or contain one matrix per pole",
        ))
        matrices = Matrix{Float64}[Matrix{Float64}(matrix) for matrix in pole_matrices]
        all(size(matrix) == (n, n) for matrix in matrices) || throw(ArgumentError(
            "MixingResult: every pole matrix must match the block basis",
        ))
        return new(
            block,
            values,
            Matrix{Float64}(vectors),
            matrices,
        )
    end
end

"""
    diagonalize_mixing_block(block; phase_anchor=1)

Diagonalize a mixing block. Eigenvector phases are fixed by making the
`phase_anchor` component non-negative when possible, so reports remain stable.
"""
function diagonalize_mixing_block(block::MixingBlock; phase_anchor::Integer = 1)
    n = length(block.basis)
    1 <= phase_anchor <= n ||
        throw(ArgumentError("phase_anchor=$phase_anchor outside 1:$n"))
    fact = eigen(Symmetric(block.matrix))
    vectors = Matrix(fact.vectors)
    for col in axes(vectors, 2)
        if vectors[phase_anchor, col] < 0
            vectors[:, col] .*= -1
        end
    end
    return MixingResult(block, collect(Float64, fact.values), vectors)
end
