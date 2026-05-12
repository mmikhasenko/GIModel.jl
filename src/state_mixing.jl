# Generic state-mixing bookkeeping and diagonalization.
#
# Public API (exported from GIModel.jl):
#   MixingMechanism, AntisymmetricSpinOrbit, TensorMixing, IsoscalarAnnihilation
#   BasisState, MixingBlock, MixingResult, diagonalize_mixing_block

"""
    MixingMechanism

Marker supertype for post-fixed-sector mixing mechanisms. These objects name
which paper mass-matrix contribution is being assigned at the comparison layer;
they are not radial solver paths.
"""
abstract type MixingMechanism end

"""Open-flavor unequal-mass `^1L_J`/`^3L_J` antisymmetric spin-orbit mixing."""
struct AntisymmetricSpinOrbit <: MixingMechanism end

"""Triplet same-`J`, different-`L` tensor mixing, e.g. `^3S_1`/`^3D_1`."""
struct TensorMixing <: MixingMechanism end

"""Self-conjugate isoscalar flavor/radial annihilation mixing."""
struct IsoscalarAnnihilation <: MixingMechanism end

"""
    BasisState(n, L_label, multiplicity, J; label="")

Pure spectroscopic basis state used by a mixing block. This intentionally
contains only quantum-number bookkeeping; radial matrix elements and physical
couplings are supplied by the block builder for each mixing mechanism.
"""
struct BasisState
    n::Int
    L_label::String
    multiplicity::Int
    J::Int
    label::String
    function BasisState(
        n::Integer,
        L_label::AbstractString,
        multiplicity::Integer,
        J::Integer;
        label::AbstractString = "",
    )
        nf = Int(n)
        mf = Int(multiplicity)
        jf = Int(J)
        lf = String(L_label)
        label_s = isempty(label) ? @sprintf("%d^%d%s_%d", nf, mf, lf, jf) : String(label)
        return new(nf, lf, mf, jf, label_s)
    end
end

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
    MixingResult(block, masses, vectors)

Eigen-decomposition of a [`MixingBlock`](@ref). Columns of `vectors` are
components in `block.basis`, ordered by ascending `masses`.
"""
struct MixingResult
    block::MixingBlock
    masses::Vector{Float64}
    vectors::Matrix{Float64}
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
