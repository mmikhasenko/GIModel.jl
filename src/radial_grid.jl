# Public API (exported from GIModel.jl): (none)

function radial_grid(ngrid::Integer, rmax::Real)
    h = rmax / (ngrid + 1)
    collect(h:h:(ngrid*h)), h
end

function p2_operator(m::Real, L::Integer, r::AbstractVector, h::Real)
    ngrid = length(r)
    diagonal = similar(r)
    offdiag = fill(-1 / h^2, ngrid - 1)
    for i in eachindex(r)
        diagonal[i] = 2 / h^2 + L * (L + 1) / r[i]^2
    end
    SymTridiagonal(diagonal, offdiag)
end

function p2_operator(
    ::Type{FiniteDifferenceBasis},
    m::Real,
    L::Integer,
    r::AbstractVector,
    h::Real,
)
    return p2_operator(m, L, r, h)
end

function p2_operator(
    ::Type{HarmonicOscillatorBasis},
    m::Real,
    L::Integer,
    r::AbstractVector,
    h::Real,
)
    # Mesh-level spin-dependent diagnostics operate on reconstructed u(r).
    # The central HO Hamiltonian projects p² inside the oscillator subspace.
    return p2_operator(m, L, r, h)
end

function p2_operator(params::GIParameters, m::Real, L::Integer, r::AbstractVector, h::Real)
    return p2_operator(basis_type(params), m, L, r, h)
end
