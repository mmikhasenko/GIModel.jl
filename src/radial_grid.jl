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

# `p2_operator` used to have one method per basis, both with the same body: the
# mesh `p²` is the mesh `p²` whichever way you intend to solve. The convenience
# form taking `params` stays, since most callers hold parameters rather than a
# mass; it too is method-free.
function p2_operator(::GIParameters, m::Real, L::Integer, r::AbstractVector, h::Real)
    return p2_operator(m, L, r, h)
end
