# Public API (exported from GIModel.jl):
#   channel_solution

function channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer;
    nlevels::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    hamiltonian, r = if kinetic == :relativistic
        relativistic_hamiltonian(params, masses, L; ngrid = ngrid, rmax = rmax)
    elseif kinetic == :nonrelativistic
        nonrelativistic_hamiltonian(params, masses, L; ngrid = ngrid, rmax = rmax)
    else
        error("unknown kinetic mode: $kinetic")
    end
    values, vectors = lowest_eigenpairs(hamiltonian, nlevels; eigensolver = eigensolver)
    if kinetic == :relativistic
        values, vectors, r
    else
        values .+ (masses.m1_GeV + masses.m2_GeV), vectors, r
    end
end

function solve_channel(args...; kwargs...)
    values, _vectors, _r = channel_solution(args...; kwargs...)
    values
end
