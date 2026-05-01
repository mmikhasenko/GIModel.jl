function channel_solution(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::Integer;
    nlevels::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    hamiltonian, r =
        if kinetic == :relativistic
            relativistic_hamiltonian(params, m1, m2, L; ngrid = ngrid, rmax = rmax)
        elseif kinetic == :nonrelativistic
            nonrelativistic_hamiltonian(params, m1, m2, L; ngrid = ngrid, rmax = rmax)
        else
            error("unknown kinetic mode: $kinetic")
        end
    values, vectors = lowest_eigenpairs(hamiltonian, nlevels; eigensolver = eigensolver)
    if kinetic == :relativistic
        values, vectors, r
    else
        values .+ (m1 + m2), vectors, r
    end
end

function solve_channel(args...; kwargs...)
    values, _vectors, _r = channel_solution(args...; kwargs...)
    values
end
