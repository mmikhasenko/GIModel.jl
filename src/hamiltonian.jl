# Public API (exported from GIModel.jl): (none — use `GIModel.fn` in tests/scripts)

function nonrelativistic_hamiltonian(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
    ngrid::Integer = solver.ngrid,
    rmax::Real = solver.rmax,
)
    mu = reduced_mass(masses)
    r, h = radial_grid(ngrid, rmax)
    vdiag = potential_diagonal(params, masses.m1_GeV, masses.m2_GeV, r)
    T = promote_type(eltype(r), typeof(mu), eltype(vdiag))
    diagonal = Vector{T}(undef, length(r))
    offdiag = fill(convert(T, -1 / (2 * mu * h^2)), ngrid - 1)
    for i in eachindex(r)
        ri = r[i]
        diagonal[i] = 1 / (mu * h^2) + L * (L + 1) / (2 * mu * ri^2) + vdiag[i]
    end
    SymTridiagonal(diagonal, offdiag), r
end

function sqrt_kinetic_matrix_from_eigen(fact, m::Real)
    fact.vectors * Diagonal(sqrt.(max.(fact.values, 0) .+ m^2)) * fact.vectors'
end

lowest_eigenpairs(hamiltonian::AbstractMatrix, nlevels::Integer) =
    lowest_eigenpairs(hamiltonian, nlevels, Val(:full))

lowest_eigenpairs(
    hamiltonian::AbstractMatrix,
    nlevels::Integer,
    ::FiniteDifferenceSolver{Kinetic,Eigensolver},
) where {Kinetic,Eigensolver} =
    lowest_eigenpairs(hamiltonian, nlevels, Val(Eigensolver))

lowest_eigenpairs(
    hamiltonian::AbstractMatrix,
    nlevels::Integer,
    ::OscillatorSolver,
) = lowest_eigenpairs(hamiltonian, nlevels, Val(:full))

function lowest_eigenpairs(
    hamiltonian::AbstractMatrix{T},
    nlevels::Integer,
    ::Val{:full},
) where {T<:Real}
    fact = eigen(Symmetric(hamiltonian))
    return collect(fact.values[1:nlevels]), Matrix(fact.vectors[:, 1:nlevels])
end

function lowest_eigenpairs(
    hamiltonian::AbstractMatrix{T},
    nlevels::Integer,
    ::Val{:krylov},
) where {T<:Real}
    values, vectors, info = eigsolve(hamiltonian, nlevels, :SR; issymmetric = true)
    length(values) >= nlevels || error(
        "Krylov eigensolver converged only $(length(values)) values for nlevels=$nlevels: $info",
    )
    order = sortperm(real.(values))[1:nlevels]
    R = typeof(float(real(zero(T))))
    selected_values = collect(R, real.(values[order]))
    selected_vectors = Matrix{R}(undef, size(hamiltonian, 1), nlevels)
    for (column, index) in enumerate(order)
        selected_vectors[:, column] .= real.(vectors[index])
    end
    return selected_values, selected_vectors
end

function appendix_a_momentum_sandwich_matrix(
    params::GIParameters,
    masses::ConstituentMasses,
    r::AbstractVector,
    p2_fact,
)
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    λ = max.(p2_fact.values, 0)
    e1 = sqrt.(λ .+ m1^2)
    e2 = sqrt.(λ .+ m2^2)
    a_diag = sqrt.(1 .+ λ ./ (e1 .* e2))
    A = p2_fact.vectors * Diagonal(a_diag) * p2_fact.vectors'
    gdiag = Diagonal([smeared_coulomb_G_closed(params, m1, m2, ri) for ri in r])
    sdiag = Diagonal([smeared_confinement_S_closed(params, m1, m2, ri) for ri in r])
    Symmetric(A * gdiag * A + sdiag)
end

function relativistic_hamiltonian(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
    ngrid::Integer = solver.ngrid,
    rmax::Real = solver.rmax,
)
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    r, h = radial_grid(ngrid, rmax)
    p2 = p2_operator(params, m1, L, r, h)
    p2_fact = eigen(p2)
    kinetic =
        sqrt_kinetic_matrix_from_eigen(p2_fact, m1) +
        sqrt_kinetic_matrix_from_eigen(p2_fact, m2)
    potential = if params.central isa AppendixAMomentumSandwich
        appendix_a_momentum_sandwich_matrix(params, masses, r, p2_fact)
    else
        Diagonal(potential_diagonal(params, m1, m2, r))
    end
    Symmetric(kinetic + potential), r
end
