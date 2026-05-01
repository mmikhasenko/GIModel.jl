# Public API (exported from GIModel.jl): (none — use `GIModel.fn` in tests/scripts)

function nonrelativistic_hamiltonian(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::Integer;
    ngrid::Integer = 900,
    rmax::Real = 24.0,
)
    mu = reduced_mass(m1, m2)
    r, h = radial_grid(ngrid, rmax)
    diagonal = similar(r)
    offdiag = fill(-1 / (2 * mu * h^2), ngrid - 1)
    vdiag = potential_diagonal(params, m1, m2, r)
    for i in eachindex(r)
        ri = r[i]
        diagonal[i] = 1 / (mu * h^2) + L * (L + 1) / (2 * mu * ri^2) + vdiag[i]
    end
    SymTridiagonal(diagonal, offdiag), r
end

function sqrt_kinetic_matrix_from_eigen(fact, m::Real)
    fact.vectors * Diagonal(sqrt.(max.(fact.values, 0) .+ m^2)) * fact.vectors'
end

function lowest_eigenpairs(
    hamiltonian::AbstractMatrix,
    nlevels::Integer;
    eigensolver::Symbol = :full,
)
    if eigensolver == :full
        fact = eigen(hamiltonian)
        return fact.values[1:nlevels], fact.vectors[:, 1:nlevels]
    elseif eigensolver == :krylov
        values, vectors, info = eigsolve(hamiltonian, nlevels, :SR; issymmetric = true)
        length(values) >= nlevels || error(
            "Krylov eigensolver converged only $(length(values)) values for nlevels=$nlevels: $info",
        )
        order = sortperm(real.(values))[1:nlevels]
        return real.(values[order]), hcat(vectors[order]...)
    else
        error("unknown eigensolver: $eigensolver")
    end
end

function appendix_a_momentum_sandwich_matrix(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector,
    p2_fact,
)
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
    m1::Real,
    m2::Real,
    L::Integer;
    ngrid::Integer = 450,
    rmax::Real = 24.0,
)
    r, h = radial_grid(ngrid, rmax)
    p2 = p2_operator(m1, L, r, h)
    p2_fact = eigen(p2)
    kinetic =
        sqrt_kinetic_matrix_from_eigen(p2_fact, m1) +
        sqrt_kinetic_matrix_from_eigen(p2_fact, m2)
    potential = if central_potential_mode(params) == :appendix_a_momentum_sandwich
        appendix_a_momentum_sandwich_matrix(params, m1, m2, r, p2_fact)
    else
        Diagonal(potential_diagonal(params, m1, m2, r))
    end
    Symmetric(kinetic + potential), r
end
