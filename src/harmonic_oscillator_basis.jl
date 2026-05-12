# Paper-style harmonic-oscillator radial expansion for basis-comparison audits.
#
# This path intentionally reuses the same uniform mesh for operator quadrature
# and for reconstructed wavefunctions. The basis is finite and diagnostic, but
# the nonlocal `f(p) V(r) f(p)` operators are assembled in the oscillator
# subspace rather than on the full FD coordinate basis.

const HO_DEFAULT_NBASIS = 24
const HO_BETA_GRID = collect(0.25:0.10:2.35)

function generalized_laguerre(n::Integer, α::Real, x::Real)
    n == 0 && return 1.0
    n == 1 && return 1.0 + α - x
    lm2 = 1.0
    lm1 = 1.0 + α - x
    for k = 2:n
        lk = ((2k - 1 + α - x) * lm1 - (k - 1 + α) * lm2) / k
        lm2 = lm1
        lm1 = lk
    end
    return lm1
end

function ho_reduced_radial(nr::Integer, L::Integer, β::Real, r::Real)
    x = (β * r)^2
    α = L + 0.5
    norm = sqrt(2 * β * gamma(nr + 1) / gamma(nr + L + 1.5))
    return norm * (β * r)^(L + 1) * exp(-0.5 * x) * generalized_laguerre(nr, α, x)
end

function ho_basis_matrix(L::Integer, β::Real, r::AbstractVector, nbasis::Integer)
    U = Matrix{Float64}(undef, length(r), nbasis)
    for j = 1:nbasis
        nr = j - 1
        for i in eachindex(r)
            U[i, j] = ho_reduced_radial(nr, L, β, float(r[i]))
        end
    end
    return U
end

function orthonormalize_physical_basis(U::AbstractMatrix, h::Real)
    F = qr(sqrt(h) .* Matrix(U))
    cols = min(size(U, 2), size(F.Q, 2))
    return Matrix(F.Q[:, 1:cols]) ./ sqrt(h)
end

function projected_matrix(U::AbstractMatrix, h::Real, A::AbstractMatrix)
    return Symmetric(h * (U' * (A * U)))
end

function projected_diagonal(U::AbstractMatrix, h::Real, values::AbstractVector)
    return Symmetric(h * (U' * (values .* U)))
end

function oscillator_kinetic_matrix(p2_basis::AbstractMatrix, m::Real)
    fact = eigen(Symmetric(p2_basis))
    λ = max.(fact.values, 0.0)
    return Symmetric(fact.vectors * Diagonal(sqrt.(λ .+ m^2)) * fact.vectors')
end

function oscillator_momentum_factor_matrix(
    p2_basis::AbstractMatrix,
    m1::Real,
    m2::Real;
    power::Real = 0.5,
)
    fact = eigen(Symmetric(p2_basis))
    λ = max.(fact.values, 0.0)
    e1 = sqrt.(λ .+ m1^2)
    e2 = sqrt.(λ .+ m2^2)
    diag = (1 .+ λ ./ (e1 .* e2)) .^ power
    return Symmetric(fact.vectors * Diagonal(diag) * fact.vectors')
end

function oscillator_hamiltonian_for_beta(
    params::GIParameters{HarmonicOscillatorBasis},
    masses::ConstituentMasses,
    L::Integer,
    r::AbstractVector,
    h::Real,
    β::Real;
    nbasis::Integer = HO_DEFAULT_NBASIS,
)
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    U = orthonormalize_physical_basis(ho_basis_matrix(L, β, r, nbasis), h)
    p2_grid = p2_operator(m1, L, r, h)
    p2_basis = projected_matrix(U, h, p2_grid)
    kinetic = oscillator_kinetic_matrix(p2_basis, m1) + oscillator_kinetic_matrix(p2_basis, m2)
    potential = if central_potential_mode(params) == :appendix_a_momentum_sandwich
        A = oscillator_momentum_factor_matrix(p2_basis, m1, m2; power = 0.5)
        g = projected_diagonal(
            U,
            h,
            [smeared_coulomb_G_closed(params, m1, m2, ri) for ri in r],
        )
        s = projected_diagonal(
            U,
            h,
            [smeared_confinement_S_closed(params, m1, m2, ri) for ri in r],
        )
        Symmetric(A * g * A + s)
    else
        projected_diagonal(U, h, potential_diagonal(params, m1, m2, r))
    end
    return Symmetric(kinetic + potential), U
end

function oscillator_beta_candidates(masses::ConstituentMasses, L::Integer)
    # The broad grid is deliberately not tuned per sector. It spans diffuse
    # light states and compact bottomonia well enough for an audit comparison.
    return HO_BETA_GRID
end

function oscillator_channel_solution(
    params::GIParameters{HarmonicOscillatorBasis},
    masses::ConstituentMasses,
    L::Integer;
    nlevels::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    nbasis::Integer = max(HO_DEFAULT_NBASIS, nlevels + 4),
)
    r, h = radial_grid(ngrid, rmax)
    best = nothing
    for β in oscillator_beta_candidates(masses, L)
        H, U = oscillator_hamiltonian_for_beta(params, masses, L, r, h, β; nbasis = nbasis)
        vals, vecs = lowest_eigenpairs(Matrix(H), nlevels; eigensolver = :full)
        # For an orthogonal set in a fixed sector, use one beta for all reported
        # levels. Following the paper's practical convention, choose the beta
        # that minimizes the last requested state rather than overfitting the
        # ground state.
        if isnothing(best) || vals[end] < best.values[end]
            best = (beta = β, values = vals, coeffs = vecs, basis = U)
        end
    end
    waves = best.basis * best.coeffs
    for col in axes(waves, 2)
        nrm = sqrt(sum(abs2, waves[:, col]) * h)
        nrm > 0 && (waves[:, col] ./= nrm)
    end
    return collect(best.values), Matrix(waves), r
end
