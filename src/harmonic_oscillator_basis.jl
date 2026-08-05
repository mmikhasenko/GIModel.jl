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

"""
    ho_p2_matrix(L, β, nbasis) -> SymTridiagonal

Exact matrix elements of `p²` in the 3D harmonic-oscillator basis — the first
ingredient of Eq. (A17), with no mesh anywhere.

The oscillator Hamiltonian `H = p²/2μ + ½μω²r²` is diagonal with eigenvalue
`(2n+L+3/2)ω`, and `r²` is tridiagonal with known elements, so

    p² = 2μH − β⁴r²,    β² = μω

is tridiagonal too:

    ⟨n|p²|n⟩   = β² (2n + L + 3/2)
    ⟨n|p²|n+1⟩ = β² √((n+1)(n + L + 3/2))

Verified against the mesh projection this is intended to replace: the difference
falls from 4.3e-3 at `(ngrid, rmax) = (450, 24)` to 7.4e-5 at `(8000, 56)`,
i.e. it is the mesh's error and not this formula's.

**Not yet wired into `oscillator_hamiltonian_for_beta`.** Substituting it there
moves the charmonium `1S` by 9.4 MeV and turns a 0.01 MeV finite-difference /
oscillator agreement into a 9.4 MeV disagreement, with the oscillator result
landing *below* the finite-difference one — the wrong side for a variational
calculation in a finite basis. Two candidate explanations were tested and
rejected: the mesh basis is orthonormal to 1.2e-15 (so quadrature error in the
basis is not it), and the QR sign convention is a uniform -1 that cancels in
the matrix. The remaining suspect is that Eq. (A17) cannot be done half
analytically: the potential side is still projected through the mesh, and
mixing an exact momentum side with an approximate position side need not be
variational. Resolving that is the position-space half of A3', not this
function.
"""
function ho_p2_matrix(L::Integer, β::Real, nbasis::Integer)
    nbasis >= 1 || throw(ArgumentError("ho_p2_matrix: nbasis must be ≥ 1"))
    β > 0 || throw(ArgumentError("ho_p2_matrix: β must be positive"))
    b2 = float(β)^2
    diagonal = [b2 * (2n + L + 1.5) for n = 0:(nbasis-1)]
    offdiag = [b2 * sqrt((n + 1) * (n + L + 1.5)) for n = 0:(nbasis-2)]
    return SymTridiagonal(diagonal, offdiag)
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
    potential = if params.central isa AppendixAMomentumSandwich
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

# The beta grid is fixed and does not adapt to the masses. If the variational
# optimum lands on the LAST candidate, the true optimum lies outside the grid:
# the basis is too diffuse to represent a state this compact, and the result is
# silently under-resolved rather than obviously wrong. Bottomonium picks
# beta = 1.55 and m_Q = 8 GeV picks 1.95, so this stays quiet for every sector
# the paper uses and first fires around m_Q ~ 15 GeV.
function _warn_if_beta_railed(best_beta::Real, candidates, masses::ConstituentMasses, L::Integer)
    best_beta == last(candidates) || return nothing
    @warn """
    Harmonic-oscillator basis railed: the optimal beta hit the top of the fixed \
    HO_BETA_GRID ($(last(candidates))), so this state is more compact than the \
    basis can represent and the result is under-resolved. Extend HO_BETA_GRID \
    for constituent masses well above bottomonium.""" m1 = masses.m1_GeV m2 =
        masses.m2_GeV L maxlog = 1
    return nothing
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
    solver::RadialSolver = RadialSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
    ngrid::Integer = solver.ngrid,
    rmax::Real = solver.rmax,
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
    _warn_if_beta_railed(best.beta, oscillator_beta_candidates(masses, L), masses, L)
    waves = physically_normalized_waves(best.basis * best.coeffs, h)
    return collect(best.values), Matrix(waves), r
end

"""
    ho_full_distorted_states(params, masses, L, V; nlevels, ngrid, rmax, nbasis)
        -> (values, waves, r)

Paper-order spin-distorted radial waves by **full diagonalization** of the
central-plus-spin Hamiltonian in the finite harmonic-oscillator subspace — the
paper's literal method. For each β candidate the spin-dependent grid operator
`V` is projected into the oscillator basis and added to the central `H`, and the
full `H + V` is diagonalized; the paper β convention (one β per sector,
minimizing the `nlevels`-th state) is applied to `H + V`. Returned waves are
physically normalized (`∫u² dr = 1`).

This differs from [`ho_first_order_distorted_states`](@ref) in resumming `V`
within the finite basis rather than truncating at first order. The two coincide
for heavy-quark spin splittings (where `V` is a small perturbation) but diverge
where `V` is large: for the light `¹S₀` nonstrange sector the huge attractive
contact term collapses the mass (`≈0.10 GeV`, matching the fine-grid FD
resummation), whereas first-order PT badly overestimates it (`≈0.28 GeV`). The
finite basis — not a perturbation order — is the mechanism: it resums less than
the fine FD grid (so the Table VII heavy gluonic ratios land below the
nonperturbative-FD overshoot) yet fully for the light pion. This is the single
treatment that harmonizes the whole Table VII audit (gluonic distortion +
light-pseudoscalar leptonic rows) on the paper's own basis.
"""
function resummed_channel_solution(
    params::GIParameters{HarmonicOscillatorBasis},
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix;
    solver::RadialSolver = RadialSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
    ngrid::Integer = solver.ngrid,
    rmax::Real = solver.rmax,
    nbasis::Integer = max(HO_DEFAULT_NBASIS, nlevels + 4),
)
    r, h = radial_grid(ngrid, rmax)
    size(V, 1) == length(r) || error("V must live on the (ngrid, rmax) mesh")
    best = nothing
    for β in oscillator_beta_candidates(masses, L)
        H, U = oscillator_hamiltonian_for_beta(params, masses, L, r, h, β; nbasis = nbasis)
        Vproj = projected_matrix(U, h, Matrix(V))
        F = eigen(Symmetric(Matrix(H) + Matrix(Vproj)))
        if isnothing(best) || F.values[nlevels] < best.values[nlevels]
            best = (beta = β, values = F.values, coeffs = F.vectors, basis = U)
        end
    end
    _warn_if_beta_railed(best.beta, oscillator_beta_candidates(masses, L), masses, L)
    waves = physically_normalized_waves(best.basis * best.coeffs, h)
    return collect(best.values[1:nlevels]), Matrix(waves[:, 1:nlevels]), r
end

"""
    ho_first_order_distorted_states(params, masses, L, V; nlevels, ngrid, rmax, nbasis)
        -> (values, waves, r)

Paper-order spin-distorted radial waves: diagonalize the central Hamiltonian in
the harmonic-oscillator subspace (paper β convention — one β per sector,
minimizing the `nlevels`-th state), then treat the spin-dependent grid operator
`V` in first-order perturbation theory within that eigenbasis:

    |n⟩₁ = |n⟩ + Σ_{k≠n} |k⟩ ⟨k|V|n⟩ / (Eₙ - Eₖ),   Mₙ = Eₙ + ⟨n|V|n⟩.

`V` is a dense operator on the same uniform mesh (e.g.
[`contact_hyperfine_operator`](@ref GIModel.contact_hyperfine_operator) for
S-waves or [`fine_structure_grid_operator`](@ref) for `³P_J`). Returned waves
are physically normalized (`∫u² dr = 1`). This is the W6-validated paper-order
treatment: resumming `V` nonperturbatively (FD or HO) overshoots the paper's
Table VII wavefunction-at-origin distortions, while this first-order form
reproduces them.
"""
function ho_first_order_distorted_states(
    params::GIParameters{HarmonicOscillatorBasis},
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix;
    solver::RadialSolver = RadialSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
    ngrid::Integer = solver.ngrid,
    rmax::Real = solver.rmax,
    nbasis::Integer = max(HO_DEFAULT_NBASIS, nlevels + 4),
)
    r, h = radial_grid(ngrid, rmax)
    size(V, 1) == length(r) || error("V must live on the (ngrid, rmax) mesh")
    best = nothing
    for β in oscillator_beta_candidates(masses, L)
        H, U = oscillator_hamiltonian_for_beta(params, masses, L, r, h, β; nbasis = nbasis)
        F = eigen(Symmetric(Matrix(H)))
        if isnothing(best) || F.values[nlevels] < best.values[nlevels]
            best = (beta = β, values = F.values, coeffs = F.vectors, basis = U)
        end
    end
    _warn_if_beta_railed(best.beta, oscillator_beta_candidates(masses, L), masses, L)
    waves = physically_normalized_waves(best.basis * best.coeffs, h)
    # V in the central eigenbasis: ⟨k|V|n⟩ = h · wₖ' V wₙ (physical normalization)
    Vkn = h .* (waves' * (Matrix(V) * waves))
    values = Float64[]
    distorted = Matrix{Float64}(undef, length(r), nlevels)
    for n = 1:nlevels
        ψ = copy(waves[:, n])
        for k in axes(waves, 2)
            k == n && continue
            denom = best.values[n] - best.values[k]
            abs(denom) < 1.0e-9 && continue
            ψ .+= (Vkn[k, n] / denom) .* waves[:, k]
        end
        ψ ./= sqrt(sum(abs2, ψ) * h)
        distorted[:, n] = ψ
        push!(values, best.values[n] + Vkn[n, n])
    end
    return values, distorted, r
end


"""
    ho_full_distorted_states(params, masses, L, V; ...)

The oscillator-basis method of [`resummed_channel_solution`](@ref), under its
original name. Kept because the Table VI and W6 audits call it directly.
"""
ho_full_distorted_states(
    params::GIParameters{HarmonicOscillatorBasis},
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix;
    kwargs...,
) = resummed_channel_solution(params, masses, L, V; kwargs...)
