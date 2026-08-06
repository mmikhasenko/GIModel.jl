# Paper-style harmonic-oscillator radial expansion — the method Godfrey & Isgur
# actually used (Appendix A, Eq. A17).
#
# On the Appendix-A central path there is NO spatial mesh in the operators:
#
#   momentum side   `ho_p2_matrix`        exact closed-form matrix elements
#   position side   `ho_operator_matrix`  generalized Gauss-Laguerre (Golub-Welsch)
#   A(p) factor     spectral function of the exact p²
#
# The uniform mesh survives only to reconstruct wavefunctions for reporting, and
# for the comparator central methods (pointwise, 1D/3D-smeared, derivative-G),
# several of which smear numerically ON the mesh and so are not closed-form
# functions of r. Those keep both sides on the mesh rather than becoming a
# hybrid: an exact kinetic operator with a mesh-projected potential is not the
# Hamiltonian of any single problem, and is not variational.

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
    ho_operator_matrix(L, β, nbasis, g; rtol=1e-12, nq=64, nq_max=4096) -> Symmetric

Matrix elements `⟨a|g(r)|b⟩ = ∫₀^∞ u_a(r) g(r) u_b(r) dr` in the 3D oscillator
basis — the **position-space** half of Eq. (A17), with no spatial mesh.

This is generalized Gauss–Laguerre quadrature in Golub–Welsch (DVR) form, and
the Jacobi matrix it needs is already [`ho_r2_matrix`](@ref): for weight
`x^(L+1/2) e^{-x}` with `x = (βr)²` the recurrence coefficients are exactly
`⟨n|r²|n⟩β²` and `⟨n|r²|n+1⟩β²`. So

    β² r²  =  Z diag(xᵢ) Zᵀ ,      rᵢ = √xᵢ / β
    ⟨a|g|b⟩ = Σᵢ Z[a,i] g(rᵢ) Z[b,i]

**Why this form and not weights-times-polynomials.** Evaluating `√wᵢ` and
`p_n(xᵢ)` separately is what destroys the accuracy: at 120 Gauss–Hermite nodes
the outer node reaches Laguerre argument ≈218, where `L₁₂^{1/2} ≈ 1.7e19` while
the weight is ≈1e-95. Their product is O(1), but forming it from those two
factors has already lost the digits before the sum starts — no amount of
compensated summation recovers them. The eigenvector entries **are** the product
`√wᵢ p_n(xᵢ)`; neither extreme ever exists.

Two exact self-checks, needing no reference data:

  - `g = 1` returns the identity (`Z` is orthogonal), to ~5e-15.
  - `g = r²` returns [`ho_r2_matrix`](@ref) — spectral reconstruction of the
    very Jacobi matrix that generated the rule — to ~1e-13.

`nq` (the quadrature size, distinct from `nbasis`) doubles until the result stops
moving by `rtol`. The default 1e-10 is well below what the eigenvalues need
(~1e-6 GeV) and is reached by every `(L, β)` in `HO_BETA_GRID`; the most diffuse
`β = 0.25` needs the most nodes, since a smooth non-polynomial `g` needs more nodes than a
polynomial one and diffuse `β` needs more than compact `β`. Verified against
`QuadGK` on the Appendix-A smeared potential: agreement 2.8e-16 to 1.2e-13.
"""
function ho_operator_matrix(
    L::Integer, β::Real, nbasis::Integer, g;
    rtol::Real = 1e-10, nq::Integer = 64, nq_max::Integer = 8192,
)
    nbasis >= 1 || throw(ArgumentError("ho_operator_matrix: nbasis must be ≥ 1"))
    β > 0 || throw(ArgumentError("ho_operator_matrix: β must be positive"))
    n = max(Int(nq), 2nbasis)
    prev = _ho_operator_matrix_at(L, β, nbasis, g, n)
    while n < nq_max
        n *= 2
        cur = _ho_operator_matrix_at(L, β, nbasis, g, n)
        scale = max(maximum(abs, cur), 1.0)
        maximum(abs, cur - prev) <= rtol * scale && return Symmetric(cur)
        prev = cur
    end
    @warn """
    ho_operator_matrix: quadrature did not reach rtol=$rtol by nq=$nq_max
    (L=$L, β=$β, nbasis=$nbasis). Diffuse β needs the most nodes; the result is
    the largest-nq value, not a converged one.
    """ maxlog = 1
    return Symmetric(prev)
end

function _ho_operator_matrix_at(L::Integer, β::Real, nbasis::Integer, g, nq::Integer)
    sqrt_x, Z = gauss_laguerre_dvr(L, nbasis, nq)
    r = sqrt_x ./ float(β)
    return Z * Diagonal([float(g(ri)) for ri in r]) * transpose(Z)
end

# The quadrature rule does not depend on β. In the dimensionless variable
# `x = (βr)²` the Jacobi matrix is `β² ho_r2_matrix(L, β, nq) = ho_r2_matrix(L, 1, nq)`
# — every β cancels — so the nodes `xᵢ` and the DVR matrix `Z` are functions of
# `(L, nq)` alone, and β enters only as the rescaling `rᵢ = √xᵢ / β`.
#
# That matters because `oscillator_channel_solution` scans 22 β candidates and
# assembles two operators (G̃ and S̃) at each, i.e. 44 requests for the same few
# decompositions per L. Memoizing on `(L, nbasis, nq)` turns the eigensolves from
# the dominant cost into a one-off; the stored slice is `nbasis × nq`, not
# `nq × nq`, so the whole cache is a few MB.
const _GAUSS_LAGUERRE_DVR = Dict{NTuple{3,Int},Tuple{Vector{Float64},Matrix{Float64}}}()

"""
    gauss_laguerre_dvr(L, nbasis, nq) -> (sqrt_x, Z)

Golub–Welsch data for the generalized Gauss–Laguerre rule with weight
`x^(L+1/2) e^{-x}`: `sqrt_x[i] = √xᵢ` at the `nq` nodes, and `Z = V[1:nbasis, :]`
the leading rows of the Jacobi eigenvectors, whose entries **are** the products
`√wᵢ p_n(xᵢ)`. See [`ho_operator_matrix`](@ref) for why that product must never
be formed from its two factors.

Depends on `β` not at all — see the note above the cache. Memoized, and the
memo is exact: same key, same `eigen` call, same bits.
"""
function gauss_laguerre_dvr(L::Integer, nbasis::Integer, nq::Integer)
    return get!(_GAUSS_LAGUERRE_DVR, (Int(L), Int(nbasis), Int(nq))) do
        # Keep the Jacobi matrix TRIDIAGONAL. Densifying it costs O(nq³) where
        # the tridiagonal solver is O(nq²), and the diffuse end of HO_BETA_GRID
        # needs nq ≈ 2048 — the difference between a fast gate and an unusable one.
        F = eigen(ho_r2_matrix(L, 1, nq))
        (sqrt.(max.(F.values, 0.0)), Matrix(F.vectors[1:min(nbasis, nq), :]))
    end
end

"""
    ho_r2_matrix(L, β, nbasis) -> SymTridiagonal

Exact matrix elements of `r²` in the 3D harmonic-oscillator basis — the
position-space companion of [`ho_p2_matrix`](@ref):

    ⟨n|r²|n⟩   =  (2n + L + 3/2) / β²
    ⟨n|r²|n+1⟩ = -√((n+1)(n + L + 3/2)) / β²

The two together are self-validating without any mesh, reference data or model
input: since `H = p²/2μ + ½μω²r²` is diagonal with eigenvalue `(2n+L+3/2)ω` and
`β² = μω`, the combination

    p²/2μ + (β⁴/2μ) r²

must come out **exactly diagonal**. Any error in either operator's magnitude,
sign or power of β breaks that cancellation, so the reconstruction tests both
and their relative normalization at once. See the `ho_r2_matrix / ho_p2_matrix`
testset.
"""
function ho_r2_matrix(L::Integer, β::Real, nbasis::Integer)
    nbasis >= 1 || throw(ArgumentError("ho_r2_matrix: nbasis must be ≥ 1"))
    β > 0 || throw(ArgumentError("ho_r2_matrix: β must be positive"))
    ib2 = 1 / float(β)^2
    diagonal = [ib2 * (2n + L + 1.5) for n = 0:(nbasis-1)]
    offdiag = [-ib2 * sqrt((n + 1) * (n + L + 1.5)) for n = 0:(nbasis-2)]
    return SymTridiagonal(diagonal, offdiag)
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

"""
    orthonormalize_physical_basis(U, h) -> Matrix

Orthonormalize mesh-sampled basis columns under the physical inner product
`integral f g dr = h * sum(f_i g_i)`, **with the phase of the input preserved**.

The phase fix is not cosmetic. LAPACK's QR assigns each column whatever sign its
algorithm produces, and those signs vary with `nbasis`, `β` and the mesh (at
`β = 0.65, nbasis = 24` they flipped at n = 10 and n = 18 — isolated, so not a
convention). That is invisible while everything is numerical, because a column
sign flip is a unitary transformation and leaves every eigenvalue and observable
alone. It becomes fatal the moment a **closed-form** matrix element is used, since
analytic formulas are written in `ho_reduced_radial`'s convention: mixing the two
put the oscillator charmonium 1S 9.4 MeV *below* the finite-difference answer,
which a variational calculation in a finite basis cannot do.

Same lesson as [`fix_annihilation_phase!`](@ref), which exists because eigenvector
signs are arbitrary and the Table III amplitudes depend on them. A basis that will
ever meet an analytic expression needs a stated phase.
"""
function orthonormalize_physical_basis(U::AbstractMatrix, h::Real)
    F = qr(sqrt(h) .* Matrix(U))
    cols = min(size(U, 2), size(F.Q, 2))
    Q = Matrix(F.Q[:, 1:cols]) ./ sqrt(h)
    for j in 1:cols
        # keep the sign of the basis function this column came from
        sum(view(Q, :, j) .* view(U, :, j)) < 0 && (Q[:, j] .*= -1)
    end
    return Q
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
    # `U` is built only to reconstruct waves back onto the reporting mesh; no
    # operator is assembled through it on the Appendix-A path below.
    U = orthonormalize_physical_basis(ho_basis_matrix(L, β, r, nbasis), h)
    if params.central isa AppendixAMomentumSandwich
        # Eq. (A17) with the paper's own ingredients, both sides exact and no
        # spatial mesh: p² has closed-form oscillator matrix elements, and the
        # smeared G̃/S̃ are integrated by generalized Gauss-Laguerre in
        # Golub-Welsch form. The momentum factor A(p) is a spectral function of
        # p², which is what makes the f(p) g(r) ordering of A17 assemble as a
        # matrix product.
        p2_basis = Symmetric(Matrix(ho_p2_matrix(L, β, nbasis)))
        kinetic = oscillator_kinetic_matrix(p2_basis, m1) +
                  oscillator_kinetic_matrix(p2_basis, m2)
        A = oscillator_momentum_factor_matrix(p2_basis, m1, m2; power = 0.5)
        g = ho_operator_matrix(L, β, nbasis, ri -> smeared_coulomb_G_closed(params, m1, m2, ri))
        s = ho_operator_matrix(L, β, nbasis, ri -> smeared_confinement_S_closed(params, m1, m2, ri))
        return Symmetric(kinetic + Symmetric(A * g * A + s)), U
    end
    # Comparator central methods (pointwise, 1D-smeared, 3D-smeared, derivative-G)
    # are not closed-form functions of r — several smear numerically ON the mesh —
    # so they keep the mesh projection. Mixing an exact kinetic operator with a
    # mesh-projected potential is not the Hamiltonian of any single problem, so
    # both sides stay on the mesh here.
    p2_grid = p2_operator(m1, L, r, h)
    p2_basis = projected_matrix(U, h, p2_grid)
    kinetic = oscillator_kinetic_matrix(p2_basis, m1) + oscillator_kinetic_matrix(p2_basis, m2)
    potential = projected_diagonal(U, h, potential_diagonal(params, m1, m2, r))
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
