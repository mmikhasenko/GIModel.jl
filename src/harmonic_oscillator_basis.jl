# Paper-style harmonic-oscillator radial expansion — the method Godfrey & Isgur
# actually used (Appendix A, Eq. A17).
#
# On the Appendix-A central path there is NO spatial mesh in the operators:
#
#   momentum side   `ho_p2_matrix`        exact closed-form matrix elements
#   position side   `ho_operator_matrix`  generalized Gauss-Laguerre (Golub-Welsch)
#   A(p) factor     spectral function of the exact p²
#
# A uniform mesh belongs only to the finite-difference comparator. Plotting code
# may explicitly sample an `OscillatorWave`, but no mesh is retained by the HO
# solver or accepted as an HO operator input.

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

**Why this form and not weights-times-polynomials.** Not because the huge×tiny
product loses digits — it does not. At 120 Gauss–Hermite nodes the outer node
reaches Laguerre argument ≈218, where `L₁₂^{1/2} ≈ 1.7e19` against a weight of
≈1e-95, but floating-point multiplication preserves relative accuracy and the
summands of the rule are all comparable in size (the scaled weight `wᵢ·e^{tᵢ²}`
stays between 0.29 and 0.74 across every node at 60), so there is no
cancellation for the digits to be lost to. Given correct weights, the
two-factor form holds 1.3e-15 on `g = 1` and 3e-13 on `g = r²` out to 200
nodes. Compensated summation is not needed and would not be the fix.

The actual hazard is narrower and lives entirely in the weights. Golub–Welsch
delivers them as `μ₀·v₁²`, the square of an eigenvector's first component. Once
the true weight is small enough that `v₁` reaches the eigensolver's noise floor,
LAPACK returns it as **exactly zero** — first zeros around a true weight of
1e-36, then 4 of 60 weights, 34 of 120, 88 of 200 — each one silently deleting a
whole node from the rule. Ordinary Gauss–Hermite users never notice, because
their integrand is negligible where those nodes sit; here it carries `e^{+tᵢ²}`
in the polynomial factor, so a deleted node costs a full-size contribution: the
`g = r²` error jumps from 3.3e-12 at 44 nodes to 2.0e-5 at 60.

The DVR form is immune because it never asks for that number. The eigenvector
entries **are** the product `√wᵢ p_n(xᵢ)`; neither extreme ever exists. Any
future rule that takes `√wᵢ` from an eigensolver inherits the trap.

Two exact self-checks, needing no reference data:

  - `g = 1` returns the identity (`Z` is orthogonal), to ~5e-15.
  - `g = r²` returns [`ho_r2_matrix`](@ref) — spectral reconstruction of the
    very Jacobi matrix that generated the rule — to ~1e-13.

`nq` (the quadrature size, distinct from `nbasis`) doubles until the result stops
moving by `rtol`. The default `1e-10` is well below what central-potential
eigenvalues need (~`1e-6` GeV) and is reached throughout `HO_BETA_GRID` for the
smooth Appendix-A central kernels. The very narrow Gaussian contact/tensor
kernels at deliberately diffuse endpoint bases use a separately certified
`1e-8`: after the A15 mass prefactor their remaining matrix uncertainty is
below the 0.1 MeV spectrum gate. Verified against `QuadGK` on the Appendix-A
smeared potential: agreement 2.8e-16 to 1.2e-13.
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
    return Z * Diagonal([g(ri) for ri in r]) * transpose(Z)
end

# The quadrature rule does not depend on β. In the dimensionless variable
# `x = (βr)²` the Jacobi matrix is `β² ho_r2_matrix(L, β, nq) = ho_r2_matrix(L, 1, nq)`
# — every β cancels — so the nodes `xᵢ` and the DVR matrix `Z` are functions of
# `(L, nq)` alone, and β enters only as the rescaling `rᵢ = √xᵢ / β`.
#
# That matters because the initial beta-bracket scan and its continuous
# refinement assemble two operators (G̃ and S̃) repeatedly for the same
# `(L, nbasis, nq)`. Memoizing the dimensionless rule makes those decompositions
# one-off work. Threading a rule store through five signatures instead would put
# this numerics-internal object in the public API.
#
# Capped by BYTES rather than entry count, because entries span `24 × 64` to
# `24 × 8192` — 256 of the small ones cost less than one of the large. The gate's
# own sweep needs 9.3 MiB; a user scanning `nbasis` (the convergence study
# `OscillatorSolver(nbasis = ...)` exists to invite) reached 67 MiB unbounded.
#
# The cap is a safety valve, not a policy: one β scan's working set is ~9 MiB
# against a 256 MiB cap, so it never fires in normal use, and when it does a
# flush costs time and nothing else. An LRU would be the better policy if the
# working set ever approached the cap — measured against a plain `Dict` here it
# was a wash (5.06 s vs 5.25 s on the 8-meson sweep, inside the noise), which is
# not worth a dependency.
const _GAUSS_LAGUERRE_DVR_MAX_BYTES = 256 * 2^20
const _GAUSS_LAGUERRE_DVR = Dict{NTuple{3,Int},Tuple{Vector{Float64},Matrix{Float64}}}()

_dvr_cache_bytes() =
    sum(sizeof(v[1]) + sizeof(v[2]) for v in values(_GAUSS_LAGUERRE_DVR); init = 0)

"""
    gauss_laguerre_dvr(L, nbasis, nq) -> (sqrt_x, Z)

Golub–Welsch data for the generalized Gauss–Laguerre rule with weight
`x^(L+1/2) e^{-x}`: `sqrt_x[i] = √xᵢ` at the `nq` nodes, and `Z = V[1:nbasis, :]`
the leading rows of the Jacobi eigenvectors, whose entries **are** the products
`√wᵢ p_n(xᵢ)`. Taking them this way is what keeps `wᵢ` itself from ever being
asked for — see [`ho_operator_matrix`](@ref) for why that matters, and for what
does and does not go wrong when it is.

Depends on `β` not at all — see the note above the cache. Memoized under a byte
cap, and the memo is exact: same key, same `eigen` call, same bits. Flushing can
only cost time, never accuracy.
"""
function gauss_laguerre_dvr(L::Integer, nbasis::Integer, nq::Integer)
    _dvr_cache_bytes() > _GAUSS_LAGUERRE_DVR_MAX_BYTES && empty!(_GAUSS_LAGUERRE_DVR)
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

Verified against the mesh projection it replaced: the difference falls from
4.3e-3 at `(ngrid, rmax) = (450, 24)` to 7.4e-5 at `(8000, 56)`, i.e. it is the
mesh's error and not this formula's.

**Wired into `oscillator_hamiltonian_for_beta`, but only atomically** — together
with an exact position side, never alone. On its own it produced a Hamiltonian
whose kinetic operator belonged to the continuum problem and whose potential
belonged to the discretized one, which is the Hamiltonian of no single problem
and not variational: the charmonium `1S` landed 9.4 MeV *below* the
finite-difference answer, the wrong side for a finite basis, and two Table VII
gluonic ratios fell out of band. Two prerequisites had to be found first. The
basis had no phase convention (see [`orthonormalize_physical_basis`](@ref)); an
earlier note here dismissed the QR signs as a uniform `-1` that cancels, which
was wrong — they flip on isolated columns (n = 10 and n = 18 at `β = 0.65,
nbasis = 24`), so they do not cancel. And the potential side was still projected
through the mesh, now [`ho_operator_matrix`](@ref).

With both sides exact the sign is right: charm sits +0.18 MeV and bottom
+0.55 MeV *above* the finite-difference result, as a variational calculation in
a finite basis must. The light sectors sit ≈1.5 MeV below, which reads the other
way — the oscillator answer is a true bound on the continuum one, so it is the
finite-difference mesh that is high where short-distance structure is hardest to
resolve. The finite-difference / oscillator tolerance widened from 1e-3 to 3e-3
at the same time, not from lost accuracy but because the two paths became
independent: while the oscillator path projected the finite-difference `p²` it
was a Galerkin restriction of that problem and inherited its discretization
error, so the two agreed artificially well.
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
    OscillatorWave(L, beta, coefficients)

A normalized radial wave expanded in harmonic-oscillator basis functions.
`coefficients[n+1]` multiplies `ho_reduced_radial(n, L, beta, r)`; `L` is the
orbital angular momentum and `beta` is a momentum scale in GeV. No spatial or
momentum mesh is stored.

A GIModel oscillator-basis solve determines the coefficients by diagonalizing
the model Hamiltonian. In that case this type represents a calculated model
wave, subject to basis convergence. Manually supplying `[1.0]` instead assumes
a single oscillator shape with a chosen width. The constructor normalizes the
supplied coefficients.
Prefer [`radial_wave`](@ref) to retrieve a solved wave.

## Example

```julia
using GIModel
assumed_wave = OscillatorWave(0, 0.5, [1.0]) # Assumed S-wave, not a GI solve.
@assert wave_norm(assumed_wave) ≈ 1
```

## Related

[`radial_wave`](@ref), [`RadialWave`](@ref), [`sample_wave`](@ref).
"""
struct OscillatorWave <: RadialWave
    L::Int
    beta::Float64
    coefficients::Vector{Float64}
    function OscillatorWave(
        L::Integer,
        beta::Real,
        coefficients::AbstractVector{<:Real},
    )
        L >= 0 || throw(ArgumentError("OscillatorWave: L must be nonnegative"))
        beta > 0 || throw(ArgumentError("OscillatorWave: beta must be positive"))
        isempty(coefficients) && throw(ArgumentError("OscillatorWave: empty coefficients"))
        coeffs = collect(Float64, coefficients)
        nrm = norm(coeffs)
        nrm > 0 || throw(ArgumentError("OscillatorWave: zero-norm coefficients"))
        return new(Int(L), float(beta), coeffs ./ nrm)
    end
end

wave_norm(w::OscillatorWave) = sum(abs2, w.coefficients)

function fix_outer_phase(w::OscillatorWave)
    # A tiny coefficient of the highest retained polynomial controls the strict
    # r -> infinity sign but not any physically occupied lobe. Locate the last
    # significant antinode of the analytic expansion on a dimensionless HO
    # interval instead. This fixes a phase; it does not discretize an operator
    # or turn the native wave into a mesh representation.
    rho_max = sqrt(4 * (length(w.coefficients) - 1) + 2w.L + 3) + 6
    values = [
        _oscillator_radial_value(w, rho / w.beta) for
        rho in range(0.0, rho_max; length = 2049)
    ]
    peak = maximum(abs, values)
    index = findlast(x -> abs(x) > 0.2peak, values)
    (isnothing(index) || values[index] >= 0) && return w
    return OscillatorWave(w.L, w.beta, -w.coefficients)
end

function _oscillator_radial_value(w::OscillatorWave, r::Real)
    return sum(
        w.coefficients[n + 1] * ho_reduced_radial(n, w.L, w.beta, r) for
        n in 0:(length(w.coefficients)-1)
    )
end

function _ho_reduced_radial_derivative(nr::Integer, L::Integer, beta::Real, r::Real)
    rho = beta * r
    x = rho^2
    alpha = L + 0.5
    norm = sqrt(2 * beta * gamma(nr + 1) / gamma(nr + L + 1.5))
    envelope = norm * rho^(L + 1) * exp(-0.5x)
    laguerre = generalized_laguerre(nr, alpha, x)
    # d L_n^alpha(x) / dx = -L_(n-1)^(alpha+1)(x).
    laguerre_derivative = iszero(nr) ? 0.0 :
                          -generalized_laguerre(nr - 1, alpha + 1, x)
    if iszero(r)
        # Only the L=0 basis has a nonzero derivative at the origin.
        return L == 0 ? norm * beta * generalized_laguerre(nr, alpha, 0.0) : 0.0
    end
    return envelope * (
        ((L + 1) / r - beta^2 * r) * laguerre +
        2beta^2 * r * laguerre_derivative
    )
end

function _oscillator_radial_derivative(w::OscillatorWave, r::Real)
    return sum(
        w.coefficients[n + 1] *
        _ho_reduced_radial_derivative(n, w.L, w.beta, r) for
        n in 0:(length(w.coefficients)-1)
    )
end

# Native HO waves are polynomial times a Gaussian. Integrating them through an
# infinite-interval variable transform needlessly evaluates the polynomial at
# enormous arguments, where high-order bases can form `Inf * 0 = NaN`. Ten
# dimensionless units beyond the classical turning radius suppress the omitted
# Gaussian tail far below the quadrature tolerances used by observables.
_oscillator_tail_rho(L::Integer, nbasis::Integer) =
    sqrt(4 * (nbasis - 1) + 2L + 3) + 10

_oscillator_coordinate_cutoff(w::OscillatorWave) =
    _oscillator_tail_rho(w.L, length(w.coefficients)) / w.beta

_oscillator_momentum_cutoff(w::OscillatorWave) =
    _oscillator_tail_rho(w.L, length(w.coefficients)) * w.beta

"""
    sample_wave(w::OscillatorWave, r) -> MeshWave

Sample an oscillator wave on a uniform radial grid `r` in GeV⁻¹ for plotting
or export. Supply at least two points and enough radial extent to contain the
wave. The returned samples are normalized on this finite grid, so this is a
plotting representation, not an independent convergence check.

## Example

```julia
using GIModel
wave = OscillatorWave(0, 0.5, [1.0])
sampled = sample_wave(wave, range(0.0, 15.0; length=301))
sampled.r, sampled.u  # plot radius against the reduced radial wave u(r)
```

A [`MeshWave`](@ref) already has `.r` and `.u`; it needs no sampling step.

## Related

[`RadialWave`](@ref), [`OscillatorWave`](@ref), [`radial_wave`](@ref),
[`physical_components`](@ref), [`wave_norm`](@ref).
"""
function sample_wave(w::OscillatorWave, r::AbstractVector{<:Real})
    samples = [_oscillator_radial_value(w, ri) for ri in r]
    h = r[2] - r[1]
    samples ./= sqrt(sum(abs2, samples) * h)
    return MeshWave(samples, r, h)
end

function radial_expect(w::OscillatorWave, f)
    op = ho_operator_matrix(w.L, w.beta, length(w.coefficients), f)
    return dot(w.coefficients, op * w.coefficients) / wave_norm(w)
end

function radial_overlap(left::OscillatorWave, right::OscillatorWave, f)
    nl, nr = sqrt(wave_norm(left)), sqrt(wave_norm(right))
    probe = f(inv(max(left.beta, right.beta)))
    if left.L == right.L && left.beta == right.beta && probe isa Real
        n = max(length(left.coefficients), length(right.coefficients))
        cl = vcat(left.coefficients, zeros(n - length(left.coefficients)))
        cr = vcat(right.coefficients, zeros(n - length(right.coefficients)))
        op = ho_operator_matrix(left.L, left.beta, n, f)
        return dot(cl, op * cr) / (nl * nr)
    end
    rmax = max(
        _oscillator_coordinate_cutoff(left),
        _oscillator_coordinate_cutoff(right),
    )
    value, _ = quadgk(
        r -> _oscillator_radial_value(left, r) * f(r) *
             _oscillator_radial_value(right, r),
        0.0,
        rmax;
        rtol = 1e-10,
    )
    return value / (nl * nr)
end


function radial_derivative_overlap(left::OscillatorWave, right::OscillatorWave, f)
    nl, nr = sqrt(wave_norm(left)), sqrt(wave_norm(right))
    rmax = max(
        _oscillator_coordinate_cutoff(left),
        _oscillator_coordinate_cutoff(right),
    )
    value, _ = quadgk(
        r -> _oscillator_radial_derivative(left, r) * f(r) *
             _oscillator_radial_value(right, r),
        0.0,
        rmax;
        rtol = 1e-10,
    )
    return value / (nl * nr)
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

function momentum_relativization_matrix(m1::Real, m2::Real, exponent::Real, p2_fact)
    lambda = max.(p2_fact.values, 0)
    e1 = sqrt.(lambda .+ m1^2)
    e2 = sqrt.(lambda .+ m2^2)
    diagonal = (m1 * m2 ./ (e1 .* e2)) .^ exponent
    return p2_fact.vectors * Diagonal(diagonal) * p2_fact.vectors'
end

gi_spin_dependent_side_exponent(epsilon::Real) = 0.5 + float(epsilon)

"""
    ho_momentum_sandwich_matrix(L, beta, nbasis, masses, epsilon, kernel)

Native harmonic-oscillator matrix for `B(p^2) kernel(r) B(p^2)`, where
`B=(m1*m2/(E1*E2))^(1/2+epsilon)`. No coordinate mesh is constructed.
"""
function ho_momentum_sandwich_matrix(
    L::Integer,
    beta::Real,
    nbasis::Integer,
    masses::ConstituentMasses,
    epsilon::Real,
    kernel;
    rtol::Real = 1e-10,
)
    p2 = Symmetric(Matrix(ho_p2_matrix(L, beta, nbasis)))
    B = momentum_relativization_matrix(
        masses.m1_GeV,
        masses.m2_GeV,
        gi_spin_dependent_side_exponent(epsilon),
        eigen(p2),
    )
    K = ho_operator_matrix(L, beta, nbasis, kernel; rtol = rtol)
    return Symmetric(B * K * B)
end

"""Native Eq. (A17) spin-independent oscillator Hamiltonian."""
function oscillator_central_matrix(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    beta::Real,
    nbasis::Integer,
)
    params.central isa AppendixAMomentumSandwich || throw(ArgumentError(
        "OscillatorSolver implements only the native Appendix-A momentum-sandwich " *
        "central Hamiltonian; use FiniteDifferenceSolver for " *
        "`$(central_potential_method(params.central))`",
    ))
    m1, m2 = masses.m1_GeV, masses.m2_GeV
    p2 = Symmetric(Matrix(ho_p2_matrix(L, beta, nbasis)))
    kinetic = oscillator_kinetic_matrix(p2, m1) + oscillator_kinetic_matrix(p2, m2)
    A = oscillator_momentum_factor_matrix(p2, m1, m2; power = 0.5)
    G = ho_operator_matrix(
        L, beta, nbasis, r -> smeared_coulomb_G_closed(params, m1, m2, r),
    )
    S = ho_operator_matrix(
        L, beta, nbasis, r -> smeared_confinement_S_closed(params, m1, m2, r),
    )
    return Symmetric(kinetic + Symmetric(A * G * A + S))
end

function oscillator_hamiltonian_for_beta(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    r::AbstractVector,
    h::Real,
    β::Real;
    nbasis::Integer = HO_DEFAULT_NBASIS,
)
    # Some convergence diagnostics still draw the analytic basis on a reporting
    # grid. The Hamiltonian itself is assembled entirely in coefficient space.
    U = orthonormalize_physical_basis(ho_basis_matrix(L, β, r, nbasis), h)
    return oscillator_central_matrix(params, masses, L, β, nbasis), U
end

const HO_REQUIRED_CONVERGED_REFINEMENTS = 2

_ho_objective(result) = last(result.values)

function _cached_ho_evaluation!(cache::Dict{Float64,T}, evaluate, beta::Real) where {T}
    beta_f = float(beta)
    return get!(cache, beta_f) do
        result = evaluate(beta_f)
        all(isfinite, result.values) || throw(ErrorException(
            "oscillator beta search produced non-finite eigenvalues at beta=$beta_f",
        ))
        result
    end
end

function _best_cached_ho_evaluation(cache)
    betas = sort!(collect(keys(cache)))
    best_beta = first(betas)
    best_result = cache[best_beta]
    for beta in Iterators.drop(betas, 1)
        result = cache[beta]
        if _ho_objective(result) < _ho_objective(best_result)
            best_beta, best_result = beta, result
        end
    end
    return (beta = best_beta, result = best_result)
end

function _local_beta_grid_index!(cache, evaluate, grid, seed)
    if isnothing(seed)
        for beta in grid
            _cached_ho_evaluation!(cache, evaluate, beta)
        end
        objectives = [_ho_objective(cache[beta]) for beta in grid]
        return argmin(objectives)
    end

    index = argmin(abs.(grid .- seed))
    index = clamp(index, 2, length(grid) - 1)
    while true
        for i in (index - 1):min(index + 1, length(grid))
            _cached_ho_evaluation!(cache, evaluate, grid[i])
        end
        center = _ho_objective(cache[grid[index]])
        left = _ho_objective(cache[grid[index - 1]])
        right = _ho_objective(cache[grid[index + 1]])
        if left < center && left <= right
            index -= 1
        elseif right < center && right < left
            index += 1
        else
            return index
        end
        1 < index < length(grid) || return index
    end
end

function _refine_oscillator_beta(evaluate, solver::OscillatorSolver; seed = nothing)
    grid = solver.beta_grid
    first_beta = isnothing(seed) ? first(grid) : clamp(float(seed), first(grid), last(grid))
    first_result = evaluate(first_beta)
    all(isfinite, first_result.values) || throw(ErrorException(
        "oscillator beta search produced non-finite eigenvalues at beta=$first_beta",
    ))
    cache = Dict{Float64,typeof(first_result)}(first_beta => first_result)
    if length(grid) == 1
        return (beta = first_beta, result = first_result)
    end

    index = _local_beta_grid_index!(cache, evaluate, grid, seed)
    if index == 1 || index == length(grid)
        edge = index == 1 ? "lower" : "upper"
        throw(ErrorException(
            "oscillator beta optimum reached the $edge beta_grid endpoint " *
            "($(grid[index]) GeV); widen the declared beta bracket",
        ))
    end

    a, b = grid[index - 1], grid[index + 1]
    inverse_phi = (sqrt(5.0) - 1.0) / 2.0
    c = b - inverse_phi * (b - a)
    d = a + inverse_phi * (b - a)
    fc = _ho_objective(_cached_ho_evaluation!(cache, evaluate, c))
    fd = _ho_objective(_cached_ho_evaluation!(cache, evaluate, d))
    while b - a > solver.beta_tolerance_GeV
        if fc <= fd
            b, d, fd = d, c, fc
            c = b - inverse_phi * (b - a)
            fc = _ho_objective(_cached_ho_evaluation!(cache, evaluate, c))
        else
            a, c, fc = c, d, fd
            d = a + inverse_phi * (b - a)
            fd = _ho_objective(_cached_ho_evaluation!(cache, evaluate, d))
        end
    end
    return _best_cached_ho_evaluation(cache)
end

function _oscillator_solution_search(
    solver::OscillatorSolver,
    L::Integer,
    nlevels::Integer,
    evaluate,
)
    nlevels >= 1 || throw(ArgumentError(
        "oscillator solve requires at least one eigenlevel",
    ))
    nlevels <= solver.nbasis || throw(ArgumentError(
        "oscillator solve requested $nlevels levels from an initial basis of " *
        "$(solver.nbasis); raise nbasis to at least nlevels",
    ))
    initial_nbasis = solver.nbasis
    initial_nbasis <= solver.max_nbasis || throw(ErrorException(
        "oscillator convergence needs at least nbasis=$initial_nbasis for $nlevels levels, " *
        "above max_nbasis=$(solver.max_nbasis)",
    ))

    beta_seed = nothing
    previous_values = nothing
    previous_waves = nothing
    consecutive = 0
    refinements = 0
    nbasis = initial_nbasis
    best = nothing
    while true
        best = _refine_oscillator_beta(
            beta -> evaluate(beta, nbasis), solver; seed = beta_seed,
        )
        waves = [
            fix_outer_phase(OscillatorWave(L, best.beta, view(best.result.vectors, :, n)))
            for n in 1:nlevels
        ]
        if !solver.converge
            certificate = OscillatorConvergence(
                :unchecked, best.beta, nbasis, nothing, nothing,
                solver.energy_tolerance_GeV, 0,
            )
            return best, waves, certificate
        end

        if !isnothing(previous_values)
            refinements += 1
            energy_delta = maximum(abs.(best.result.values .- previous_values))
            overlap_defect = maximum(
                max(0.0, 1.0 - min(1.0, abs(radial_overlap(
                    previous_waves[n], waves[n], _ -> 1.0,
                )))) for n in 1:nlevels
            )
            consecutive = energy_delta <= solver.energy_tolerance_GeV ?
                          consecutive + 1 : 0
            if consecutive >= HO_REQUIRED_CONVERGED_REFINEMENTS
                certificate = OscillatorConvergence(
                    :converged,
                    best.beta,
                    nbasis,
                    energy_delta,
                    overlap_defect,
                    solver.energy_tolerance_GeV,
                    refinements,
                )
                return best, waves, certificate
            end
        end

        nbasis == solver.max_nbasis && break
        previous_values = copy(best.result.values)
        previous_waves = waves
        beta_seed = best.beta
        nbasis = min(nbasis + solver.basis_step, solver.max_nbasis)
    end
    final_delta = isnothing(previous_values) || isnothing(best) ? "not measured" :
                  "$(1000 * maximum(abs.(best.result.values .- previous_values))) MeV"
    throw(ErrorException(
        "oscillator basis did not converge through max_nbasis=$(solver.max_nbasis); " *
        "final requested-level change was $final_delta, tolerance is " *
        "$(1000 * solver.energy_tolerance_GeV) MeV",
    ))
end

function oscillator_channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer;
    solver::OscillatorSolver = OscillatorSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    evaluate = function (β, basis_size)
        H = oscillator_central_matrix(params, masses, L, β, basis_size)
        vals, vecs = lowest_eigenpairs(Matrix(H), nlevels, solver)
        return (values = vals, vectors = vecs)
    end
    best, waves, convergence = _oscillator_solution_search(
        solver, L, nlevels, evaluate,
    )
    return ChannelRadialSolution(
        best.result.values, waves; convergence = convergence,
    )
end

# A matrix without representation metadata is a finite-difference mesh
# operator. HO deliberately has no adapter for it.
function _resummed_channel_solution(
    solver::OscillatorSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix,
    nlevels::Integer,
)
    throw(ArgumentError(
        "OscillatorSolver cannot consume a mesh operator; use " *
        "fixed_channel_solution to assemble the native fixed-(L,S,J) Hamiltonian",
    ))
end
