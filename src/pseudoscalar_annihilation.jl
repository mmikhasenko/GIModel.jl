# Isoscalar pseudoscalar annihilation diagnostics and mixing blocks.
#
# Public API (exported from GIModel.jl):
#   isoscalar_pseudoscalar_annihilation_solution

const PSEUDOSCALAR_PERTURBATIVE_COEFF = 2π / 3 * (log(2) - 1)

abstract type PseudoscalarAnnihilationModel end

abstract type PseudoscalarSmearingScheme end

"""Legacy FD proxy: coordinate-space origin times one averaged relativistic factor."""
struct FDOriginP2Smearing <: PseudoscalarSmearingScheme end

"""FD evaluation of the Eq. (17) S-wave momentum integral."""
struct FDMomentumIntegralSmearing <: PseudoscalarSmearingScheme
    npoints::Int
end
FDMomentumIntegralSmearing() = FDMomentumIntegralSmearing(900)

"""Calibrated rank-one control matching the digitized GI Fig. 5/Table III masses."""
struct CalibratedP1Annihilation <: PseudoscalarAnnihilationModel end

"""Paper Eq. (18a) pseudoscalar annihilation model."""
struct PaperP1Annihilation{S<:PseudoscalarSmearingScheme} <: PseudoscalarAnnihilationModel
    smearing::S
end
PaperP1Annihilation() = PaperP1Annihilation(FDMomentumIntegralSmearing())

"""Paper Eq. (18b) mass-dependent pseudoscalar annihilation model."""
struct PaperP2Annihilation{S<:PseudoscalarSmearingScheme} <: PseudoscalarAnnihilationModel
    smearing::S
end
PaperP2Annihilation() = PaperP2Annihilation(FDMomentumIntegralSmearing())

struct PseudoscalarAnnihilationBasisInput
    label::String
    constituent_mass_GeV::Float64
    diagonal_GeV::Float64
    radial::RadialWaveOnUniformMesh
end

function pseudoscalar_annihilation_basis_input(
    label::AbstractString,
    constituent_mass_GeV::Real,
    diagonal_GeV::Real,
    radial::RadialWaveOnUniformMesh,
)
    return PseudoscalarAnnihilationBasisInput(
        String(label),
        float(constituent_mass_GeV),
        float(diagonal_GeV),
        radial,
    )
end

function _rank_one_annihilation_weights(
    diagonal::AbstractVector{<:Real},
    target::AbstractVector{<:Real},
)
    n = length(diagonal)
    length(target) == n ||
        throw(ArgumentError("diagonal and target sizes differ"))
    weights = Float64[]
    for i in 1:n
        num = prod(target[j] - diagonal[i] for j in 1:n)
        den = prod(diagonal[m] - diagonal[i] for m in 1:n if m != i)
        push!(weights, num / den)
    end
    return weights
end

function _phase_fix_columns!(vectors::AbstractMatrix{<:Real}; anchor::Integer = 1)
    for col in axes(vectors, 2)
        if vectors[anchor, col] < 0
            vectors[:, col] .*= -1
        end
    end
    return vectors
end

function _phase_fix_by_largest_component!(vectors::AbstractMatrix{<:Real})
    for col in axes(vectors, 2)
        anchor = argmax(abs.(vectors[:, col]))
        if vectors[anchor, col] < 0
            vectors[:, col] .*= -1
        end
    end
    return vectors
end

function _alpha_s_mass_scale(mass_GeV::Real)
    q = max(float(mass_GeV), 1.0e-9)
    return alpha_s_q(q)
end

function _reduced_p2_expectation(radial::RadialWaveOnUniformMesh)
    u = radial.u
    h = radial.h
    norm = max(sum(abs2, u) * h, eps(Float64))
    prev = 0.0
    accum = 0.0
    for ui in u
        accum += (ui - prev)^2 / h
        prev = ui
    end
    accum += prev^2 / h
    return accum / norm
end

function _s0_smearing_factor(::FDOriginP2Smearing, input::PseudoscalarAnnihilationBasisInput)
    wave = input.radial
    origin_R = abs(wave.u[1] / wave.r[1])
    p2 = _reduced_p2_expectation(wave)
    rel = input.constituent_mass_GeV / sqrt(input.constituent_mass_GeV^2 + p2)
    return origin_R * rel / sqrt(4π)
end

function _j0(x::Real)
    abs(x) < 1.0e-8 && return 1.0 - x^2 / 6
    return sin(x) / x
end

function _spherical_bessel_j(L::Integer, x::Real)
    L == 0 && return _j0(x)
    if abs(x) < 1.0e-4
        # j_L(x) ≈ x^L / (2L+1)!! for small argument.
        dfact = prod(1:2:(2L+1))
        return x^L / dfact
    end
    jm1 = _j0(x)
    j = sin(x) / x^2 - cos(x) / x
    for l in 1:(L-1)
        jp1 = (2l + 1) / x * j - jm1
        jm1 = j
        j = jp1
    end
    return j
end

function _momentum_radial_wave(radial::RadialWaveOnUniformMesh, p::Real, L::Integer)
    accum = 0.0
    for k in eachindex(radial.r)
        accum += radial.r[k] * radial.u[k] * _spherical_bessel_j(L, p * radial.r[k])
    end
    return sqrt(2 / π) * accum * radial.h
end

_momentum_radial_swave(radial::RadialWaveOnUniformMesh, p::Real) =
    _momentum_radial_wave(radial, p, 0)

"""
    _sL_smearing_factor(scheme, input, L)

Eq. (17) smearing of the wavefunction at the origin for orbital `L`:
`S_L(Ψ) = (2π)^(-3/2) ∫ d³p (4π)^(-1/2) Φ(p) [p/E]^L (m/E)` with `Φ(p)` the
normalized radial momentum wavefunction obtained from the `j_L` transform.
"""
function _sL_smearing_factor(
    scheme::FDMomentumIntegralSmearing,
    input::PseudoscalarAnnihilationBasisInput,
    L::Integer,
)
    wave = input.radial
    mass = input.constituent_mass_GeV
    npoints = max(scheme.npoints, 32)
    pmax = π / wave.h
    dp = pmax / (npoints - 1)
    accum = 0.0
    for k in 1:npoints
        p = (k - 1) * dp
        weight = (k == 1 || k == npoints) ? 0.5 : 1.0
        Φ = _momentum_radial_wave(wave, p, L)
        E = sqrt(mass^2 + p^2)
        accum += weight * p^2 * Φ * (p / E)^L * mass / E
    end
    return sqrt(2 / π) * accum * dp / sqrt(4π)
end

_s0_smearing_factor(scheme::FDMomentumIntegralSmearing, input::PseudoscalarAnnihilationBasisInput) =
    _sL_smearing_factor(scheme, input, 0)

"""
    fix_annihilation_phase!(vecs, r)

Enforce the GI annihilation phase convention on radial eigenvector columns:
the momentum-space wave at the origin, `Φ(0) ∝ ∫ r u(r) dr`, is positive for
every level. For nodeless ground states this coincides with `u(r_min) > 0`;
for radially excited states the outer lobe dominates the integral, so the
sign alternates relative to the small-`r` convention. Table III amplitude
signs follow this convention: with it, all four literal-P1 pseudoscalar
eigenvectors match the published sign structure, while the small-`r`
convention flips the `2 ns`/`2 ss` couplings and doubles the amplitude RMS.
"""
function fix_annihilation_phase!(vecs::AbstractMatrix{<:Real}, r::AbstractVector{<:Real})
    for col in eachcol(vecs)
        accum = 0.0
        for k in eachindex(r)
            accum += r[k] * col[k]
        end
        accum < 0 && (col .*= -1)
    end
    return vecs
end

function _flavor_coherence_factor(input::PseudoscalarAnnihilationBasisInput)
    label = lowercase(input.label)
    return (occursin("ns", label) || occursin("n nbar", label)) ? sqrt(2.0) : 1.0
end

function _annihilation_overlap_factor(
    scheme::PseudoscalarSmearingScheme,
    input::PseudoscalarAnnihilationBasisInput,
)
    return _flavor_coherence_factor(input) * _s0_smearing_factor(scheme, input)
end

function _annihilation_overlap_factor(
    scheme::FDMomentumIntegralSmearing,
    input::PseudoscalarAnnihilationBasisInput,
    L::Integer,
)
    return _flavor_coherence_factor(input) * _sL_smearing_factor(scheme, input, L)
end

function _paper_p1_bracket(
    params::GIParameters,
    left::PseudoscalarAnnihilationBasisInput,
    right::PseudoscalarAnnihilationBasisInput,
)
    mi = left.constituent_mass_GeV
    mj = right.constituent_mass_GeV
    alpha_i = _alpha_s_mass_scale(left.diagonal_GeV)
    alpha_j = _alpha_s_mass_scale(right.diagonal_GeV)
    nonperturbative = params.annihilation_p1_A_np *
                      exp(-(mi^2 + mj^2) / params.annihilation_p1_m_eta^2)
    perturbative = PSEUDOSCALAR_PERTURBATIVE_COEFF * alpha_i * alpha_j / π^2
    return nonperturbative + perturbative
end

function _paper_p2_bracket(
    params::GIParameters,
    left::PseudoscalarAnnihilationBasisInput,
    right::PseudoscalarAnnihilationBasisInput,
    pole_mass_GeV::Real,
)
    mi = left.constituent_mass_GeV
    mj = right.constituent_mass_GeV
    M = max(float(pole_mass_GeV), 1.0e-9)
    M0 = params.annihilation_p2_M0
    alpha = _alpha_s_mass_scale(M)
    nonperturbative =
        params.annihilation_p2_A_np *
        (1 - (M / M0)^4) *
        exp(-(mi^2 + mj^2) / M0^2 - M^4 / (4M0^4))
    perturbative = PSEUDOSCALAR_PERTURBATIVE_COEFF * (alpha / π)^2
    return nonperturbative + perturbative
end

function _paper_annihilation_matrix(
    model::PaperP1Annihilation,
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput},
)
    n = length(basis)
    matrix = Matrix(Diagonal([state.diagonal_GeV for state in basis]))
    factors = [_annihilation_overlap_factor(model.smearing, state) for state in basis]
    for j in 1:n, i in 1:n
        mi = basis[i].constituent_mass_GeV
        mj = basis[j].constituent_mass_GeV
        matrix[j, i] += 4π *
                        _paper_p1_bracket(params, basis[j], basis[i]) *
                        factors[j] *
                        factors[i] / (mj * mi)
    end
    return matrix
end

function _paper_annihilation_matrix(
    model::PaperP2Annihilation,
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput};
    pole_mass_GeV::Real,
)
    n = length(basis)
    matrix = Matrix(Diagonal([state.diagonal_GeV for state in basis]))
    factors = [_annihilation_overlap_factor(model.smearing, state) for state in basis]
    for j in 1:n, i in 1:n
        mi = basis[i].constituent_mass_GeV
        mj = basis[j].constituent_mass_GeV
        matrix[j, i] += 4π *
                        _paper_p2_bracket(params, basis[j], basis[i], pole_mass_GeV) *
                        factors[j] *
                        factors[i] / (mj * mi)
    end
    return matrix
end

function _basis_states_from_inputs(basis::AbstractVector{PseudoscalarAnnihilationBasisInput})
    return [
        BasisState(i <= 2 ? 1 : 2, "S", 1, 0; label = state.label) for
        (i, state) in enumerate(basis)
    ]
end

"""
    isoscalar_pseudoscalar_annihilation_solution(diagonal; targets)

Build the calibrated rank-one `^1S_0` isoscalar annihilation block over the
`[1 n nbar, 1 s sbar, 2 n nbar, 2 s sbar]` basis, targeting the four masses in
`targets` (GeV). The digitized GI Fig. 5 values live in the GIPaper comparison
package — this function only performs the calibration. It deliberately labels
the operation as an annihilation block rather than changing contact hyperfine
or central potential parameters.
"""
function isoscalar_pseudoscalar_annihilation_solution(
    ::CalibratedP1Annihilation,
    diagonal::AbstractVector{<:Real};
    targets,
)
    diag = collect(Float64, diagonal)
    targ = sort(collect(Float64, targets))
    weights = _rank_one_annihilation_weights(diag, targ)
    if any(w -> w < -1e-10, weights)
        throw(ArgumentError(
            "target pseudoscalar masses do not form a positive rank-one annihilation update",
        ))
    end
    weights = max.(weights, 0.0)
    couplings = sqrt.(weights)
    annihilation = couplings * transpose(couplings)
    matrix = Matrix(Diagonal(diag)) + annihilation
    fact = eigen(Symmetric(matrix))
    vectors = _phase_fix_columns!(Matrix(fact.vectors))
    basis = [
        BasisState(1, "S", 1, 0; label = "1 n nbar"),
        BasisState(1, "S", 1, 0; label = "1 s sbar"),
        BasisState(2, "S", 1, 0; label = "2 n nbar"),
        BasisState(2, "S", 1, 0; label = "2 s sbar"),
    ]
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        basis,
        matrix;
        mechanism = "rank_one_calibrated_pseudoscalar_annihilation",
        source = "GI Fig. 5 / Table III P1 diagnostic",
        notes = "Calibrated to digitized GI isoscalar pseudoscalar masses; separate PaperP1Annihilation and PaperP2Annihilation methods implement Eq. (18a,b).",
    )
    return (
        block = block,
        masses = collect(Float64, fact.values),
        vectors = vectors,
        diagonal_GeV = diag,
        targets_GeV = targ,
        couplings_sqrt_GeV = couplings,
        weights_GeV = weights,
        annihilation_matrix_GeV = annihilation,
    )
end

"""
    isoscalar_general_annihilation_solution(params, basis; amplitude_A, L, multiplicity, J, smearing)

General non-pseudoscalar Eq. (16) annihilation block in the channel
`^{2S+1}L_J`. The matrix element from `q_i q̄_i → q_j q̄_j` is

`4π(2L+1) A(^{2S+1}L_J) [α_s(M_j²)α_s(M_i²)/π²]^{n/2} S_L(Ψ_j)S_L(Ψ_i)/(m_i m_j)`

with `n = 2` for `C = +` and `n = 3` for `C = −` (`C = (−1)^{L+S}`) and the
Eq. (17) smearing `S_L`. `basis` should be built with HO radial wavefunctions
for paper-consistent wavefunction-at-origin scale.
"""
function isoscalar_general_annihilation_solution(
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput};
    amplitude_A::Real,
    L::Integer,
    multiplicity::Integer,
    J::Integer,
    smearing::FDMomentumIntegralSmearing = FDMomentumIntegralSmearing(),
)
    n = length(basis)
    S = multiplicity == 3 ? 1 : 0
    c_even = iseven(L + S)
    n_gluons = c_even ? 2 : 3
    matrix = Matrix(Diagonal([state.diagonal_GeV for state in basis]))
    factors = [_annihilation_overlap_factor(smearing, state, L) for state in basis]
    for j in 1:n, i in 1:n
        mi = basis[i].constituent_mass_GeV
        mj = basis[j].constituent_mass_GeV
        alpha_i = alpha_s_q(basis[i].diagonal_GeV)
        alpha_j = alpha_s_q(basis[j].diagonal_GeV)
        matrix[j, i] += 4π * (2L + 1) * amplitude_A *
                         (alpha_i * alpha_j / π^2)^(n_gluons / 2) *
                         factors[j] * factors[i] / (mj * mi)
    end
    fact = eigen(Symmetric(matrix))
    vectors = _phase_fix_by_largest_component!(Matrix(fact.vectors))
    basis_states = [BasisState(1, "S", multiplicity, J; label = state.label) for state in basis]
    channel = @sprintf("^%d%s_%d", multiplicity, L_LABELS[L], J)
    block = MixingBlock(
        "isoscalar $channel annihilation",
        basis_states,
        matrix;
        mechanism = "general_eq16_annihilation",
        source = "GI Eq. (16) with Table III A($channel)=$(amplitude_A), n=$(n_gluons) gluons",
        notes = "Uses HO radial waves for paper-consistent wavefunction-at-origin scale.",
    )
    return (
        block = block,
        masses = collect(Float64, fact.values),
        vectors = vectors,
        diagonal_GeV = [state.diagonal_GeV for state in basis],
        annihilation_matrix_GeV = matrix - Diagonal([state.diagonal_GeV for state in basis]),
    )
end

"""
    isoscalar_general_s1_solution(params, basis; smearing)

`^3S_1` special case of [`isoscalar_general_annihilation_solution`](@ref) with
the Table III amplitude `A(^3S_1)` from the parameters (three-gluon bracket,
`C = −`).
"""
function isoscalar_general_s1_solution(
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput};
    smearing::FDMomentumIntegralSmearing = FDMomentumIntegralSmearing(),
)
    return isoscalar_general_annihilation_solution(
        params,
        basis;
        amplitude_A = params.annihilation_s1_A,
        L = 0,
        multiplicity = 3,
        J = 1,
        smearing = smearing,
    )
end

function isoscalar_pseudoscalar_annihilation_solution(
    diagonal::AbstractVector{<:Real};
    targets,
)
    return isoscalar_pseudoscalar_annihilation_solution(
        CalibratedP1Annihilation(),
        diagonal;
        targets = targets,
    )
end

function isoscalar_pseudoscalar_annihilation_solution(
    model::PaperP1Annihilation,
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput},
)
    matrix = _paper_annihilation_matrix(model, params, basis)
    fact = eigen(Symmetric(matrix))
    vectors = _phase_fix_by_largest_component!(Matrix(fact.vectors))
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        _basis_states_from_inputs(basis),
        matrix;
        mechanism = "paper_p1_pseudoscalar_annihilation",
        source = "GI Eq. (16) with Eq. (18a) and Table III P1 constants",
        notes = "Uses cached FD radial waves in the Eq. (17) S-wave momentum integral.",
    )
    return (
        block = block,
        masses = collect(Float64, fact.values),
        vectors = vectors,
        diagonal_GeV = [state.diagonal_GeV for state in basis],
        annihilation_matrix_GeV = matrix - Diagonal([state.diagonal_GeV for state in basis]),
    )
end

function isoscalar_pseudoscalar_annihilation_solution(
    model::PaperP2Annihilation,
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput};
    maxiter::Integer = 80,
    tol::Real = 1.0e-10,
)
    n = length(basis)
    starts = sort([state.diagonal_GeV for state in basis])
    masses = zeros(Float64, n)
    vectors = zeros(Float64, n, n)
    matrices = Vector{Matrix{Float64}}(undef, n)
    for level in 1:n
        M = starts[level]
        matrix = _paper_annihilation_matrix(model, params, basis; pole_mass_GeV = M)
        fact = eigen(Symmetric(matrix))
        for _ in 1:maxiter
            newM = fact.values[level]
            abs(newM - M) <= tol && break
            M = 0.5 * (M + newM)
            matrix = _paper_annihilation_matrix(model, params, basis; pole_mass_GeV = M)
            fact = eigen(Symmetric(matrix))
        end
        masses[level] = fact.values[level]
        vectors[:, level] .= fact.vectors[:, level]
        matrices[level] = matrix
    end
    _phase_fix_by_largest_component!(vectors)
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        _basis_states_from_inputs(basis),
        matrices[1];
        mechanism = "paper_p2_pseudoscalar_annihilation",
        source = "GI Eq. (16) with Eq. (18b) and Table III P2 constants",
        notes = "Each pole is solved as a fixed point of the mass-dependent Eq. (18b) matrix using the Eq. (17) S-wave momentum integral, so eigenvectors are not expected to be orthogonal.",
    )
    return (
        block = block,
        masses = masses,
        vectors = vectors,
        diagonal_GeV = [state.diagonal_GeV for state in basis],
        annihilation_matrices_GeV = matrices,
    )
end
