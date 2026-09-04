# Isoscalar pseudoscalar annihilation diagnostics and mixing blocks.
#
# Public API (exported from GIModel.jl):
#   isoscalar_pseudoscalar_annihilation_solution

const PSEUDOSCALAR_PERTURBATIVE_COEFF = 2π / 3 * (log(2) - 1)

abstract type PseudoscalarAnnihilationModel end

abstract type PseudoscalarSmearingScheme end

"""Numerical evaluation of the Eq. (17) momentum integral."""
struct MomentumIntegralSmearing <: PseudoscalarSmearingScheme
    npoints::Int
end
MomentumIntegralSmearing() = MomentumIntegralSmearing(900)

"""Calibrated rank-one control matching the digitized GI Fig. 5/Table III masses."""
struct CalibratedP1Annihilation <: PseudoscalarAnnihilationModel end

"""Paper Eq. (18a) pseudoscalar annihilation model."""
struct PaperP1Annihilation{S<:PseudoscalarSmearingScheme} <: PseudoscalarAnnihilationModel
    smearing::S
end
PaperP1Annihilation() = PaperP1Annihilation(MomentumIntegralSmearing())

"""Paper Eq. (18b) mass-dependent pseudoscalar annihilation model."""
struct PaperP2Annihilation{S<:PseudoscalarSmearingScheme} <: PseudoscalarAnnihilationModel
    smearing::S
end
PaperP2Annihilation() = PaperP2Annihilation(MomentumIntegralSmearing())

struct PseudoscalarAnnihilationBasisInput
    basis::BasisState
    constituent_mass_GeV::Float64
    diagonal_GeV::Float64
    radial_components::Vector{Tuple{Float64,RadialWave}}
    # Whether this channel is the coherent (u ubar + d dbar)/sqrt(2) combination,
    # which carries a sqrt(2) amplitude factor relative to a single flavor. This
    # is a property of the flavor state, so it is stated, not guessed.
    isoscalar_coherent::Bool
end

Base.getproperty(input::PseudoscalarAnnihilationBasisInput, name::Symbol) =
    name === :label ? getfield(input, :basis).label : getfield(input, name)
Base.propertynames(::PseudoscalarAnnihilationBasisInput) = (
    :basis, :label, :constituent_mass_GeV, :diagonal_GeV,
    :radial_components, :isoscalar_coherent,
)

function pseudoscalar_annihilation_basis_input(
    basis::BasisState,
    constituent_mass_GeV::Real,
    diagonal_GeV::Real,
    radial::RadialWave;
    isoscalar_coherent::Bool,
)
    return PseudoscalarAnnihilationBasisInput(
        basis,
        float(constituent_mass_GeV),
        float(diagonal_GeV),
        Tuple{Float64,RadialWave}[(1.0, radial)],
        isoscalar_coherent,
    )
end

function pseudoscalar_annihilation_basis_input(
    basis::BasisState,
    constituent_mass_GeV::Real,
    diagonal_GeV::Real,
    radial_components::AbstractVector{<:Tuple};
    isoscalar_coherent::Bool,
)
    isempty(radial_components) && throw(ArgumentError(
        "pseudoscalar annihilation input needs at least one radial component",
    ))
    components = Tuple{Float64,RadialWave}[
        (float(coefficient), wave) for (coefficient, wave) in radial_components
    ]
    return PseudoscalarAnnihilationBasisInput(
        basis, float(constituent_mass_GeV), float(diagonal_GeV),
        components, isoscalar_coherent,
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

function _alpha_s_mass_scale(mass_GeV::Real)
    q = max(float(mass_GeV), 1.0e-9)
    return alpha_s_q(q)
end

"""
    _sL_smearing_factor(scheme, input, L)

Eq. (17) smearing of the wavefunction at the origin for orbital `L`:
`S_L(Ψ) = (2π)^(-3/2) ∫ d³p (4π)^(-1/2) Φ(p) [p/E]^L (m/E)` with `Φ(p)` the
normalized radial momentum wavefunction obtained from the `j_L` transform.
"""
function _sL_smearing_factor(
    scheme::MomentumIntegralSmearing,
    input::PseudoscalarAnnihilationBasisInput,
    L::Integer,
)
    mass = input.constituent_mass_GeV
    kernel = p -> begin
        E = sqrt(mass^2 + p^2)
        (p / E)^L * mass / E
    end
    factor = sum(input.radial_components) do (coefficient, radial)
        mw = radial isa MeshWave ?
             momentum_wave(radial, L; pmax = π / radial.h,
                           npoints = max(scheme.npoints, 32)) :
             momentum_wave(radial, L)
        coefficient * momentum_functional(mw, kernel)
    end
    return sqrt(2 / π) / sqrt(4π) * factor
end

_s0_smearing_factor(scheme::MomentumIntegralSmearing, input::PseudoscalarAnnihilationBasisInput) =
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

_flavor_coherence_factor(input::PseudoscalarAnnihilationBasisInput) =
    input.isoscalar_coherent ? sqrt(2.0) : 1.0

function _annihilation_overlap_factor(
    scheme::PseudoscalarSmearingScheme,
    input::PseudoscalarAnnihilationBasisInput,
)
    return _flavor_coherence_factor(input) * _s0_smearing_factor(scheme, input)
end

function _annihilation_overlap_factor(
    scheme::MomentumIntegralSmearing,
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
    nonperturbative = params.annihilation.p1_A_np *
                      exp(-(mi^2 + mj^2) / params.annihilation.p1_m_eta^2)
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
    M0 = params.annihilation.p2_M0
    alpha = _alpha_s_mass_scale(M)
    nonperturbative =
        params.annihilation.p2_A_np *
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
    return [state.basis for state in basis]
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
    basis = nothing,
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
    basis_states = isnothing(basis) ? [
        BasisState(1, "S", 1, 0; label = "1 n nbar", flavors = (:q, :q)),
        BasisState(1, "S", 1, 0; label = "1 s sbar", flavors = (:s, :s)),
        BasisState(2, "S", 1, 0; label = "2 n nbar", flavors = (:q, :q)),
        BasisState(2, "S", 1, 0; label = "2 s sbar", flavors = (:s, :s)),
    ] : BasisState[state for state in basis]
    length(basis_states) == length(diag) || throw(ArgumentError(
        "calibrated annihilation basis size does not match diagonal",
    ))
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        basis_states,
        matrix;
        mechanism = "rank_one_calibrated_pseudoscalar_annihilation",
        source = "GI Fig. 5 / Table III P1 diagnostic",
        notes = "Calibrated to digitized GI isoscalar pseudoscalar masses; separate PaperP1Annihilation and PaperP2Annihilation methods implement Eq. (18a,b).",
    )
    return MixingResult(block, fact.values, vectors)
end

"""
    isoscalar_general_annihilation_solution(params, basis; amplitude_A, L, multiplicity, J, smearing)

General non-pseudoscalar Eq. (16) annihilation block in the channel
`^{2S+1}L_J`. The matrix element from `q_i q̄_i → q_j q̄_j` is

`4π(2L+1) A(^{2S+1}L_J) [α_s(M_j²)α_s(M_i²)/π²]^{n/2} S_L(Ψ_j)S_L(Ψ_i)/(m_i m_j)`

with `n = 2` for `C = +` and `n = 3` for `C = −` (`C = (−1)^{L+S}`) and the
Eq. (17) smearing `S_L`. Each input retains its solver-native radial wave;
`momentum_wave` dispatches on that representation. Select an HO solver in the
calling spectrum calculation when reproducing the paper algorithm.
"""
function isoscalar_general_annihilation_solution(
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput};
    amplitude_A::Real,
    L::Integer,
    multiplicity::Integer,
    J::Integer,
    smearing::MomentumIntegralSmearing = MomentumIntegralSmearing(),
)
    n = length(basis)
    expected_L = L_LABELS[Int(L)]
    all(
        state.basis.L_label == expected_L &&
        state.basis.multiplicity == multiplicity &&
        state.basis.J == J for state in basis
    ) || throw(ArgumentError(
        "general annihilation basis does not match requested ^$(multiplicity)$(expected_L)_$(J) channel",
    ))
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
    vectors = _phase_fix_columns!(Matrix(fact.vectors))
    basis_states = _basis_states_from_inputs(basis)
    channel = @sprintf("^%d%s_%d", multiplicity, L_LABELS[L], J)
    block = MixingBlock(
        "isoscalar $channel annihilation",
        basis_states,
        matrix;
        mechanism = "general_eq16_annihilation",
        source = "GI Eq. (16) with Table III A($channel)=$(amplitude_A), n=$(n_gluons) gluons",
        notes = "Uses solver-native radial waves through the representation-neutral Eq. (17) momentum interface.",
    )
    return MixingResult(block, fact.values, vectors)
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
    smearing::MomentumIntegralSmearing = MomentumIntegralSmearing(),
)
    return isoscalar_general_annihilation_solution(
        params,
        basis;
        amplitude_A = params.annihilation.s1_A,
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
    vectors = _phase_fix_columns!(Matrix(fact.vectors))
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        _basis_states_from_inputs(basis),
        matrix;
        mechanism = "paper_p1_pseudoscalar_annihilation",
        source = "GI Eq. (16) with Eq. (18a) and Table III P1 constants",
        notes = "Uses the cached solver-native radial waves through the Eq. (17) momentum interface.",
    )
    return MixingResult(block, fact.values, vectors)
end

function isoscalar_pseudoscalar_annihilation_solution(
    model::PaperP2Annihilation,
    params::GIParameters,
    basis::AbstractVector{PseudoscalarAnnihilationBasisInput};
    maxiter::Integer = 80,
    tol::Real = 1.0e-10,
)
    n = length(basis)
    maxiter >= 1 || throw(ArgumentError("PaperP2Annihilation: maxiter must be positive"))
    tol > 0 || throw(ArgumentError("PaperP2Annihilation: tol must be positive"))
    starts = sort([state.diagonal_GeV for state in basis])
    masses = zeros(Float64, n)
    vectors = zeros(Float64, n, n)
    matrices = Vector{Matrix{Float64}}(undef, n)
    for level in 1:n
        M = starts[level]
        matrix = _paper_annihilation_matrix(model, params, basis; pole_mass_GeV = M)
        fact = eigen(Symmetric(matrix))
        converged = false
        for _ in 1:maxiter
            newM = fact.values[level]
            if abs(newM - M) <= tol
                converged = true
                break
            end
            M = 0.5 * (M + newM)
            matrix = _paper_annihilation_matrix(model, params, basis; pole_mass_GeV = M)
            fact = eigen(Symmetric(matrix))
        end
        converged || throw(ErrorException(
            "PaperP2Annihilation pole $level did not converge in $maxiter iterations",
        ))
        masses[level] = fact.values[level]
        vectors[:, level] .= fact.vectors[:, level]
        matrices[level] = matrix
    end
    _phase_fix_columns!(vectors)
    block = MixingBlock(
        "isoscalar ^1S_0 annihilation",
        _basis_states_from_inputs(basis),
        matrices[1];
        mechanism = "paper_p2_pseudoscalar_annihilation",
        source = "GI Eq. (16) with Eq. (18b) and Table III P2 constants",
        notes = "Each pole is solved as a fixed point of the mass-dependent Eq. (18b) matrix using the Eq. (17) S-wave momentum integral, so eigenvectors are not expected to be orthogonal.",
    )
    return MixingResult(block, masses, vectors, matrices)
end
