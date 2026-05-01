module GIModel

using LinearAlgebra
using Printf
using TOML
using CSV
using KrylovKit: eigsolve
using SpecialFunctions: erf

export GIParameters,
    ReferenceState,
    load_parameters,
    load_reference_spectrum,
    solve_sector,
    compare_sector,
    write_residual_report,
    parse_quark_masses,
    reduced_mass,
    channel_solution,
    fine_structure_split,
    LdotS,
    tensor_triplet_LJ,
    spin_dot,
    central_potential_mode,
    central_potential_values,
    CentralPotentialPath,
    central_potential_path

const ALPHA_COEFFS = (0.25, 0.15, 0.20)
const ALPHA_GAMMAS = (0.5, sqrt(10.0) / 2, sqrt(1000.0) / 2)
const L_SYMBOLS = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)

struct GIParameters
    masses::Dict{String, Float64}
    b::Float64
    c::Float64
    sigma0::Float64
    smearing_s::Float64
    appendix_a_smearing::Bool
    appendix_a_derivative_g::Bool
    appendix_a_closed_form::Bool
    appendix_a_momentum_sandwich::Bool
    contact_momentum_sandwich::Bool
    epsilon_c::Float64
    epsilon_t::Float64
    epsilon_so_vector::Float64
    epsilon_so_scalar::Float64
    fine_structure_momentum_sandwich::Bool
    fine_structure_smeared_kernels::Bool
    fine_structure::Bool
    k_spin_orbit::Float64
    k_tensor::Float64
    coulomb_1d_smear::Bool
end

struct ReferenceState
    sector::String
    quark_content::String
    composition::String
    n::Int
    multiplicity::Int
    L::String
    J::Int
    mass_GeV::Float64
    confidence::String
end

function load_parameters(path::AbstractString)
    raw = TOML.parsefile(path)
    masses = Dict(
        "u" => raw["masses"]["m_ud_avg_MeV"] / 1000,
        "d" => raw["masses"]["m_ud_avg_MeV"] / 1000,
        "q" => raw["masses"]["m_ud_avg_MeV"] / 1000,
        "s" => raw["masses"]["m_s_MeV"] / 1000,
        "c" => raw["masses"]["m_c_MeV"] / 1000,
        "b" => raw["masses"]["m_b_MeV"] / 1000,
    )
    rf = get(raw, "relativistic_factors", nothing)
    eps_c = isnothing(rf) ? 0.0 : get(rf, "epsilon_c", 0.0)
    eps_t = isnothing(rf) ? 0.0 : get(rf, "epsilon_t", 0.0)
    eps_v = isnothing(rf) ? 0.0 : get(rf, "epsilon_so_vector", 0.0)
    eps_s = isnothing(rf) ? 0.0 : get(rf, "epsilon_so_scalar", 0.0)
    fs = get(raw, "fine_structure", nothing)
    fine_on = isnothing(fs) ? true : get(fs, "enabled", true)
    k_so = isnothing(fs) ? 0.5 : get(fs, "k_spin_orbit", 0.5)
    k_tn = isnothing(fs) ? 0.4 : get(fs, "k_tensor", 0.4)
    GIParameters(
        masses,
        raw["potential"]["b_GeV2"],
        raw["potential"]["c_MeV"] / 1000,
        raw["relativistic_smearing"]["sigma0_GeV"],
        raw["relativistic_smearing"]["s"],
        get(raw["potential"], "appendix_a_smearing", false),
        get(raw["potential"], "appendix_a_derivative_g", false),
        get(raw["potential"], "appendix_a_closed_form", false),
        get(raw["potential"], "appendix_a_momentum_sandwich", false),
        get(raw["relativistic_factors"], "contact_momentum_sandwich", false),
        float(eps_c),
        float(eps_t),
        float(eps_v),
        float(eps_s),
        get(raw["relativistic_factors"], "fine_structure_momentum_sandwich", false),
        get(raw["relativistic_factors"], "fine_structure_smeared_kernels", false),
        fine_on,
        float(k_so),
        float(k_tn),
        get(raw["potential"], "coulomb_1d_smear", false),
    )
end

function load_reference_spectrum(path::AbstractString)
    states = ReferenceState[]
    for row in CSV.File(path)
        push!(
            states,
            ReferenceState(
                String(row.sector),
                String(row.quark_content),
                String(row.composition_raw),
                Int(row.n),
                Int(row.multiplicity),
                String(row.L),
                Int(row.J),
                Float64(row.mass_GeV),
                String(row.confidence),
            ),
        )
    end
    states
end

# Godfrey-Isgur parameterizes α_s(r) as a sum of error functions; use the
# library erf while keeping these thin wrappers to centralize derivative formulas.
function gi_erf(x::Real)
    erf(float(x))
end

function gi_erf_prime(x::Real)
    z = float(x)
    2 / sqrt(π) * exp(-z^2)
end

function gi_erf_second(x::Real)
    z = float(x)
    -4z / sqrt(π) * exp(-z^2)
end

alpha_s_r(r::Real) = sum(a * gi_erf(g * r) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

function central_potential(r::Real, params::GIParameters)
    params.b * r - (4 / 3) * alpha_s_r(r) / r + params.c
end

"""Coulomb piece G(r) = -4 α_s / (3 r) from the text; independent of the linear + constant term S(r) = b r + c."""
function static_coulomb_G(r::Real, params::GIParameters)
    ri = max(float(r), 1.0e-12)
    -(4 / 3) * alpha_s_r(ri) / ri
end

function static_confinement_S(r::Real, params::GIParameters)
    params.b * float(r) + params.c
end

#
# 3D isotropic Gaussian smearing of a spherically symmetric radial function V (|r|),
# Appendix A, Eqs. (A7)–(A8), PDF p. 36. Same σ as contact_smearing_sigma (A9), Table II.
# R > 0: one-dimensional form from the angle-integrated convolution (e.g. difference of
# Gaussians). R → 0: direct radial integral with isotropic 3D Gaussian.
#
function smear_3d_radial(v::AbstractVector{<:Real}, r::AbstractVector{<:Real}, σ::Real)
    n = length(r)
    n == 0 && return eltype(r)[]
    n == 1 && return v
    h = r[2] - r[1]
    w = fill(h, n)
    w[1] = h / 2
    w[n] = h / 2
    σf = max(float(σ), 1.0e-12)
    R0 = 0.25 * h
    out = similar(r, Float64)
    for i in eachindex(r)
        R = r[i]
        s = 0.0
        if R < R0
            for j in eachindex(r)
                rp = r[j]
                # ρ(r) = σ^3 / π^(3/2) exp(-σ^2 r^2), with σ in GeV and r in GeV^-1.
                ρ = σf^3 / (π^(3 / 2)) * exp(-(σf * rp)^2)
                s += 4 * π * rp^2 * w[j] * ρ * v[j]
            end
        else
            for j in eachindex(r)
                rp = r[j]
                # Angle-integrated convolution of a 3D isotropic Gaussian:
                # f̃(R) = (σ / (√π R)) ∫ dr' r' [e^{-σ^2 (R-r')^2} - e^{-σ^2 (R+r')^2}] f(r')
                pre = σf / (sqrt(π) * R)
                s += w[j] * pre * rp * (exp(-(σf * (R - rp))^2) - exp(-(σf * (R + rp))^2)) * v[j]
            end
        end
        out[i] = s
    end
    return out
end

# Experimental. Not the GI (A12)–(A13) smeared potential used in the paper: (A7)–(A8)
# with the Table II width applied to the pointwise G and S (Eqs. (11)–(13) orient.)
# is numerically uncontrolled on a fixed radial line when $\sigma$ is O(1): the 3D
# convolution weights the large-$r$ region by volume and can remove the $1/r$ well.
# Enable only for research; production defaults keep `appendix_a_smearing = false`.
function smeared_central_values(params::GIParameters, m1::Real, m2::Real, r::AbstractVector{<:Real})
    σ = contact_smearing_sigma(params, m1, m2)
    n = length(r)
    n < 2 && return [central_potential(ri, params) for ri in r]
    h = r[2] - r[1]
    rmax0 = r[end]
    # Kernel tail: exp(-(σ Δr)^2) at Δr = 8/σ gives exp(-64), effectively zero.
    n_tail = σ > 0 ? max(0, Int(ceil(8 / (σ * h)))) : 0
    r_ext = n_tail > 0 ? vcat(r, collect(range(rmax0 + h, rmax0 + n_tail * h; step = h))) : r
    g0 = [static_coulomb_G(ri, params) for ri in r_ext]
    s0 = [static_confinement_S(ri, params) for ri in r_ext]
    vsum = smear_3d_radial(g0, r_ext, σ) .+ smear_3d_radial(s0, r_ext, σ)
    return vsum[1:n]
end

function smeared_coulomb_G_closed(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = float(r)
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    τs = map(γ -> 1 / sqrt(1 / σ^2 + 1 / γ^2), ALPHA_GAMMAS)
    if abs(ri) < 1.0e-8
        return -sum(8 * α * τ / (3 * sqrt(π)) for (α, τ) in zip(ALPHA_COEFFS, τs))
    end
    -sum(4 * α * gi_erf(τ * ri) / (3 * ri) for (α, τ) in zip(ALPHA_COEFFS, τs))
end

function smeared_confinement_S_closed(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = float(r)
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    if abs(ri) < 1.0e-8
        return 2 * params.b / (sqrt(π) * σ) + params.c
    end
    z = σ * ri
    bracket = exp(-z^2) / (sqrt(π) * z) + (1 + 1 / (2 * z^2)) * gi_erf(z)
    params.b * ri * bracket + params.c
end

function appendix_a_closed_central_values(params::GIParameters, m1::Real, m2::Real, r::AbstractVector{<:Real})
    [
        smeared_coulomb_G_closed(params, m1, m2, ri) +
        smeared_confinement_S_closed(params, m1, m2, ri) for ri in r
    ]
end

include("radial_1d_coulomb_smear.jl")
include("appendix_a_derivative_potential.jl")

function reduced_mass(m1::Real, m2::Real)
    m1 * m2 / (m1 + m2)
end

function radial_grid(ngrid::Integer, rmax::Real)
    h = rmax / (ngrid + 1)
    collect(h:h:(ngrid * h)), h
end

function p2_operator(m::Real, L::Integer, r::AbstractVector, h::Real)
    ngrid = length(r)
    diagonal = similar(r)
    offdiag = fill(-1 / h^2, ngrid - 1)
    for i in eachindex(r)
        diagonal[i] = 2 / h^2 + L * (L + 1) / r[i]^2
    end
    SymTridiagonal(diagonal, offdiag)
end

function central_potential_mode(params::GIParameters)::Symbol
    if params.appendix_a_momentum_sandwich
        return :appendix_a_momentum_sandwich
    elseif params.appendix_a_closed_form
        return :appendix_a_closed_form
    elseif params.appendix_a_derivative_g
        return :appendix_a_derivative_g
    elseif params.appendix_a_smearing
        return :appendix_a_3d_a7a8
    elseif params.coulomb_1d_smear
        return :coulomb_1d
    end
    return :pointwise
end

function central_potential_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector;
    mode::Symbol = central_potential_mode(params),
)
    if mode == :pointwise
        return [central_potential(ri, params) for ri in r]
    elseif mode == :appendix_a_3d_a7a8
        return smeared_central_values(params, m1, m2, r)
    elseif mode == :coulomb_1d
        return coulomb_1d_smeared_central_values(params, m1, m2, r)
    elseif mode == :appendix_a_derivative_g
        return appendix_a_derivative_central_values(params, m1, m2, r)
    elseif mode == :appendix_a_closed_form
        return appendix_a_closed_central_values(params, m1, m2, r)
    elseif mode == :appendix_a_momentum_sandwich
        return appendix_a_closed_central_values(params, m1, m2, r)
    else
        error("unknown central potential mode: $mode")
    end
end

function potential_diagonal(params::GIParameters, m1::Real, m2::Real, r::AbstractVector)
    return central_potential_values(params, m1, m2, r)
end

function nonrelativistic_hamiltonian(params::GIParameters, m1::Real, m2::Real, L::Integer; ngrid::Integer = 900, rmax::Real = 24.0)
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

function lowest_eigenpairs(hamiltonian::AbstractMatrix, nlevels::Integer; eigensolver::Symbol = :full)
    if eigensolver == :full
        fact = eigen(hamiltonian)
        return fact.values[1:nlevels], fact.vectors[:, 1:nlevels]
    elseif eigensolver == :krylov
        values, vectors, info = eigsolve(hamiltonian, nlevels, :SR; issymmetric = true)
        length(values) >= nlevels ||
            error("Krylov eigensolver converged only $(length(values)) values for nlevels=$nlevels: $info")
        order = sortperm(real.(values))[1:nlevels]
        return real.(values[order]), hcat(vectors[order]...)
    else
        error("unknown eigensolver: $eigensolver")
    end
end

function appendix_a_momentum_sandwich_matrix(params::GIParameters, m1::Real, m2::Real, r::AbstractVector, p2_fact)
    λ = max.(p2_fact.values, 0)
    e1 = sqrt.(λ .+ m1^2)
    e2 = sqrt.(λ .+ m2^2)
    a_diag = sqrt.(1 .+ λ ./ (e1 .* e2))
    A = p2_fact.vectors * Diagonal(a_diag) * p2_fact.vectors'
    gdiag = Diagonal([smeared_coulomb_G_closed(params, m1, m2, ri) for ri in r])
    sdiag = Diagonal([smeared_confinement_S_closed(params, m1, m2, ri) for ri in r])
    Symmetric(A * gdiag * A + sdiag)
end

function relativistic_hamiltonian(params::GIParameters, m1::Real, m2::Real, L::Integer; ngrid::Integer = 450, rmax::Real = 24.0)
    r, h = radial_grid(ngrid, rmax)
    p2 = p2_operator(m1, L, r, h)
    p2_fact = eigen(p2)
    kinetic = sqrt_kinetic_matrix_from_eigen(p2_fact, m1) + sqrt_kinetic_matrix_from_eigen(p2_fact, m2)
    potential =
        if central_potential_mode(params) == :appendix_a_momentum_sandwich
            appendix_a_momentum_sandwich_matrix(params, m1, m2, r, p2_fact)
        else
            Diagonal(potential_diagonal(params, m1, m2, r))
        end
    Symmetric(kinetic + potential), r
end

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

function contact_smearing_sigma(params::GIParameters, m1::Real, m2::Real)
    # Appendix A (A9), PDF p. 36–37: universal σ(m1,m2) built from Table II σ0 and s.
    # We keep the paper's symmetric mass combinations explicit:
    #   mass_factor   = 4 m1 m2 / (m1 + m2)^2
    #   reduced_twice = 2 m1 m2 / (m1 + m2) = 2 μ
    # so σ^2 = σ0^2 * (1/2 + 1/2 * mass_factor^4) + s^2 * reduced_twice^2.
    mass_factor = 4 * m1 * m2 / (m1 + m2)^2
    reduced_twice = 2 * m1 * m2 / (m1 + m2)
    sqrt(params.sigma0^2 * (0.5 + 0.5 * mass_factor^4) + params.smearing_s^2 * reduced_twice^2)
end

function spin_dot(multiplicity::Integer)
    S = (multiplicity - 1) / 2
    0.5 * (S * (S + 1) - 1.5)
end

"""
3D normalized Gaussian regulator for a contact delta, with σ in GeV and r in GeV⁻¹.

This returns δ_σ(r) such that ∫ d³r δ_σ(r) = 1, i.e.
  4π ∫₀^∞ r² δ_σ(r) dr = 1.
"""
function delta_sigma_3d(r::Real, σ::Real)
    σ = float(σ)
    ri = float(r)
    σ > 0 || return 0.0
    return σ^3 / (π^(3 / 2)) * exp(-(σ * ri)^2)
end

function momentum_relativization_matrix(m1::Real, m2::Real, exponent::Real, p2_fact)
    λ = max.(p2_fact.values, 0)
    e1 = sqrt.(λ .+ m1^2)
    e2 = sqrt.(λ .+ m2^2)
    diag = (m1 * m2 ./ (e1 .* e2)) .^ exponent
    p2_fact.vectors * Diagonal(diag) * p2_fact.vectors'
end

function euclidean_expectation(vector::AbstractVector, operator::AbstractMatrix)
    v = collect(float.(vector))
    norm2 = sum(abs2, v)
    norm2 <= 0.0 && return 0.0
    dot(v, operator * v) / norm2
end

"""
First-order smeared contact hyperfine shift for S-waves.

Convention: the solver eigenvector is treated as the reduced radial wavefunction
`u(r)` on a uniform mesh with physical normalization `∫|u|² dr = 1`. For an
S-wave, `ψ(r) = u(r) / r · Y₀₀` and a 3D-normalized regulator `δ_σ(r)` satisfies
`∫ d³r δ_σ(r) = 1`. Therefore

`⟨α_s(r) δ_σ(r)⟩ = ∫ |u(r)|² α_s(r) δ_σ(r) dr`

with no extra `4π` factor.
"""
function contact_hyperfine_shift(params::GIParameters, m1::Real, m2::Real, L::String, multiplicity::Integer, vector::AbstractVector, r::AbstractVector)
    L == "S" || return 0.0
    multiplicity in (1, 3) || return 0.0
    length(r) >= 2 || return 0.0
    sigma = contact_smearing_sigma(params, m1, m2)
    h = r[2] - r[1]
    expectation = radial_expect_udr(
        vector,
        r,
        h,
        (ri, i) -> begin
            alpha_s_r(ri) * delta_sigma_3d(ri, sigma)
        end,
    )
    (1.0 + params.epsilon_c) * (32 * π / (9 * m1 * m2)) * expectation *
    spin_dot(multiplicity)
end

function contact_hyperfine_shift_momentum_sandwich(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::String,
    multiplicity::Integer,
    vector::AbstractVector,
    r::AbstractVector,
)
    L == "S" || return 0.0
    multiplicity in (1, 3) || return 0.0
    length(r) >= 2 || return 0.0
    h = r[2] - r[1]
    p2_fact = eigen(p2_operator(m1, 0, r, h))
    side_exponent = 0.25 + 0.5 * params.epsilon_c
    B = momentum_relativization_matrix(m1, m2, side_exponent, p2_fact)
    sigma = contact_smearing_sigma(params, m1, m2)
    kernel = Diagonal([alpha_s_r(ri) * delta_sigma_3d(ri, sigma) for ri in r])
    expectation = euclidean_expectation(vector, Symmetric(B * kernel * B))
    (32 * π / (9 * m1 * m2)) * expectation * spin_dot(multiplicity)
end

function contact_hyperfine_shift_active(params::GIParameters, args...)
    if params.contact_momentum_sandwich
        return contact_hyperfine_shift_momentum_sandwich(params, args...)
    end
    return contact_hyperfine_shift(params, args...)
end

include("masses_from_content.jl")
include("spin_fine_structure.jl")
include("appendix_a_status.jl")

function solve_sector(
    params::GIParameters,
    flavor::String;
    maxn::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    m = params.masses[flavor]
    results = Dict{Tuple{Int, String}, Float64}()
    for (symbol, L) in L_SYMBOLS
        levels = solve_channel(
            params, m, m, L;
            nlevels = maxn, ngrid = ngrid, rmax = rmax, kinetic = kinetic, eigensolver = eigensolver,
        )
        for n in 1:length(levels)
            results[(n, symbol)] = levels[n]
        end
    end
    results
end

function compare_sector(
    params::GIParameters,
    reference::Vector{ReferenceState},
    flavor::String;
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
    contact_hyperfine::Bool = true,
    use_fine_structure::Bool = true,
)
    m_fallback = params.masses[flavor]
    function masses_for(s::ReferenceState)
        try
            return parse_quark_masses(params, String(s.sector), String(s.quark_content))
        catch
            return m_fallback, m_fallback
        end
    end
    channel_cache = Dict{Tuple{Float64, Float64, String}, Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}}}()
    for state in reference
        m1, m2 = masses_for(state)
        cache_key = (round(m1, sigdigits = 12), round(m2, sigdigits = 12), state.L)
        if !haskey(channel_cache, cache_key)
            Lval = L_SYMBOLS[state.L]
            channel_cache[cache_key] = channel_solution(
                params, m1, m2, Lval;
                nlevels = 6, ngrid = ngrid, rmax = rmax, kinetic = kinetic, eigensolver = eigensolver,
            )
        end
    end
    rows = NamedTuple[]
    for state in reference
        m1, m2 = masses_for(state)
        cache_key = (round(m1, sigdigits = 12), round(m2, sigdigits = 12), state.L)
        haskey(channel_cache, cache_key) || continue
        values, vectors, r = channel_cache[cache_key]
        state.n <= length(values) || continue
        h = r[2] - r[1]
        central = values[state.n]
        contact_shift = 0.0
        spin_orbit_vector_shift = 0.0
        spin_orbit_thomas_shift = 0.0
        spin_orbit_shift = 0.0
        tensor_shift = 0.0
        fine_structure_shift = 0.0
        fine_structure_mass_convention = "disabled"
        if contact_hyperfine
            contact_shift = contact_hyperfine_shift_active(params, m1, m2, state.L, state.multiplicity, vectors[:, state.n], r)
        end
        if use_fine_structure && params.fine_structure
            comp = fine_structure_components(
                params, m1, m2, state.L, state.multiplicity, state.J,
                collect(vectors[:, state.n]), collect(r), h;
                enabled = true, k_spin_orbit = params.k_spin_orbit, k_tensor = params.k_tensor,
            )
            spin_orbit_vector_shift = comp.spin_orbit_vector
            spin_orbit_thomas_shift = comp.spin_orbit_thomas
            spin_orbit_shift = comp.spin_orbit
            tensor_shift = comp.tensor
            fine_structure_shift = comp.total
            fine_structure_mass_convention = isapprox(m1, m2; rtol = 0.0, atol = 0.0) ?
                                             "equal_mass" :
                                             "unequal_mass_equal_share_LdotS"
        end
        predicted = central + contact_shift + fine_structure_shift
        push!(
            rows,
            (
                sector = state.sector,
                state = state.composition,
                L = state.L,
                n = state.n,
                J = state.J,
                multiplicity = state.multiplicity,
                m1_GeV = m1,
                m2_GeV = m2,
                fine_structure_mass_convention = fine_structure_mass_convention,
                reference_GeV = state.mass_GeV,
                central_GeV = central,
                contact_shift_GeV = contact_shift,
                spin_orbit_vector_shift_GeV = spin_orbit_vector_shift,
                spin_orbit_thomas_shift_GeV = spin_orbit_thomas_shift,
                spin_orbit_shift_GeV = spin_orbit_shift,
                tensor_shift_GeV = tensor_shift,
                fine_structure_shift_GeV = fine_structure_shift,
                predicted_GeV = predicted,
                residual_MeV = 1000 * (predicted - state.mass_GeV),
                confidence = state.confidence,
            ),
        )
    end
    rows
end

function write_residual_report(
    path::AbstractString,
    title::AbstractString,
    rows;
    kinetic::Symbol = :relativistic,
    contact_hyperfine::Bool = true,
    appendix_a_smearing::Bool = false,
    appendix_a_derivative_g::Bool = false,
    appendix_a_closed_form::Bool = false,
    appendix_a_momentum_sandwich::Bool = false,
    contact_momentum_sandwich::Bool = false,
    fine_structure_momentum_sandwich::Bool = false,
    fine_structure_smeared_kernels::Bool = false,
    coulomb_1d_smear::Bool = false,
    use_fine_structure::Bool = true,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io)
        hyperfine_note =
            if contact_hyperfine && contact_momentum_sandwich
                "with GI momentum-sandwiched smeared S-wave contact hyperfine"
            elseif contact_hyperfine
                "with smeared S-wave contact hyperfine"
            else
                "without S-wave contact hyperfine"
            end
        fs_note =
            if use_fine_structure && fine_structure_momentum_sandwich && fine_structure_smeared_kernels
                " first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches and smeared-G/S derivative kernels; "
            elseif use_fine_structure && fine_structure_momentum_sandwich
                " first-order L·S (vector+Thomas) and OGE-tensor with GI momentum-factor sandwiches; "
            elseif use_fine_structure && fine_structure_smeared_kernels
                " first-order L·S (vector+Thomas) and OGE-tensor with smeared-G/S derivative kernels; "
            elseif use_fine_structure
                " first-order L·S (vector+Thomas) and OGE-tensor; "
            else
                " no first-order L·S/tensor; "
            end
        central_note = if appendix_a_momentum_sandwich
            "closed-form GI G̃(r), S̃(r), plus central Coulomb momentum sandwich `A(p)G̃A(p)` on the FD p² eigenbasis, "
        elseif appendix_a_closed_form
            "closed-form Gaussian-smeared GI G̃(r) and S̃(r), without the central Coulomb momentum sandwich, "
        elseif appendix_a_derivative_g
            "Appendix-A derivative proxy for G(r), `G + ∇²G/(4σ²)`, with pointwise S(r); not the full audited (A12)–(A13) expansion, "
        elseif appendix_a_smearing
            "experimental (A7)–(A8)-style 3D isotropic smearing of pointwise Coulomb G and confinement S (Table II σ₀, s), "
        elseif coulomb_1d_smear
            "1D Gaussian renormalization of G(r) only (pointwise S); same σ as contact (A9); not the full (A12)–(A13) expansion, "
        else
            "pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), "
        end
        println(
            io,
            "Model: finite-difference + `$kinetic` kinetic, $hyperfine_note,$fs_note",
            "GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`;",
            " ",
            central_note,
            "see `src/GIModel/`.",
        )
        println(io)
        println(io, "| state | reference GeV | baseline GeV | residual MeV | confidence |")
        println(io, "|---|---:|---:|---:|---|")
        for row in rows
            label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
            println(
                io,
                @sprintf(
                    "| `%s` | %.3f | %.3f | %+7.1f | %s |",
                    label,
                    row.reference_GeV,
                    row.predicted_GeV,
                    row.residual_MeV,
                    row.confidence,
                ),
            )
        end

        if !isempty(rows) && hasproperty(rows[1], :central_GeV)
            println(io)
            println(io, "## Contribution Breakdown (diagnostic)")
            println(io)
            println(io, "All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).")
            println(io)
            println(io, "| state | central GeV | contact MeV | L·S(vec) MeV | L·S(Thomas) MeV | L·S total MeV | tensor MeV | total shift MeV | predicted GeV |")
            println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for row in rows
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                total_shift = row.contact_shift_GeV + row.fine_structure_shift_GeV
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.3f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %.3f |",
                        label,
                        row.central_GeV,
                        1000 * row.contact_shift_GeV,
                        1000 * row.spin_orbit_vector_shift_GeV,
                        1000 * row.spin_orbit_thomas_shift_GeV,
                        1000 * row.spin_orbit_shift_GeV,
                        1000 * row.tensor_shift_GeV,
                        1000 * total_shift,
                        row.predicted_GeV,
                    ),
                )
            end
        end

        if use_fine_structure &&
           !isempty(rows) &&
           hasproperty(rows[1], :fine_structure_mass_convention) &&
           any(!isapprox(row.m1_GeV, row.m2_GeV; rtol = 0.0, atol = 0.0) for row in rows)
            println(io)
            println(io, "## Fine-Structure Mass Convention (audit note)")
            println(io)
            println(io, "Fine structure is currently implemented in terms of total `L·S` and a symmetric mass prefactor; this is exact for equal-mass `q\\bar q` but only a diagnostic convention for unequal masses (antisymmetric spin–orbit and mixing are not yet implemented).")
            println(io)
            println(io, "| state | m1 GeV | m2 GeV | convention |")
            println(io, "|---|---:|---:|---|")
            for row in rows
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                println(io, @sprintf("| `%s` | %.6f | %.6f | `%s` |", label, row.m1_GeV, row.m2_GeV, row.fine_structure_mass_convention))
            end
        end

        residuals = [abs(row.residual_MeV) for row in rows]
        if !isempty(residuals)
            println(io)
            println(io, @sprintf("Mean absolute residual: %.1f MeV.", sum(residuals) / length(residuals)))
            println(io, @sprintf("Max absolute residual: %.1f MeV.", maximum(residuals)))
        end
        groups = Dict{Tuple{Int, String}, Vector{eltype(rows)}}()
        for row in rows
            push!(get!(groups, (row.n, row.L), eltype(rows)[]), row)
        end
        println(io)
        println(io, "## Spin-Averaged Diagnostics")
        println(io)
        println(io, "Weighted by `2J+1` within each available `(n, L)` group.")
        println(io)
        println(io, "| multiplet | states | reference GeV | baseline GeV | residual MeV |")
        println(io, "|---|---:|---:|---:|---:|")
        for key in sort(collect(keys(groups)); by = x -> (x[2], x[1]))
            group = groups[key]
            weights = [2 * row.J + 1 for row in group]
            weight_sum = sum(weights)
            ref = sum(w * row.reference_GeV for (w, row) in zip(weights, group)) / weight_sum
            pred = sum(w * row.predicted_GeV for (w, row) in zip(weights, group)) / weight_sum
            println(io, @sprintf("| `%d%s` | %d | %.3f | %.3f | %+7.1f |", key[1], key[2], length(group), ref, pred, 1000 * (pred - ref)))
        end
        println(io)
        if appendix_a_momentum_sandwich
            println(io, "The spin-independent central path is the current GI reproduction candidate. Remaining heavy-quarkonium residuals should be read mainly as spin-dependent/operator-ordering and extraction-audit targets.")
        else
            println(io, "This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.")
        end
    end
end

end
