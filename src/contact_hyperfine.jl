# Public API (exported from GIModel.jl):
#   spin_dot

contact_smearing_sigma(params::GIParameters, m::ConstituentMasses) =
    contact_smearing_sigma(params, m.m1_GeV, m.m2_GeV)

function contact_smearing_sigma(params::GIParameters, m1::Real, m2::Real)
    # Appendix A (A9), PDF p. 36–37: universal σ(m1,m2) built from Table II σ0 and s.
    # We keep the paper's symmetric mass combinations explicit:
    #   mass_factor   = 4 m1 m2 / (m1 + m2)^2
    #   reduced_twice = 2 m1 m2 / (m1 + m2) = 2 μ
    # so σ^2 = σ0^2 * (1/2 + 1/2 * mass_factor^4) + s^2 * reduced_twice^2.
    mass_factor = 4 * m1 * m2 / (m1 + m2)^2
    reduced_twice = 2 * m1 * m2 / (m1 + m2)
    sqrt(
        params.sigma0^2 * (0.5 + 0.5 * mass_factor^4) +
        params.smearing_s^2 * reduced_twice^2,
    )
end

"""
    spin_dot(multiplicity)
    spin_dot(multiplet::FineStructureMultiplet)

`⟨S₁·S₂⟩` for a `q q̄` pair with total-spin multiplicity `2S+1`:
`-3/4` for the singlet (`multiplicity = 1`), `+1/4` for the triplet (`3`).
"""
function spin_dot(multiplicity::Integer)
    S = (multiplicity - 1) / 2
    0.5 * (S * (S + 1) - 1.5)
end

spin_dot(multiplet::FineStructureMultiplet) = spin_dot(multiplet.multiplicity)

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

"""
Appendix A's post-A14 prescription places
`(m1*m2/(E1*E2))^(1/2 + epsilon_i)` on each side of a spin-dependent potential.
With `epsilon_i = 0`, the two-sided product turns the
nonrelativistic `1/(m1*m2)` strength into `1/(E1*E2)`.
"""
gi_spin_dependent_side_exponent(epsilon::Real) = 0.5 + float(epsilon)

function euclidean_expectation(vector::AbstractVector, operator::AbstractMatrix)
    v = collect(float.(vector))
    norm2 = sum(abs2, v)
    norm2 <= 0.0 && return 0.0
    dot(v, operator * v) / norm2
end

function contact_hyperfine_operator(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    r::AbstractVector,
)
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    n = length(r)
    if L != "S" || !(multiplicity in (1, 3)) || n < 2
        return Symmetric(zeros(Float64, n, n))
    end
    h = r[2] - r[1]
    p2_fact = eigen(p2_operator(params, m1, 0, r, h))
    side_exponent = gi_spin_dependent_side_exponent(params.epsilon_c)
    B = momentum_relativization_matrix(m1, m2, side_exponent, p2_fact)
    sigma = contact_smearing_sigma(params, masses)
    kernel = Diagonal([alpha_s_r(ri) * delta_sigma_3d(ri, sigma) for ri in r])
    strength = (32 * π / (9 * m1 * m2)) * spin_dot(multiplicity)
    return Symmetric(strength * (B * kernel * B))
end

function _contact_hyperfine_shift_diagonal(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    vector::AbstractVector,
    r::AbstractVector,
)
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    L == "S" || return 0.0
    multiplicity in (1, 3) || return 0.0
    length(r) >= 2 || return 0.0
    sigma = contact_smearing_sigma(params, masses)
    h = r[2] - r[1]
    expectation = radial_expect_udr(
        vector,
        r,
        h,
        (ri, i) -> begin
            alpha_s_r(ri) * delta_sigma_3d(ri, sigma)
        end,
    )
    (1.0 + params.epsilon_c) *
    (32 * π / (9 * m1 * m2)) *
    expectation *
    spin_dot(multiplicity)
end

function _contact_hyperfine_shift_momentum_sandwich_diagonal(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    vector::AbstractVector,
    r::AbstractVector,
)
    L == "S" || return 0.0
    multiplicity in (1, 3) || return 0.0
    length(r) >= 2 || return 0.0
    operator = contact_hyperfine_operator(params, masses, L, multiplicity, r)
    euclidean_expectation(vector, operator)
end

function contact_hyperfine_nonperturbative_levels(
    params::GIParameters{FiniteDifferenceBasis},
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    r::AbstractVector,
    nlevels::Integer,
)
    if !params.contact_momentum_sandwich || L != "S" || !(multiplicity in (1, 3)) || length(r) < 2
        return Float64[]
    end
    h = r[2] - r[1]
    rmax = h * (length(r) + 1)
    hamiltonian, rebuilt_r =
        relativistic_hamiltonian(params, masses, 0; ngrid = length(r), rmax = rmax)
    length(rebuilt_r) == length(r) || error("rebuilt S-wave grid changed length")
    operator = contact_hyperfine_operator(params, masses, L, multiplicity, rebuilt_r)
    levels, _vectors =
        lowest_eigenpairs(Symmetric(Matrix(hamiltonian) + Matrix(operator)), nlevels)
    return levels
end

function contact_hyperfine_nonperturbative_levels(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    r::AbstractVector,
    nlevels::Integer,
)
    return Float64[]
end

"""
First-order smeared contact hyperfine shift for S-waves.

Convention: the solver eigenvector is treated as the reduced radial wavefunction
`u(r)` on a uniform mesh with physical normalization `∫|u|² dr = 1`. For an
S-wave, `ψ(r) = u(r) / r · Y₀₀` and a 3D-normalized regulator `δ_σ(r)` satisfies
`∫ d³r δ_σ(r) = 1`. Therefore

`⟨α_s(r) δ_σ(r)⟩ = ∫ |u(r)|² α_s(r) δ_σ(r) dr`

with no extra `4π` factor.

Uses only [`FineStructureMultiplet`](@ref).`L_label` and `.multiplicity`; `.J` is unused (same multiplet object as fine-structure).
"""
function contact_hyperfine_shift(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    wave::RadialWaveOnUniformMesh,
)
    return _contact_hyperfine_shift_diagonal(
        params,
        masses,
        multiplet.L_label,
        multiplet.multiplicity,
        wave.u,
        wave.r,
    )
end

"""Convenience: same as [`contact_hyperfine_shift`](@ref)`(params, ConstituentMasses(m1, m2), ...)`."""
function contact_hyperfine_shift(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::AbstractString,
    multiplicity::Integer,
    vector::AbstractVector,
    r::AbstractVector,
)
    return _contact_hyperfine_shift_diagonal(
        params,
        ConstituentMasses(m1, m2),
        L,
        multiplicity,
        vector,
        r,
    )
end

function contact_hyperfine_shift_momentum_sandwich(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    wave::RadialWaveOnUniformMesh,
)
    return _contact_hyperfine_shift_momentum_sandwich_diagonal(
        params,
        masses,
        multiplet.L_label,
        multiplet.multiplicity,
        wave.u,
        wave.r,
    )
end

function contact_hyperfine_shift_momentum_sandwich(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::AbstractString,
    multiplicity::Integer,
    vector::AbstractVector,
    r::AbstractVector,
)
    return _contact_hyperfine_shift_momentum_sandwich_diagonal(
        params,
        ConstituentMasses(m1, m2),
        L,
        multiplicity,
        vector,
        r,
    )
end

function contact_hyperfine_shift_active(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    wave::RadialWaveOnUniformMesh,
)
    if params.contact_momentum_sandwich
        return contact_hyperfine_shift_momentum_sandwich(params, masses, multiplet, wave)
    end
    return contact_hyperfine_shift(params, masses, multiplet, wave)
end

function contact_hyperfine_shift_active(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::AbstractString,
    multiplicity::Integer,
    vector::AbstractVector,
    r::AbstractVector,
)
    mm = ConstituentMasses(m1, m2)
    if params.contact_momentum_sandwich
        return _contact_hyperfine_shift_momentum_sandwich_diagonal(
            params,
            mm,
            L,
            multiplicity,
            vector,
            r,
        )
    end
    return _contact_hyperfine_shift_diagonal(params, mm, L, multiplicity, vector, r)
end
