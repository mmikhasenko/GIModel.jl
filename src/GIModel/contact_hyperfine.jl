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
