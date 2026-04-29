# First-order color-magnetic + Thomas (scalar confinement) spin–orbit, plus OGE-tensor
# in the structure of the paper, Eqs. (3)–(7) (text), using Table II ε in (A10).
# Radial integrals: we treat the FD eigenvector as the reduced Schrödinger radial
# wavefunction u(r) sampled on a uniform mesh. The physical normalization is
#   ∫ |u(r)|² dr = 1
# (no extra 4π factor; the spherical-harmonic angular integral is already unity for
# normalized Y_{LM}). On a uniform mesh with spacing h, the discrete proxy is
#   ∑ |uᵢ|² h = 1.
#
# A global scale k_spin_orbit / k_tensor bridges the present FD + semirelativistic path to
# the large HO-basis results in the original paper; defaults are in parameters.toml.

function alpha_s_prime_r(r::Real)
    r = float(r)
    s = 0.0
    for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
        s += a * g * erf_approx_prime(g * r)
    end
    return s
end

function dV_coul_central_dr(r::Real, params::GIParameters)
    ri = max(r, 1.0e-9)
    α = alpha_s_r(ri)
    αp = alpha_s_prime_r(ri)
    return 4.0 * α / (3.0 * ri^2) - 4.0 * αp / (3.0 * ri)
end

function tensor_triplet_LJ(L::Int, J::Int, S::Int)
    S == 1 || return 0.0
    L <= 0 && return 0.0
    J == L - 1 && return -2.0 * (L + 1) / (2L - 1)
    J == L && return 2.0
    J == L + 1 && return -2.0 * L / (2L + 3)
    return 0.0
end

function LdotS(L::Int, S::Int, J::Int)
    0.5 * (J * (J + 1) - L * (L + 1) - S * (S + 1))
end

function physical_u_norm(r::AbstractVector{<:Real}, h::Real, u::AbstractVector{<:Real})
    length(r) == length(u) || throw(ArgumentError("physical_u_norm: length(r) != length(u)"))
    isfinite(float(h)) && h > 0 || throw(ArgumentError("physical_u_norm: invalid mesh spacing h=$h"))
    s = 0.0
    for i in eachindex(r)
        s += abs2(float(u[i])) * h
    end
    s <= 0.0 && return 0.0
    return 1.0 / sqrt(s)
end

function smeared_r_inv(params::GIParameters, m1::Real, m2::Real, r::Real, p::Int)
    # Helper for Appendix-A-style "relativized" regulators of 1/r^p singularities.
    #
    # Convention audit: Table II σ(m1,m2) has units of GeV, while the FD mesh uses
    # r in GeV⁻¹, so the corresponding smear length scale in r-space is 1/σ.
    #
    # This helper is currently unused by the active diagnostics, but keeping it
    # unit-consistent avoids accidentally reintroducing a factor-of-σ^2 bug when
    # wiring it into future Appendix A operators.
    σ = contact_smearing_sigma(params, m1, m2)
    ℓ = 1.0 / max(float(σ), 1.0e-12) # smear length in GeV⁻¹
    rs2 = float(r)^2 + ℓ^2
    if p == 1
        return 1.0 / sqrt(max(rs2, 1.0e-20))
    elseif p == 2
        return 1.0 / max(rs2, 1.0e-20)
    elseif p == 3
        rs = sqrt(max(rs2, 1.0e-20))
        return 1.0 / (rs2 * rs)
    end
    return 0.0
end

function radial_expect_udr(
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real,
    f::F,
) where {F<:Function}
    length(r) == length(u) || throw(ArgumentError("radial_expect_udr: length(r) != length(u)"))
    n = physical_u_norm(r, h, u)
    n == 0.0 && return 0.0
    s = 0.0
    for i in eachindex(r)
        ui = n * float(u[i])
        s += abs2(ui) * h * f(float(r[i]), i)
    end
    return s
end

function fine_structure_components(
    params::GIParameters,
    m1::Real,
    m2::Real,
    Ls::String,
    multiplicity::Int,
    J::Int,
    u::Vector{Float64},
    r::Vector{Float64},
    h::Real;
    enabled::Bool = true,
    k_spin_orbit::Real = 1.0,
    k_tensor::Real = 1.0,
)
    !enabled && return (spin_orbit = 0.0, tensor = 0.0, total = 0.0)
    Ln = L_SYMBOLS[Ls]
    S = (multiplicity - 1) ÷ 2
    (Ln == 0 || S < 0) && return (spin_orbit = 0.0, tensor = 0.0, total = 0.0)
    S == 1 || return (spin_orbit = 0.0, tensor = 0.0, total = 0.0)

    m1, m2 = float(m1), float(m2)
    inv2 = 0.25 * (1.0 / m1^2 + 1.0 / m2^2)
    Ivp = radial_expect_udr(u, r, h, (ri, i) -> (1.0 / max(ri, 1.0e-8)) * dV_coul_central_dr(ri, params))
    I1 = radial_expect_udr(u, r, h, (ri, i) -> 1.0 / max(ri, 1.0e-8))
    # For L>0 the FD radial wave function suppresses the origin. Using the same
    # broad Gaussian width as the S-wave contact term over-damps tensor
    # splittings; the full GI tensor term should come from derivatives of the
    # smeared G(r), but this unsmeared alpha_s/r^3 proxy is a better diagnostic
    # until that Appendix A operator is implemented.
    Its = radial_expect_udr(u, r, h, (ri, i) -> alpha_s_r(ri) / max(ri, 1.0e-8)^3)
    ls = LdotS(Ln, 1, J)
    vec_term = (1.0 + params.epsilon_so_vector) * Ivp
    thomas_term = (1.0 + params.epsilon_so_scalar) * params.b * I1
    spin_orbit = k_spin_orbit * inv2 * ls * (3 * vec_term - thomas_term)
    tensor = (1.0 + params.epsilon_t) * k_tensor * (1.0 / (3.0 * m1 * m2)) * Its * tensor_triplet_LJ(Ln, J, 1)
    return (spin_orbit = spin_orbit, tensor = tensor, total = spin_orbit + tensor)
end

function fine_structure_split(
    params::GIParameters,
    m1::Real,
    m2::Real,
    Ls::String,
    multiplicity::Int,
    J::Int,
    u::Vector{Float64},
    r::Vector{Float64},
    h::Real;
    enabled::Bool = true,
    k_spin_orbit::Real = 1.0,
    k_tensor::Real = 1.0,
)
    fine_structure_components(
        params,
        m1,
        m2,
        Ls,
        multiplicity,
        J,
        u,
        r,
        h;
        enabled = enabled,
        k_spin_orbit = k_spin_orbit,
        k_tensor = k_tensor,
    ).total
end
