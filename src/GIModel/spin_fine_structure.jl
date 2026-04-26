# First-order color-magnetic + Thomas (scalar confinement) spin–orbit, plus OGE-tensor
# in the structure of the paper, Eqs. (3)–(7) (text), using Table II ε in (A10).
# Radial integrals: physical normalization ∫ 4π u² dr = 1 for the reduced Schrödinger
# radial u(r) on the same mesh as the FD solver.
#
# A global scale k_spin_orbit / k_tensor bridges the present FD + semirelativistic path to
# the large HO-basis results in the original paper; defaults are in parameters.toml.

function alpha_s_prime_r(r::Real)
    r = float(r)
    s = 0.0
    for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
        z = g * r
        e = exp(-z^2)
        s += a * 2.0 * g / (sqrt(π) * 1) * e
    end
    return s
end

function dV_coul_central_dr(r::Real, params::GIParameters)
    ri = max(r, 1.0e-9)
    α = alpha_s_r(ri)
    αp = alpha_s_prime_r(ri)
    return 4.0 * αp / (3.0 * ri) - 4.0 * α / (3.0 * ri^2)
end

function tensor_f_LJ_barnes(L::Int, J::Int, S::Int)
    S == 1 || return 0.0
    (L, J) == (1, 0) && return -4.0
    (L, J) == (1, 1) && return 1.0
    (L, J) == (1, 2) && return -0.2
    (L, J) == (2, 1) && return -1.0
    (L, J) == (2, 2) && return 0.5
    (L, J) == (2, 3) && return -0.1
    (L, J) == (3, 2) && return 0.4
    (L, J) == (3, 3) && return -0.2
    (L, J) == (3, 4) && return 0.1
    (L, J) == (4, 3) && return 0.2
    (L, J) == (4, 4) && return -0.1
    (L, J) == (4, 5) && return 0.05
    return 0.0
end

function LdotS(L::Int, S::Int, J::Int)
    0.5 * (J * (J + 1) - L * (L + 1) - S * (S + 1))
end

function physical_u_norm(r::Vector{Float64}, h::Real, u::Vector{Float64})
    s = 0.0
    for i in eachindex(r)
        s += 4.0 * π * abs2(u[i]) * h
    end
    s <= 0.0 && return 0.0
    return 1.0 / sqrt(s)
end

function smeared_r_inv(params::GIParameters, m1::Real, m2::Real, r::Real, p::Int)
    σ = contact_smearing_sigma(params, m1, m2)
    rs2 = r^2 + σ^2
    if p == 1
        return 1.0 / max(r, 1.0e-8)
    elseif p == 2
        return 1.0 / rs2
    elseif p == 3
        return 1.0 / (rs2 * sqrt(max(rs2, 1.0e-20)))
    end
    return 0.0
end

function radial_expect_4piudr(
    u::Vector{Float64},
    r::Vector{Float64},
    h::Real,
    f::F,
) where {F<:Function}
    n = physical_u_norm(r, h, u)
    n == 0.0 && return 0.0
    s = 0.0
    for i in eachindex(r)
        ui = n * u[i]
        s += 4.0 * π * abs2(ui) * h * f(r[i], i)
    end
    return s
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
    !enabled && return 0.0
    Ln = L_SYMBOLS[Ls]
    S = (multiplicity - 1) ÷ 2
    (Ln == 0 || S < 0) && return 0.0
    if S == 0
        return 0.0
    end
    m1, m2 = float(m1), float(m2)
    inv2 = 0.25 * (1.0 / m1^2 + 1.0 / m2^2)
    Ivp = radial_expect_4piudr(u, r, h, (ri, i) -> (1.0 / max(ri, 1.0e-8)) * dV_coul_central_dr(ri, params))
    I1 = radial_expect_4piudr(u, r, h, (ri, i) -> 1.0 / max(ri, 1.0e-8))
    Its = radial_expect_4piudr(
        u, r, h,
        (ri, i) -> alpha_s_r(ri) * smeared_r_inv(params, m1, m2, ri, 3),
    )
    ls = LdotS(Ln, 1, J)
    vec_term = (1.0 + params.epsilon_so_vector) * Ivp
    thomas_term = (1.0 + params.epsilon_so_scalar) * params.b * I1
    delta_so = k_spin_orbit * inv2 * ls * (vec_term + thomas_term)
    tq = (1.0 + params.epsilon_t) * k_tensor * (1.0 / (3.0 * m1 * m2)) * Its * tensor_f_LJ_barnes(Ln, J, 1)
    return delta_so + tq
end
