# First-order color-magnetic + Thomas-precession spin–orbit, plus OGE-tensor
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
#
# Public API (exported from GIModel.jl):
#   fine_structure_components, fine_structure_split, spin_orbit_mixing_components,
#   same_j_mixing, LdotS, tensor_triplet_LJ

function alpha_s_prime_r(r::Real)
    r = float(r)
    s = 0.0
    for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
        s += a * g * erf_prime(g * r)
    end
    return s
end

function alpha_s_second_r(r::Real)
    r = float(r)
    s = 0.0
    for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
        s += a * g^2 * erf_second(g * r)
    end
    return s
end

function coulomb_G_running(r::Real)
    ri = max(float(r), 1.0e-9)
    -(4.0 / 3.0) * alpha_s_r(ri) / ri
end

function coulomb_G_prime_running(r::Real)
    # d/dr [-(4/3) α_s(r) / r] = (4/3) [ α_s(r)/r² - α_s'(r)/r ].
    ri = max(float(r), 1.0e-9)
    α = alpha_s_r(ri)
    αp = alpha_s_prime_r(ri)
    (4.0 / 3.0) * (α / ri^2 - αp / ri)
end

function coulomb_G_second_running(r::Real)
    # d²/dr² [-(4/3) α_s(r) / r] = (4/3) [ 2 α_s'(r)/r² - 2 α_s(r)/r³ - α_s''(r)/r ].
    ri = max(float(r), 1.0e-9)
    α = alpha_s_r(ri)
    αp = alpha_s_prime_r(ri)
    αpp = alpha_s_second_r(ri)
    (4.0 / 3.0) * (2.0 * αp / ri^2 - 2.0 * α / ri^3 - αpp / ri)
end

function tensor_kernel_coulomb_running(r::Real)
    # Tensor kernel built from the Coulomb piece G(r) = -4 α_s(r) / (3 r):
    #   K(r) = (1/r) dG/dr - d²G/dr².
    ri = max(float(r), 1.0e-9)
    (1.0 / ri) * coulomb_G_prime_running(ri) - coulomb_G_second_running(ri)
end

smeared_coulomb_G_prime_closed(params::GIParameters, m::ConstituentMasses, r::Real) =
    smeared_coulomb_G_prime_closed(params, m.m1_GeV, m.m2_GeV, r)

function smeared_coulomb_G_prime_closed(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = max(float(r), 1.0e-7)
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    s = 0.0
    for (α, γ) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
        τ = 1 / sqrt(1 / σ^2 + 1 / γ^2)
        e = erf(τ * ri)
        ep = 2 * τ / sqrt(π) * exp(-(τ * ri)^2)
        s += -(4 * α / 3) * (ep / ri - e / ri^2)
    end
    return s
end

smeared_coulomb_G_second_closed(params::GIParameters, m::ConstituentMasses, r::Real) =
    smeared_coulomb_G_second_closed(params, m.m1_GeV, m.m2_GeV, r)

function smeared_coulomb_G_second_closed(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = max(float(r), 1.0e-7)
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    s = 0.0
    for (α, γ) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
        τ = 1 / sqrt(1 / σ^2 + 1 / γ^2)
        e = erf(τ * ri)
        ep = 2 * τ / sqrt(π) * exp(-(τ * ri)^2)
        fpp = -2 * τ^2 * ep - 2 * ep / ri^2 + 2 * e / ri^3
        s += -(4 * α / 3) * fpp
    end
    return s
end

tensor_kernel_smeared_coulomb(params::GIParameters, m::ConstituentMasses, r::Real) =
    tensor_kernel_smeared_coulomb(params, m.m1_GeV, m.m2_GeV, r)

function tensor_kernel_smeared_coulomb(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = max(float(r), 1.0e-7)
    (1.0 / ri) * smeared_coulomb_G_prime_closed(params, m1, m2, ri) -
    smeared_coulomb_G_second_closed(params, m1, m2, ri)
end

smeared_confinement_S_prime_closed(params::GIParameters, m::ConstituentMasses, r::Real) =
    smeared_confinement_S_prime_closed(params, m.m1_GeV, m.m2_GeV, r)

function smeared_confinement_S_prime_closed(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::Real,
)
    ri = float(r)
    abs(ri) < 1.0e-7 && return 0.0
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    z = σ * ri
    expz = exp(-z^2)
    h = ri + 1 / (2 * σ^2 * ri)
    hp = 1 - 1 / (2 * σ^2 * ri^2)
    params.b *
    ((-2 * σ * ri / sqrt(π)) * expz + hp * erf(z) + h * (2 * σ / sqrt(π)) * expz)
end

function dV_coul_central_dr(r::Real, params::GIParameters)
    # Central derivative for V_G(r) = G(r) = -4 α_s(r) / (3 r).
    # Kept as a thin wrapper to reduce divergence risk between spin–orbit and tensor kernels.
    ri = max(r, 1.0e-9)
    return coulomb_G_prime_running(ri)
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
    length(r) == length(u) ||
        throw(ArgumentError("physical_u_norm: length(r) != length(u)"))
    isfinite(float(h)) && h > 0 ||
        throw(ArgumentError("physical_u_norm: invalid mesh spacing h=$h"))
    # Convention guardrail: our expectation-value machinery assumes the solver eigenvector is sampled
    # on a *uniform* r mesh and that callers pass the correct mesh spacing `h`. A mismatch silently
    # rescales ⟨f(r)⟩ integrals and is a common source of normalization/convention drift.
    if length(r) >= 2
        hf = float(h)
        rprev = float(r[1])
        isfinite(rprev) || throw(ArgumentError("physical_u_norm: non-finite r[1]=$(r[1])"))
        for i = 2:length(r)
            ri = float(r[i])
            isfinite(ri) ||
                throw(ArgumentError("physical_u_norm: non-finite r[$i]=$(r[i])"))
            Δ = ri - rprev
            (Δ > 0 && isapprox(Δ, hf; rtol = 1e-10, atol = 1e-12)) || throw(
                ArgumentError(
                    "physical_u_norm: non-uniform r mesh or h mismatch at i=$i (Δr=$Δ, h=$hf)",
                ),
            )
            rprev = ri
        end
    end
    s = 0.0
    for i in eachindex(r)
        s += abs2(float(u[i])) * h
    end
    s <= 0.0 && return 0.0
    return 1.0 / sqrt(s)
end

function radial_expect_udr(
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real,
    f::F,
) where {F<:Function}
    length(r) == length(u) ||
        throw(ArgumentError("radial_expect_udr: length(r) != length(u)"))
    n = physical_u_norm(r, h, u)
    n == 0.0 && return 0.0
    s = 0.0
    for i in eachindex(r)
        ui = n * float(u[i])
        s += abs2(ui) * h * f(float(r[i]), i)
    end
    return s
end

function radial_expect_momentum_sandwich(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real,
    epsilon::Real,
    f::F,
) where {F<:Function}
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    length(r) == length(u) ||
        throw(ArgumentError("radial_expect_momentum_sandwich: length(r) != length(u)"))
    length(r) >= 2 || return 0.0
    # Reuse the same uniform-mesh convention guard as the diagonal expectation path.
    physical_u_norm(r, h, u)
    p2_fact = eigen(p2_operator(m1, L, r, h))
    side_exponent = 0.25 + 0.5 * epsilon
    B = momentum_relativization_matrix(m1, m2, side_exponent, p2_fact)
    kernel = Diagonal([f(float(ri), i) for (i, ri) in enumerate(r)])
    return euclidean_expectation(u, Symmetric(B * kernel * B))
end

function fine_structure_components(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    radial::RadialWaveOnUniformMesh;
    enabled::Bool = true,
    k_spin_orbit::Real = 1.0,
    k_tensor::Real = 1.0,
)
    Ls = multiplet.L_label
    multiplicity = multiplet.multiplicity
    J = multiplet.J
    u = radial.u
    r = radial.r
    h = radial.h
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    !enabled && return (
        I_cm = 0.0,
        I_tp = 0.0,
        I_tk = 0.0,
        spin_orbit_vector = 0.0,
        spin_orbit_thomas = 0.0,
        spin_orbit = 0.0,
        tensor = 0.0,
        total = 0.0,
    )
    Ln = L_SYMBOLS[Ls]
    S = (multiplicity - 1) ÷ 2
    (Ln == 0 || S < 0) && return (
        I_cm = 0.0,
        I_tp = 0.0,
        I_tk = 0.0,
        spin_orbit_vector = 0.0,
        spin_orbit_thomas = 0.0,
        spin_orbit = 0.0,
        tensor = 0.0,
        total = 0.0,
    )
    S == 1 || return (
        I_cm = 0.0,
        I_tp = 0.0,
        I_tk = 0.0,
        spin_orbit_vector = 0.0,
        spin_orbit_thomas = 0.0,
        spin_orbit = 0.0,
        tensor = 0.0,
        total = 0.0,
    )

    m1, m2 = float(m1), float(m2)
    # Paper Eq. (6) uses α_s(r)/r^3 directly (no α_s' term). Eq. (7) Thomas-precession
    # term uses (1/2r) dH_conf/dr, which *does* include α_s'(r) via d/dr[-α_s(r)/r].
    # For unequal masses the exact operator splits into symmetric and antisymmetric
    # spin–orbit pieces. We currently keep only the symmetric L·S contraction:
    #   L·(S_i/m_i^2 + S_j/m_j^2)  →  (1/2)(1/m1^2 + 1/m2^2) L·S
    #   L·[(1/mi+1/mj)(S_i/mi + S_j/mj)]
    #     → (1/2)(1/m1^2 + 1/m2^2 + 2/(m1 m2)) L·S
    inv2_tp = 0.5 * (1.0 / m1^2 + 1.0 / m2^2)
    inv2_cm = 0.5 * (1.0 / m1^2 + 1.0 / m2^2 + 2.0 / (m1 * m2))
    expect_kernel(epsilon, f) =
        params.fine_structure_momentum_sandwich ?
        radial_expect_momentum_sandwich(params, masses, Ln, u, r, h, epsilon, f) :
        radial_expect_udr(u, r, h, f)
    Icm = expect_kernel(
        params.epsilon_so_vector,
        (ri, i) -> begin
            r0 = max(ri, 1.0e-8)
            if params.fine_structure_smeared_kernels
                (1.0 / r0) * smeared_coulomb_G_prime_closed(params, masses, r0)
            else
                (4.0 / 3.0) * alpha_s_r(r0) / r0^3
            end
        end,
    )
    Itp = expect_kernel(
        params.epsilon_so_scalar,
        (ri, i) -> begin
            r0 = max(ri, 1.0e-8)
            if params.fine_structure_smeared_kernels
                (1.0 / (2.0 * r0)) * (
                    smeared_confinement_S_prime_closed(params, masses, r0) +
                    smeared_coulomb_G_prime_closed(params, masses, r0)
                )
            else
                (1.0 / (2.0 * r0)) * (params.b + dV_coul_central_dr(r0, params))
            end
        end,
    )
    # The active research path uses the smeared Coulomb tensor kernel from
    # derivatives of G~(r). The legacy branch keeps the pointwise running-Coulomb
    # kernel, including the α_s'(r) and α_s''(r) pieces.
    Itk = expect_kernel(
        params.epsilon_t,
        (ri, i) ->
            params.fine_structure_smeared_kernels ?
            tensor_kernel_smeared_coulomb(params, masses, ri) :
            tensor_kernel_coulomb_running(ri),
    )
    ls = LdotS(Ln, 1, J)
    vec_term =
        params.fine_structure_momentum_sandwich ? Icm :
        (1.0 + params.epsilon_so_vector) * Icm
    thomas_term =
        params.fine_structure_momentum_sandwich ? Itp :
        (1.0 + params.epsilon_so_scalar) * Itp
    spin_orbit_vector = k_spin_orbit * inv2_cm * ls * vec_term
    spin_orbit_thomas = k_spin_orbit * (-inv2_tp) * ls * thomas_term
    spin_orbit = spin_orbit_vector + spin_orbit_thomas
    # Coulomb-limit check: for G(r) = -4 α_s / (3 r) with constant α_s,
    #   (1/r dG/dr - d²G/dr²) = 4 α_s / r³
    # so the tensor prefactor reduces to 4/(3 m1 m2) times ⟨α_s / r³⟩.
    tensor_scale = params.fine_structure_momentum_sandwich ? 1.0 : (1.0 + params.epsilon_t)
    tensor =
        tensor_scale *
        k_tensor *
        (1.0 / (3.0 * m1 * m2)) *
        Itk *
        tensor_triplet_LJ(Ln, J, 1)
    return (
        I_cm = Icm,
        I_tp = Itp,
        I_tk = Itk,
        spin_orbit_vector = spin_orbit_vector,
        spin_orbit_thomas = spin_orbit_thomas,
        spin_orbit = spin_orbit,
        tensor = tensor,
        total = spin_orbit + tensor,
    )
end

function spin_orbit_radial_integrals(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real,
)
    expect_kernel(epsilon, f) =
        params.fine_structure_momentum_sandwich ?
        radial_expect_momentum_sandwich(params, masses, L, u, r, h, epsilon, f) :
        radial_expect_udr(u, r, h, f)
    Icm = expect_kernel(
        params.epsilon_so_vector,
        (ri, i) -> begin
            r0 = max(ri, 1.0e-8)
            if params.fine_structure_smeared_kernels
                (1.0 / r0) * smeared_coulomb_G_prime_closed(params, masses, r0)
            else
                (4.0 / 3.0) * alpha_s_r(r0) / r0^3
            end
        end,
    )
    Itp = expect_kernel(
        params.epsilon_so_scalar,
        (ri, i) -> begin
            r0 = max(ri, 1.0e-8)
            if params.fine_structure_smeared_kernels
                (1.0 / (2.0 * r0)) * (
                    smeared_confinement_S_prime_closed(params, masses, r0) +
                    smeared_coulomb_G_prime_closed(params, masses, r0)
                )
            else
                (1.0 / (2.0 * r0)) * (params.b + dV_coul_central_dr(r0, params))
            end
        end,
    )
    vec_term =
        params.fine_structure_momentum_sandwich ? Icm :
        (1.0 + params.epsilon_so_vector) * Icm
    thomas_term =
        params.fine_structure_momentum_sandwich ? Itp :
        (1.0 + params.epsilon_so_scalar) * Itp
    return (I_cm = Icm, I_tp = Itp, vec_term = vec_term, thomas_term = thomas_term)
end

"""
    spin_orbit_mixing_components(params, masses, L_label, radial; ...)

Diagnostic antisymmetric spin-orbit matrix element for the same-`J` basis
`|n ^1L_L>` / `|n ^3L_L>`. The angular convention is
`<^1L_L| L·(S1-S2) |^3L_L> = sqrt(L(L+1))`; the sign of the reported mixing
angle is therefore tied to the basis ordering used in [`same_j_mixing`](@ref).
"""
function spin_orbit_mixing_components(
    params::GIParameters,
    masses::ConstituentMasses,
    L_label::AbstractString,
    radial::RadialWaveOnUniformMesh;
    enabled::Bool = true,
    k_spin_orbit::Real = 1.0,
)
    L = L_SYMBOLS[String(L_label)]
    if !enabled || L == 0
        return (
            I_cm = 0.0,
            I_tp = 0.0,
            vector = 0.0,
            thomas = 0.0,
            total = 0.0,
            angular = 0.0,
        )
    end
    m1 = float(masses.m1_GeV)
    m2 = float(masses.m2_GeV)
    radial_terms = spin_orbit_radial_integrals(params, masses, L, radial.u, radial.r, radial.h)
    angular = sqrt(L * (L + 1.0))
    inv2_asym = 0.5 * (1.0 / m1^2 - 1.0 / m2^2)
    vector = k_spin_orbit * inv2_asym * angular * radial_terms.vec_term
    thomas = k_spin_orbit * (-inv2_asym) * angular * radial_terms.thomas_term
    return (
        I_cm = radial_terms.I_cm,
        I_tp = radial_terms.I_tp,
        vector = vector,
        thomas = thomas,
        total = vector + thomas,
        angular = angular,
    )
end

"""
    same_j_mixing(singlet_mass, triplet_mass, offdiag)

Diagonalize the `(^1L_L, ^3L_L)` same-`J` mass matrix. The returned angle
uses the Fig. 9 convention
`low = cos(theta) * singlet + sin(theta) * triplet`.
"""
function same_j_mixing(singlet_mass::Real, triplet_mass::Real, offdiag::Real)
    matrix = Symmetric([float(singlet_mass) float(offdiag); float(offdiag) float(triplet_mass)])
    fact = eigen(matrix)
    low_vec = fact.vectors[:, 1]
    if low_vec[1] < 0
        low_vec = -low_vec
    end
    theta = atan(low_vec[2], low_vec[1])
    return (
        matrix = Matrix(matrix),
        masses = fact.values,
        vectors = fact.vectors,
        theta_rad = theta,
        theta_deg = theta * 180 / π,
    )
end

function fine_structure_split(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    radial::RadialWaveOnUniformMesh;
    enabled::Bool = true,
    k_spin_orbit::Real = 1.0,
    k_tensor::Real = 1.0,
)
    return fine_structure_components(
        params,
        masses,
        multiplet,
        radial;
        enabled = enabled,
        k_spin_orbit = k_spin_orbit,
        k_tensor = k_tensor,
    ).total
end

function fine_structure_components(
    params::GIParameters,
    masses::ConstituentMasses,
    Ls::AbstractString,
    multiplicity::Integer,
    J::Integer,
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real;
    enabled::Bool = true,
    k_spin_orbit::Real = 1.0,
    k_tensor::Real = 1.0,
)
    return fine_structure_components(
        params,
        masses,
        FineStructureMultiplet(Ls, multiplicity, J),
        RadialWaveOnUniformMesh(u, r, h);
        enabled = enabled,
        k_spin_orbit = k_spin_orbit,
        k_tensor = k_tensor,
    )
end

"""Convenience: same as [`fine_structure_components`](@ref)`(params, ConstituentMasses(m1, m2), ...)`."""
function fine_structure_components(
    params::GIParameters,
    m1::Real,
    m2::Real,
    Ls::AbstractString,
    multiplicity::Integer,
    J::Integer,
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real;
    kwargs...,
)
    return fine_structure_components(params, ConstituentMasses(m1, m2), Ls, multiplicity, J, u, r, h; kwargs...)
end

function fine_structure_split(
    params::GIParameters,
    masses::ConstituentMasses,
    Ls::AbstractString,
    multiplicity::Integer,
    J::Integer,
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real;
    kwargs...,
)
    return fine_structure_split(
        params,
        masses,
        FineStructureMultiplet(Ls, multiplicity, J),
        RadialWaveOnUniformMesh(u, r, h);
        kwargs...,
    )
end

"""Convenience: same as [`fine_structure_split`](@ref)`(params, ConstituentMasses(m1, m2), ...)`."""
function fine_structure_split(
    params::GIParameters,
    m1::Real,
    m2::Real,
    Ls::AbstractString,
    multiplicity::Integer,
    J::Integer,
    u::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real;
    kwargs...,
)
    return fine_structure_split(
        params,
        ConstituentMasses(m1, m2),
        Ls,
        multiplicity,
        J,
        u,
        r,
        h;
        kwargs...,
    )
end
