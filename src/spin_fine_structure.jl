# First-order color-magnetic + Thomas-precession spin–orbit, plus OGE-tensor
# in the structure of the paper, Eqs. (3)–(7) (text), using Table II ε in the
# post-A14 Appendix-A momentum-factor prescription.
# Physics-level radial integrals consume the native `RadialWave` interface.
# FD assembly uses reduced-wave samples with ∫|u(r)|²dr=1; HO assembly uses
# coefficients with the same physical normalization.
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
    params.potential.b *
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

function tensor_triplet_offdiag_sameJ(J::Int, S::Int)
    S == 1 || return 0.0
    J <= 0 && return 0.0
    return 6.0 * sqrt(J * (J + 1.0)) / (2J + 1)
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

function radial_cross_expect_udr(
    u_left::AbstractVector{<:Real},
    u_right::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real,
    f::F,
) where {F<:Function}
    length(r) == length(u_left) == length(u_right) ||
        throw(ArgumentError("radial_cross_expect_udr: vector lengths differ"))
    n_left = physical_u_norm(r, h, u_left)
    n_right = physical_u_norm(r, h, u_right)
    (n_left == 0.0 || n_right == 0.0) && return 0.0
    s = 0.0
    for i in eachindex(r)
        s += (n_left * float(u_left[i])) * (n_right * float(u_right[i])) * h *
             f(float(r[i]), i)
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
    p2_fact = eigen(p2_operator(params, m1, L, r, h))
    side_exponent = gi_spin_dependent_side_exponent(epsilon)
    B = momentum_relativization_matrix(m1, m2, side_exponent, p2_fact)
    kernel = Diagonal([f(float(ri), i) for (i, ri) in enumerate(r)])
    return euclidean_expectation(u, Symmetric(B * kernel * B))
end

radial_expect_momentum_sandwich(
    params::GIParameters, masses::ConstituentMasses, L::Integer,
    wave::MeshWave, epsilon::Real, f,
) = radial_expect_momentum_sandwich(
    params, masses, L, wave.u, wave.r, wave.h, epsilon, f,
)

function _ho_expansion_value(L::Integer, beta::Real, coefficients, r::Real)
    return sum(
        coefficients[n + 1] * ho_reduced_radial(n, L, beta, r) for
        n in 0:(length(coefficients)-1)
    )
end

function radial_expect_momentum_sandwich(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    wave::OscillatorWave,
    epsilon::Real,
    f,
)
    L == wave.L || throw(ArgumentError(
        "radial_expect_momentum_sandwich: L=$L does not match wave L=$(wave.L)",
    ))
    p2_fact = eigen(Symmetric(Matrix(ho_p2_matrix(L, wave.beta, length(wave.coefficients)))))
    B = momentum_relativization_matrix(
        masses.m1_GeV,
        masses.m2_GeV,
        gi_spin_dependent_side_exponent(epsilon),
        p2_fact,
    )
    transformed = B * wave.coefficients
    rmax = _oscillator_coordinate_cutoff(wave)
    value, _ = quadgk(
        r -> _ho_expansion_value(L, wave.beta, transformed, r)^2 * f(r, 0),
        0.0,
        rmax;
        rtol = 1e-9,
    )
    return value / wave_norm(wave)
end

function radial_cross_expect_momentum_sandwich(
    params::GIParameters,
    masses::ConstituentMasses,
    L_left::Integer,
    u_left::AbstractVector{<:Real},
    L_right::Integer,
    u_right::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    h::Real,
    epsilon::Real,
    f::F,
) where {F<:Function}
    length(r) == length(u_left) == length(u_right) ||
        throw(ArgumentError("radial_cross_expect_momentum_sandwich: vector lengths differ"))
    # Uniform-mesh guardrails (return discarded; see physical_u_norm).
    physical_u_norm(r, h, u_left)
    physical_u_norm(r, h, u_right)
    # Normalization-invariant cross expectation. Unlike the diagonal path,
    # which divides by ⟨u|u⟩ via `euclidean_expectation`, the bare cross product
    # dot(u_left, B K B u_right) scales with the norms of *both* inputs. The FD
    # solver returns Euclidean-normalized eigenvectors (‖u‖₂ = 1) while the HO
    # path returns physically-normalized reconstructions (‖u‖₂ = 1/√h), so an
    # unnormalized cross element inflated the HO off-diagonal mixing by ≈ 1/h.
    # Dividing by ‖u_left‖₂ ‖u_right‖₂ makes the element basis-independent and
    # matches the diagonal Rayleigh-quotient convention (the explicit `h` in the
    # non-sandwich `radial_cross_expect_udr` cancels to the same expression).
    nl = sqrt(dot(u_left, u_left))
    nr = sqrt(dot(u_right, u_right))
    (nl == 0.0 || nr == 0.0) && return 0.0
    p2_left = eigen(p2_operator(params, masses.m1_GeV, L_left, r, h))
    p2_right = eigen(p2_operator(params, masses.m1_GeV, L_right, r, h))
    side_exponent = gi_spin_dependent_side_exponent(epsilon)
    B_left = momentum_relativization_matrix(masses.m1_GeV, masses.m2_GeV, side_exponent, p2_left)
    B_right = momentum_relativization_matrix(masses.m1_GeV, masses.m2_GeV, side_exponent, p2_right)
    kernel = Diagonal([f(float(ri), i) for (i, ri) in enumerate(r)])
    return dot(u_left, B_left * kernel * B_right * u_right) / (nl * nr)
end

radial_cross_expect_momentum_sandwich(
    params::GIParameters, masses::ConstituentMasses,
    L_left::Integer, left::MeshWave, L_right::Integer, right::MeshWave,
    epsilon::Real, f,
) = begin
    length(left.r) == length(right.r) &&
        all(isapprox.(left.r, right.r; rtol = 1e-10, atol = 1e-12)) ||
        throw(ArgumentError("radial cross expectation requires a shared mesh"))
    radial_cross_expect_momentum_sandwich(
        params, masses, L_left, left.u, L_right, right.u,
        left.r, left.h, epsilon, f,
    )
end

function radial_cross_expect_momentum_sandwich(
    params::GIParameters,
    masses::ConstituentMasses,
    L_left::Integer,
    left::OscillatorWave,
    L_right::Integer,
    right::OscillatorWave,
    epsilon::Real,
    f,
)
    L_left == left.L && L_right == right.L || throw(ArgumentError(
        "radial_cross_expect_momentum_sandwich: orbital labels do not match waves",
    ))
    side_exponent = gi_spin_dependent_side_exponent(epsilon)
    left_p2 = eigen(Symmetric(Matrix(ho_p2_matrix(
        L_left, left.beta, length(left.coefficients),
    ))))
    right_p2 = eigen(Symmetric(Matrix(ho_p2_matrix(
        L_right, right.beta, length(right.coefficients),
    ))))
    B_left = momentum_relativization_matrix(
        masses.m1_GeV, masses.m2_GeV, side_exponent, left_p2,
    )
    B_right = momentum_relativization_matrix(
        masses.m1_GeV, masses.m2_GeV, side_exponent, right_p2,
    )
    c_left = B_left * left.coefficients
    c_right = B_right * right.coefficients
    rmax = max(
        _oscillator_coordinate_cutoff(left),
        _oscillator_coordinate_cutoff(right),
    )
    value, _ = quadgk(
        r -> _ho_expansion_value(L_left, left.beta, c_left, r) *
             _ho_expansion_value(L_right, right.beta, c_right, r) * f(r, 0),
        0.0,
        rmax;
        rtol = 1e-9,
    )
    return value / sqrt(wave_norm(left) * wave_norm(right))
end

_mass_pair(m1::Real, m2::Real) = ConstituentMasses(float(m1), float(m2))

function _vector_so_kernel(params, pair, r)
    r0 = max(float(r), 1e-8)
    return params.factors.fine_structure_smeared_kernels ?
           smeared_coulomb_G_prime_closed(params, pair, r0) / r0 :
           (4 / 3) * alpha_s_r(r0) / r0^3
end

function _scalar_so_kernel(params, pair, r)
    r0 = max(float(r), 1e-8)
    return params.factors.fine_structure_smeared_kernels ?
           smeared_confinement_S_prime_closed(params, pair, r0) / r0 :
           params.potential.b / r0
end

"""
    fine_structure_radial_kernels(params, masses, r)

Return the six local radial kernels appearing in Appendix A15-A16 before
their momentum-dependent factors are applied. The fields are `vector_11`,
`vector_22`, `vector_12`, `scalar_11`, `scalar_22`, and `tensor_12`.

These are the coordinate-space middle operators in the GI sandwiches
`B(p²) K(r) B(p²)`, not complete potentials: spin-angular coefficients and
the explicit mass denominators in A15-A16 are deliberately not folded in.
Broadcast this scalar method over a radial grid when plotting a profile.
"""
function fine_structure_radial_kernels(
    params::GIParameters,
    masses::ConstituentMasses,
    r::Real,
)
    m1, m2 = masses.m1_GeV, masses.m2_GeV
    pair11, pair22 = _mass_pair(m1, m1), _mass_pair(m2, m2)
    return (
        vector_11 = _vector_so_kernel(params, pair11, r),
        vector_22 = _vector_so_kernel(params, pair22, r),
        vector_12 = _vector_so_kernel(params, masses, r),
        scalar_11 = _scalar_so_kernel(params, pair11, r),
        scalar_22 = _scalar_so_kernel(params, pair22, r),
        tensor_12 = params.factors.fine_structure_smeared_kernels ?
                    tensor_kernel_smeared_coulomb(params, masses, r) :
                    tensor_kernel_coulomb_running(r),
    )
end

function _spin_expectation(params, pair, L, wave, epsilon, kernel)
    value = params.factors.fine_structure_momentum_sandwich ?
            radial_expect_momentum_sandwich(
                params, pair, L, wave, epsilon, (r, _) -> kernel(r),
            ) : radial_expect(wave, kernel)
    return params.factors.fine_structure_momentum_sandwich ? value : (1 + epsilon) * value
end

function _spin_cross_expectation(params, pair, Lleft, left, Lright, right, epsilon, kernel)
    value = params.factors.fine_structure_momentum_sandwich ?
            radial_cross_expect_momentum_sandwich(
                params, pair, Lleft, left, Lright, right, epsilon, (r, _) -> kernel(r),
            ) : radial_overlap(left, right, kernel)
    return params.factors.fine_structure_momentum_sandwich ? value : (1 + epsilon) * value
end

function spin_orbit_radial_integrals(params, masses, L, left::RadialWave, right::RadialWave)
    m1, m2 = masses.m1_GeV, masses.m2_GeV
    pair11, pair22 = _mass_pair(m1, m1), _mass_pair(m2, m2)
    cross(pair, eps, kernel) = _spin_cross_expectation(
        params, pair, L, left, L, right, eps, kernel,
    )
    return (
        vector_11 = cross(pair11, params.factors.epsilon_so_vector,
                          r -> _vector_so_kernel(params, pair11, r)),
        vector_22 = cross(pair22, params.factors.epsilon_so_vector,
                          r -> _vector_so_kernel(params, pair22, r)),
        scalar_11 = cross(pair11, params.factors.epsilon_so_scalar,
                          r -> _scalar_so_kernel(params, pair11, r)),
        scalar_22 = cross(pair22, params.factors.epsilon_so_scalar,
                          r -> _scalar_so_kernel(params, pair22, r)),
    )
end

spin_orbit_radial_integrals(params, masses, L, wave::RadialWave) =
    spin_orbit_radial_integrals(params, masses, L, wave, wave)

function _zero_fine_components()
    return (
        I_vector_11 = 0.0, I_vector_22 = 0.0, I_vector_12 = 0.0,
        I_scalar_11 = 0.0, I_scalar_22 = 0.0, I_tk = 0.0,
        spin_orbit_vector = 0.0, spin_orbit_thomas = 0.0,
        spin_orbit = 0.0, tensor = 0.0, total = 0.0,
    )
end

function _zero_fine_matrices(n::Integer)
    zero_matrix = Symmetric(zeros(Float64, n, n))
    return (
        spin_orbit_vector = zero_matrix,
        spin_orbit_thomas = zero_matrix,
        spin_orbit = zero_matrix,
        tensor = zero_matrix,
        total = zero_matrix,
    )
end

# Diagonal L·S and S12 vanish for L=0 or total spin S=0. This is
# angular algebra, unlike the smeared contact interaction, which acts in both.
diagonal_fine_structure_active(params, L, multiplicity) =
    params.fine_structure.enabled && L > 0 && multiplicity == 3

# One angular/mass definition for matrices and scalar expectations.
function _fine_structure_algebra(masses, L, J, G11, G22, G12, S11, S22, T12)
    m1, m2 = masses.m1_GeV, masses.m2_GeV
    ls = LdotS(Int(L), 1, Int(J))
    vector = ls * (G11 / (4m1^2) + G22 / (4m2^2) + G12 / (m1 * m2))
    thomas = -ls * (S11 / (4m1^2) + S22 / (4m2^2))
    tensor = tensor_triplet_LJ(Int(L), Int(J), 1) * T12 / (12m1 * m2)
    return (spin_orbit_vector = vector, spin_orbit_thomas = thomas,
            spin_orbit = vector + thomas, tensor = tensor,
            total = vector + thomas + tensor)
end

# A15-A16 operator algebra is independent of the radial representation.  FD
# and HO differ only in how the six radial sandwich matrices below are built;
# keeping their angular/mass assembly here prevents the two numerical paths
# from acquiring different physics through copy-and-paste edits.
function _assemble_fine_structure_matrices(
    masses::ConstituentMasses,
    L::Integer,
    J::Integer,
    G11::AbstractMatrix,
    G22::AbstractMatrix,
    G12::AbstractMatrix,
    S11::AbstractMatrix,
    S22::AbstractMatrix,
    T12::AbstractMatrix,
)
    matrices = (G11, G22, G12, S11, S22, T12)
    size0 = size(G11)
    size0[1] == size0[2] || throw(ArgumentError(
        "fine-structure radial matrices must be square; got $size0",
    ))
    all(size(matrix) == size0 for matrix in matrices) || throw(ArgumentError(
        "fine-structure radial matrices must have one common size",
    ))

    return map(Symmetric, _fine_structure_algebra(
        masses, L, J, G11, G22, G12, S11, S22, T12,
    ))
end

function fine_structure_components(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    radial::RadialWave;
    enabled::Bool = true,
)
    L = L_SYMBOLS[multiplet.L_label]
    (!enabled || !diagonal_fine_structure_active(params, L, multiplet.multiplicity)) &&
        return _zero_fine_components()

    so = spin_orbit_radial_integrals(params, masses, L, radial)
    I12 = _spin_expectation(
        params, masses, L, radial, params.factors.epsilon_so_vector,
        r -> _vector_so_kernel(params, masses, r),
    )
    Itk = _spin_expectation(
        params, masses, L, radial, params.factors.epsilon_t,
        r -> params.factors.fine_structure_smeared_kernels ?
             tensor_kernel_smeared_coulomb(params, masses, r) :
             tensor_kernel_coulomb_running(r),
    )
    contributions = _fine_structure_algebra(
        masses, L, multiplet.J, so.vector_11, so.vector_22, I12,
        so.scalar_11, so.scalar_22, Itk,
    )
    return merge((
        I_vector_11 = so.vector_11, I_vector_22 = so.vector_22,
        I_vector_12 = I12, I_scalar_11 = so.scalar_11,
        I_scalar_22 = so.scalar_22, I_tk = Itk,
    ), contributions)
end

"""
    fine_structure_grid_matrices(params, masses, J, r, h; L = 1, multiplicity = 3)

Triplet `³L_J` spin-orbit and tensor potentials as dense operators on the
uniform radial mesh `r`, term-by-term identical to Appendix A15-A16. The
returned fields match
[`ho_fine_structure_matrices`](@ref): `spin_orbit_vector`,
`spin_orbit_thomas`, `spin_orbit`, `tensor`, and `total`.

Requires the smeared momentum-sandwich prescription; throws otherwise.
"""
function fine_structure_grid_matrices(
    params::GIParameters,
    masses::ConstituentMasses,
    J::Integer,
    r::AbstractVector,
    h::Real;
    L::Integer = 1,
    multiplicity::Integer = 3,
)
    diagonal_fine_structure_active(params, L, multiplicity) ||
        return _zero_fine_matrices(length(r))
    params.factors.fine_structure_momentum_sandwich &&
        params.factors.fine_structure_smeared_kernels ||
        error("fine_structure_grid_operator requires the Appendix-A smeared momentum-sandwich path")
    m1 = masses.m1_GeV
    m2 = masses.m2_GeV
    n = length(r)
    p2_fact = eigen(p2_operator(params, m1, L, r, h))
    side(pair, eps) = momentum_relativization_matrix(
        pair.m1_GeV,
        pair.m2_GeV,
        gi_spin_dependent_side_exponent(eps),
        p2_fact,
    )
    pair11, pair22 = _mass_pair(m1, m1), _mass_pair(m2, m2)
    sandwich(pair, eps, kernel) = begin
        B = side(pair, eps)
        B * Diagonal([kernel(ri) for ri in r]) * B
    end
    G11 = sandwich(pair11, params.factors.epsilon_so_vector,
                   ri -> _vector_so_kernel(params, pair11, ri))
    G22 = sandwich(pair22, params.factors.epsilon_so_vector,
                   ri -> _vector_so_kernel(params, pair22, ri))
    G12 = sandwich(masses, params.factors.epsilon_so_vector,
                   ri -> _vector_so_kernel(params, masses, ri))
    S11 = sandwich(pair11, params.factors.epsilon_so_scalar,
                   ri -> _scalar_so_kernel(params, pair11, ri))
    S22 = sandwich(pair22, params.factors.epsilon_so_scalar,
                   ri -> _scalar_so_kernel(params, pair22, ri))
    T12 = sandwich(masses, params.factors.epsilon_t,
                   ri -> tensor_kernel_smeared_coulomb(params, masses, ri))
    return _assemble_fine_structure_matrices(
        masses, L, J, G11, G22, G12, S11, S22, T12,
    )
end

"""Summed FD spin-orbit plus tensor operator; see [`fine_structure_grid_matrices`](@ref)."""
fine_structure_grid_operator(args...; kwargs...) =
    fine_structure_grid_matrices(args...; kwargs...).total

"""
    ho_fine_structure_matrices(params, masses, L, multiplicity, J, beta, nbasis)

Native HO matrices for the vector spin-orbit, scalar/Thomas spin-orbit, and
tensor terms of one fixed `(L,S,J)` sector. Each radial kernel is enclosed by
its own exact spectral momentum factor; no sampled grid operator is used.
"""
function ho_fine_structure_matrices(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    multiplicity::Integer,
    J::Integer,
    beta::Real,
    nbasis::Integer,
)
    if !diagonal_fine_structure_active(params, L, multiplicity)
        return _zero_fine_matrices(nbasis)
    end
    params.factors.fine_structure_momentum_sandwich &&
        params.factors.fine_structure_smeared_kernels || throw(ArgumentError(
        "native HO fine structure requires smeared Appendix-A momentum sandwiches",
    ))

    m1, m2 = masses.m1_GeV, masses.m2_GeV
    pair11, pair22 = _mass_pair(m1, m1), _mass_pair(m2, m2)
    sandwich(pair, eps, kernel) = ho_momentum_sandwich_matrix(
        L, beta, nbasis, pair, eps, kernel,
    )
    G11 = sandwich(pair11, params.factors.epsilon_so_vector,
                   r -> _vector_so_kernel(params, pair11, r))
    G22 = sandwich(pair22, params.factors.epsilon_so_vector,
                   r -> _vector_so_kernel(params, pair22, r))
    G12 = sandwich(masses, params.factors.epsilon_so_vector,
                   r -> _vector_so_kernel(params, masses, r))
    S11 = sandwich(pair11, params.factors.epsilon_so_scalar,
                   r -> _scalar_so_kernel(params, pair11, r))
    S22 = sandwich(pair22, params.factors.epsilon_so_scalar,
                   r -> _scalar_so_kernel(params, pair22, r))
    T12 = ho_momentum_sandwich_matrix(
        L, beta, nbasis, masses, params.factors.epsilon_t,
        r -> tensor_kernel_smeared_coulomb(params, masses, r),
        # The beta scan deliberately visits diffuse endpoint bases even for
        # bottomonium. At that irrelevant endpoint the narrow tensor kernel
        # reaches 1e-8 matrix accuracy by nq=8192; its A15 prefactor suppresses
        # the residual far below the 0.1 MeV eigenvalue convergence target.
        rtol = 1e-8,
    )
    return _assemble_fine_structure_matrices(
        masses, L, J, G11, G22, G12, S11, S22, T12,
    )
end

"""
    spin_orbit_mixing_components(params, masses, L_label, radial; ...)

Antisymmetric spin-orbit matrix element for the same-`J` basis
`|n ^1L_L>` / `|n ^3L_L>`. The angular convention is
`<^1L_L| L·(S1-S2) |^3L_L> = sqrt(L(L+1))`; the sign of the reported mixing
angle is therefore tied to the basis ordering used in [`same_j_mixing`](@ref).
"""
function spin_orbit_mixing_components(
    params::GIParameters,
    masses::ConstituentMasses,
    L_label::AbstractString,
    radial::RadialWave;
    kwargs...,
)
    return spin_orbit_mixing_components(
        params, masses, L_label, radial, radial; kwargs...,
    )
end

function spin_orbit_mixing_components(
    params::GIParameters,
    masses::ConstituentMasses,
    L_label::AbstractString,
    radial_left::RadialWave,
    radial_right::RadialWave;
    enabled::Bool = true,
)
    L = L_SYMBOLS[String(L_label)]
    if !enabled || L == 0
        return (
            I_vector_11 = 0.0,
            I_vector_22 = 0.0,
            I_scalar_11 = 0.0,
            I_scalar_22 = 0.0,
            vector = 0.0,
            thomas = 0.0,
            total = 0.0,
            angular = 0.0,
        )
    end
    m1 = float(masses.m1_GeV)
    m2 = float(masses.m2_GeV)
    radial_terms = spin_orbit_radial_integrals(
        params, masses, L, radial_left, radial_right,
    )
    angular = sqrt(L * (L + 1.0))
    vector = angular * (radial_terms.vector_11 / (4m1^2) -
                        radial_terms.vector_22 / (4m2^2))
    thomas = angular * (-radial_terms.scalar_11 / (4m1^2) +
                        radial_terms.scalar_22 / (4m2^2))
    return (
        I_vector_11 = radial_terms.vector_11,
        I_vector_22 = radial_terms.vector_22,
        I_scalar_11 = radial_terms.scalar_11,
        I_scalar_22 = radial_terms.scalar_22,
        vector = vector,
        thomas = thomas,
        total = vector + thomas,
        angular = angular,
    )
end

function tensor_mixing_components(
    params::GIParameters,
    masses::ConstituentMasses,
    radial_left::RadialWave,
    radial_right::RadialWave,
    J::Integer;
    enabled::Bool = true,
)
    Jn = Int(J)
    L_left = Jn - 1
    L_right = Jn + 1
    if !enabled || Jn <= 0
        return (I_tk = 0.0, angular = 0.0, total = 0.0)
    end
    kernel =
        (ri, i) ->
            params.factors.fine_structure_smeared_kernels ?
            tensor_kernel_smeared_coulomb(params, masses, ri) :
            tensor_kernel_coulomb_running(ri)
    Itk =
        params.factors.fine_structure_momentum_sandwich ?
        radial_cross_expect_momentum_sandwich(
            params,
            masses,
            L_left,
            radial_left,
            L_right,
            radial_right,
            params.factors.epsilon_t,
            kernel,
        ) :
        radial_overlap(radial_left, radial_right, ri -> kernel(ri, 0))
    angular = tensor_triplet_offdiag_sameJ(Jn, 1)
    total = Itk * angular / (12masses.m1_GeV * masses.m2_GeV)
    return (I_tk = Itk, angular = angular, total = total)
end

"""
    same_j_mixing(singlet_mass, triplet_mass, offdiag)

Diagonalize the `(^1L_L, ^3L_L)` same-`J` mass matrix. The returned angle
uses the Fig. 9 convention
`low = cos(theta) * singlet + sin(theta) * triplet`.
"""
function same_j_mixing(
    singlet_mass::Real,
    triplet_mass::Real,
    offdiag::Real;
    basis::AbstractVector{BasisState} = [
        BasisState(1, "L", 1, 0; label = "^1L_L"),
        BasisState(1, "L", 3, 0; label = "^3L_L"),
    ],
)
    length(basis) == 2 || throw(ArgumentError("same_j_mixing requires two basis states"))
    block = MixingBlock(
        "same-J ^1L_L/^3L_L",
        basis,
        [float(singlet_mass) float(offdiag); float(offdiag) float(triplet_mass)];
        mechanism = "antisymmetric_spin_orbit",
    )
    result = diagonalize_mixing_block(block)
    low_vec = result.vectors[:, 1]
    if low_vec[1] < 0
        low_vec = -low_vec
    end
    theta = atan(low_vec[2], low_vec[1])
    return (
        result = result,
        block = block,
        matrix = block.matrix,
        masses = result.masses,
        vectors = result.vectors,
        theta_rad = theta,
        theta_deg = theta * 180 / π,
    )
end

function fine_structure_split(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    radial::RadialWave;
    enabled::Bool = true,
)
    return fine_structure_components(
        params,
        masses,
        multiplet,
        radial;
        enabled = enabled,
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
)
    return fine_structure_components(
        params,
        masses,
        FineStructureMultiplet(Ls, multiplicity, J),
        MeshWave(u, r, h);
        enabled = enabled,
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
        MeshWave(u, r, h);
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
