# Public API (exported from GIModel.jl):
#   spin_dot, contact_smearing_sigma

raw"""
    contact_smearing_sigma(params, m1, m2) -> Float64
    contact_smearing_sigma(params, masses::ConstituentMasses) -> Float64

The Appendix A (A9) universal smearing width ``\\sigma(m_1, m_2)`` in GeV,

    σ² = σ₀² (1/2 + 1/2 [4 m₁m₂/(m₁+m₂)²]⁴) + s² (2 m₁m₂/(m₁+m₂))²

built from the Table II inputs `σ₀` and `s`. Together with the relativistic
weight ``(m_1 m_2 / E_1 E_2)^{1/2 + \\epsilon_i}`` this is **the entire route by
which quark mass enters the model** — the potential parameters (`b`, `c`) and
[`alpha_s_q`](@ref) never see a mass or a flavor. σ grows monotonically with the
constituent masses, which is why the smeared contact term (and with it the
hyperfine splitting) collapses toward heavy quarkonium.
"""
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
        params.smearing.sigma0^2 * (0.5 + 0.5 * mass_factor^4) +
        params.smearing.s^2 * reduced_twice^2,
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

raw"""
    smeared_contact_kernel(params, masses, r)

The radial contact kernel obtained from the Laplacian of the Appendix-A
smeared Coulomb potential.  With

``\widetilde G(r)=-\sum_k 4\alpha_k\,\mathrm{erf}(\tau_k r)/(3r)``

Eq. (A15) gives a contact density proportional to
``\sum_k \alpha_k\delta_{\tau_k}(r)``, where
``\tau_k^{-2}=\gamma_k^{-2}+\sigma_{12}^{-2}``.  This is not the product
``\alpha_s(r)\delta_{\sigma_{12}}(r)``: that older approximation both changed
the operator and introduced a non-analytic ``sqrt(x)`` dependence into the HO
quadrature variable ``x=(\beta r)^2``.
"""
function smeared_contact_kernel(
    params::GIParameters,
    masses::ConstituentMasses,
    r::Real,
)
    σ = contact_smearing_sigma(params, masses)
    return sum(
        α * delta_sigma_3d(r, inv(sqrt(inv(σ^2) + inv(γ^2)))) for
        (α, γ) in zip(ALPHA_COEFFS, ALPHA_GAMMAS)
    )
end

"""
Appendix A's post-A14 prescription places
`(m1*m2/(E1*E2))^(1/2 + epsilon_i)` on each side of a spin-dependent potential.
With `epsilon_i = 0`, the two-sided product turns the
nonrelativistic `1/(m1*m2)` strength into `1/(E1*E2)`.
"""
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
    side_exponent = gi_spin_dependent_side_exponent(params.factors.epsilon_c)
    B = momentum_relativization_matrix(m1, m2, side_exponent, p2_fact)
    kernel = Diagonal([smeared_contact_kernel(params, masses, ri) for ri in r])
    strength = (32 * π / (9 * m1 * m2)) * spin_dot(multiplicity)
    return Symmetric(strength * (B * kernel * B))
end

"""
    ho_contact_hyperfine_matrix(params, masses, L, multiplicity, beta, nbasis)

Native oscillator-basis contact-hyperfine matrix from Appendix A. The Gaussian
contact kernel, running coupling, and both relativistic momentum factors are
assembled directly in the HO basis.
"""
function ho_contact_hyperfine_matrix(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    multiplicity::Integer,
    beta::Real,
    nbasis::Integer,
)
    (L == 0 && multiplicity in (1, 3)) ||
        return Symmetric(zeros(Float64, nbasis, nbasis))
    params.factors.contact_momentum_sandwich || throw(ArgumentError(
        "native HO contact requires the Appendix-A momentum-sandwich prescription",
    ))
    radial = ho_momentum_sandwich_matrix(
        L,
        beta,
        nbasis,
        masses,
        params.factors.epsilon_c,
        r -> smeared_contact_kernel(params, masses, r),
        # The adaptive beta scan includes deliberately diffuse endpoints. A
        # 1e-8 radial-matrix tolerance keeps their narrow heavy-quark Gaussian
        # finite while remaining safely below the solver's 0.1 MeV energy gate
        # after the A15 mass prefactor; the accepted beta is tighter in practice.
        rtol = 1e-8,
    )
    strength = (32pi / (9 * masses.m1_GeV * masses.m2_GeV)) *
               spin_dot(multiplicity)
    return Symmetric(strength * radial)
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
    h = r[2] - r[1]
    expectation = radial_expect_udr(
        vector,
        r,
        h,
        (ri, i) -> begin
            smeared_contact_kernel(params, masses, ri)
        end,
    )
    (1.0 + params.factors.epsilon_c) *
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

"""
    resummed_channel_solution(params, masses, L, V; solver, nlevels)
        -> ChannelRadialSolution

Finite-difference diagnostic that diagonalizes the central Hamiltonian plus a
dense operator `V` on the solver's uniform mesh. Native HO callers use
[`fixed_channel_solution`](@ref); a bare mesh matrix is not an HO operator.

The returned wave has physical normalization `∫u² dr = 1`.
"""
function resummed_channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix;
    solver::RadialSolver = FiniteDifferenceSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    solver isa FiniteDifferenceSolver || throw(ArgumentError(
        "resummed_channel_solution accepts an FD mesh operator only; " *
        "use fixed_channel_solution for native HO assembly",
    ))
    return _resummed_channel_solution(
        solver, params, masses, L, V, nlevels)
end

function _resummed_channel_solution(
    solver::FiniteDifferenceSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix,
    nlevels::Integer,
)
    r, h = radial_grid(solver.ngrid, solver.rmax)
    size(V, 1) == length(r) || error("V must live on the (ngrid, rmax) mesh")
    hamiltonian, _ = relativistic_hamiltonian(params, masses, L; solver = solver)
    values, vectors = lowest_eigenpairs(
        Symmetric(Matrix(hamiltonian) + Matrix(V)), nlevels, solver)
    waves = physically_normalized_waves(Matrix(vectors), h)
    return ChannelRadialSolution(values, waves, r)
end

function contact_hyperfine_nonperturbative_levels(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    r::AbstractVector,
    nlevels::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    solver isa FiniteDifferenceSolver || throw(ArgumentError(
        "the overload carrying r is FD-only; omit r for native solver dispatch",
    ))
    if !params.factors.contact_momentum_sandwich || L != "S" || !(multiplicity in (1, 3)) || length(r) < 2
        return Float64[]
    end
    solution = _resummed_contact_solve(params, masses, L, multiplicity, r, nlevels, solver)
    return solution.eigenvalues_GeV
end

# Representation-free entry used by observables and audits.
function contact_hyperfine_nonperturbative_levels(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    nlevels::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    solution = contact_hyperfine_nonperturbative_states(
        params, masses, L, multiplicity, nlevels; solver = solver,
    )
    return isnothing(solution) ? Float64[] : solution.eigenvalues_GeV
end

"""
    contact_hyperfine_nonperturbative_states(params, masses, L, multiplicity, r, nlevels)
        -> Union{Nothing,ChannelRadialSolution}

Like [`contact_hyperfine_nonperturbative_levels`](@ref) but retains the native
radial waves of the S-wave Hamiltonian with the contact-hyperfine operator added
non-perturbatively. The singlet/triplet split of these waves is what makes the
`^1S_0` (e.g. `pi`) more compact than the `^3S_1` (e.g. `rho`) and drives the
Eq. (20)/(21) realistic-factor ratios. Returns `nothing` when the
non-perturbative contact path is inactive.
"""
function contact_hyperfine_nonperturbative_states(
    params::GIParameters,
    masses::ConstituentMasses,
    L::AbstractString,
    multiplicity::Integer,
    r::AbstractVector,
    nlevels::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    solver isa FiniteDifferenceSolver || throw(ArgumentError(
        "the overload carrying r is FD-only; omit r for native solver dispatch",
    ))
    if !params.factors.contact_momentum_sandwich || L != "S" || !(multiplicity in (1, 3)) || length(r) < 2
        return nothing
    end
    return _resummed_contact_solve(params, masses, L, multiplicity, r, nlevels, solver)
end

function contact_hyperfine_nonperturbative_states(
    params::GIParameters,
    masses::ConstituentMasses,
    L,
    multiplicity::Integer,
    nlevels::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    if !params.factors.contact_momentum_sandwich || L != "S" ||
       !(multiplicity in (1, 3))
        return nothing
    end
    J = multiplicity == 1 ? 0 : 1
    return fixed_channel_solution(
        params,
        masses,
        FineStructureMultiplet(String(L), multiplicity, J);
        solver = solver,
        terms = SpinTerms(
            contact_hyperfine = true,
            fine_structure = false,
            same_j_spin_orbit = false,
            tensor = false,
        ),
        nlevels = nlevels,
    )
end

# FD-only body for callers that already own an explicit operator grid.
function _resummed_contact_solve(
    params, masses, L, multiplicity, r, nlevels, solver::FiniteDifferenceSolver,
)
    h = r[2] - r[1]
    rmax = h * (length(r) + 1)
    rebuilt_r, _ = radial_grid(length(r), rmax)
    length(rebuilt_r) == length(r) || error("rebuilt S-wave grid changed length")
    V = contact_hyperfine_operator(params, masses, L, multiplicity, rebuilt_r)
    return resummed_channel_solution(
        params, masses, 0, Matrix(V);
        solver = with_mesh(solver, length(r), rmax), nlevels = nlevels,
    )
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
    wave::RadialWave,
)
    multiplet.L_label == "S" || return 0.0
    multiplet.multiplicity in (1, 3) || return 0.0
    expectation = radial_expect(
        wave,
        r -> smeared_contact_kernel(params, masses, r),
    )
    return (1.0 + params.factors.epsilon_c) *
           (32π / (9 * masses.m1_GeV * masses.m2_GeV)) *
           expectation * spin_dot(multiplet.multiplicity)
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
    wave::RadialWave,
)
    multiplet.L_label == "S" || return 0.0
    multiplet.multiplicity in (1, 3) || return 0.0
    expectation = radial_expect_momentum_sandwich(
        params,
        masses,
        0,
        wave,
        params.factors.epsilon_c,
        (r, _) -> smeared_contact_kernel(params, masses, r),
    )
    return (32π / (9 * masses.m1_GeV * masses.m2_GeV)) *
           expectation * spin_dot(multiplet.multiplicity)
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
    wave::RadialWave,
)
    if params.factors.contact_momentum_sandwich
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
    if params.factors.contact_momentum_sandwich
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
