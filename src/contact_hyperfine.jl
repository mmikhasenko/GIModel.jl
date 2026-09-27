# Public API (exported from GIModel.jl):
#   spin_dot, contact_smearing_sigma

@doc raw"""
    contact_smearing_sigma(params, m1, m2) -> Float64
    contact_smearing_sigma(params, masses::ConstituentMasses) -> Float64

The Appendix A (A9) universal smearing width ``\\sigma(m_1, m_2)`` in GeV,

    σ² = σ₀² (1/2 + 1/2 [4 m₁m₂/(m₁+m₂)²]⁴) + s² (2 m₁m₂/(m₁+m₂))²

built from the Table II inputs `σ₀` and `s`. It specifies the mass dependence
of the smearing; the interaction also carries explicit mass denominators and
relativistic momentum factors. The potential parameters (`b`, `c`) and
[`alpha_s_q`](@ref) are flavor independent.
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

@doc raw"""
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

"""
Internal definition of the smeared spin-spin interaction. Orbital momentum is
an input to its representation, never a selection rule: A15 acts for every L.
The explicit `sandwich` override supports the historical local approximation
without duplicating its strength or kernel in solvers and observables.
"""
struct ContactHyperfine{P<:GIParameters}
    params::P
    masses::ConstituentMasses
    multiplicity::Int
    sandwich::Bool
    function ContactHyperfine(params, masses, multiplicity;
                              sandwich = params.factors.contact_momentum_sandwich)
        multiplicity in (1, 3) || throw(ArgumentError(
            "contact hyperfine requires q-qbar spin multiplicity 1 or 3",
        ))
        new{typeof(params)}(params, masses, multiplicity, sandwich)
    end
end

contact_kernel(op::ContactHyperfine, r) = smeared_contact_kernel(op.params, op.masses, r)
contact_strength(op::ContactHyperfine) =
    (32pi / (9 * op.masses.m1_GeV * op.masses.m2_GeV)) *
    spin_dot(op.multiplicity) * (op.sandwich ? 1.0 : 1 + op.params.factors.epsilon_c)

# Representation adapters contain numerical operations only. Neither may
# suppress an orbital sector or supply an implicit L=0.
function contact_matrix(op::ContactHyperfine, L::Integer, r::AbstractVector)
    L >= 0 || throw(ArgumentError("contact matrix requires L >= 0"))
    length(r) >= 2 || throw(ArgumentError("contact matrix requires at least two mesh points"))
    h = r[2] - r[1]
    physical_u_norm(r, h, ones(length(r))) # validates the uniform radial mesh
    K = Diagonal([contact_kernel(op, ri) for ri in r])
    if !op.sandwich
        return Symmetric(Matrix(contact_strength(op) * K))
    end
    masses = op.masses
    p2_fact = eigen(p2_operator(op.params, masses.m1_GeV, L, r, h))
    B = momentum_relativization_matrix(
        masses.m1_GeV, masses.m2_GeV,
        gi_spin_dependent_side_exponent(op.params.factors.epsilon_c), p2_fact,
    )
    return Symmetric(contact_strength(op) * (B * K * B))
end

function contact_matrix(op::ContactHyperfine, L::Integer, beta::Real, nbasis::Integer)
    L >= 0 || throw(ArgumentError("contact matrix requires L >= 0"))
    kernel = r -> contact_kernel(op, r)
    # Diffuse beta-scan endpoints need 1e-8 matrix accuracy; the accepted
    # basis is tighter in practice. Solving and reporting use the same rule.
    radial = op.sandwich ? ho_momentum_sandwich_matrix(
        L, beta, nbasis, op.masses, op.params.factors.epsilon_c, kernel; rtol = 1e-8,
    ) : ho_operator_matrix(L, beta, nbasis, kernel; rtol = 1e-8)
    return Symmetric(contact_strength(op) * radial)
end

function contact_expectation(op::ContactHyperfine, L::Integer, wave::MeshWave)
    return euclidean_expectation(wave.u, contact_matrix(op, L, wave.r))
end

function contact_expectation(op::ContactHyperfine, L::Integer, wave::OscillatorWave)
    L == wave.L || throw(ArgumentError("contact expectation: sector L does not match wave L"))
    return euclidean_expectation(
        wave.coefficients, contact_matrix(op, L, wave.beta, length(wave.coefficients)),
    )
end

# Compatibility entry points delegate to the same definition and adapters.
function contact_hyperfine_operator(params::GIParameters, masses::ConstituentMasses,
                                    L::AbstractString, multiplicity::Integer, r::AbstractVector)
    return contact_matrix(ContactHyperfine(params, masses, multiplicity; sandwich = true),
                          L_SYMBOLS[L], r)
end

function ho_contact_hyperfine_matrix(params::GIParameters, masses::ConstituentMasses,
                                     L::Integer, multiplicity::Integer, beta::Real, nbasis::Integer)
    return contact_matrix(ContactHyperfine(params, masses, multiplicity), L, beta, nbasis)
end

function _contact_hyperfine_shift_diagonal(params::GIParameters, masses::ConstituentMasses,
                                          L::AbstractString, multiplicity::Integer,
                                          vector::AbstractVector, r::AbstractVector)
    return contact_expectation(ContactHyperfine(params, masses, multiplicity; sandwich = false),
                               L_SYMBOLS[L], MeshWave(vector, r))
end

function _contact_hyperfine_shift_momentum_sandwich_diagonal(
    params::GIParameters, masses::ConstituentMasses, L::AbstractString,
    multiplicity::Integer, vector::AbstractVector, r::AbstractVector,
)
    return contact_expectation(ContactHyperfine(params, masses, multiplicity; sandwich = true),
                               L_SYMBOLS[L], MeshWave(vector, r))
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
    solution = contact_hyperfine_nonperturbative_states(
        params, masses, L, multiplicity, r, nlevels; solver,
    )
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
    return solution.eigenvalues_GeV
end

"""
    contact_hyperfine_nonperturbative_states(params, masses, L, multiplicity, r, nlevels)
        -> ChannelRadialSolution

Like `contact_hyperfine_nonperturbative_levels` but retains the native
radial waves of the fixed-L Hamiltonian with the contact-hyperfine operator added
non-perturbatively. The singlet/triplet split of these waves is what makes the
`^1S_0` (e.g. `pi`) more compact than the `^3S_1` (e.g. `rho`) and drives the
Eq. (20)/(21) realistic-factor ratios. All orbital partial waves use the
same fixed-sector solve. Invalid spins or meshes throw instead of returning
an empty result that could trigger a silent fallback.
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
    length(r) >= 2 || throw(ArgumentError("contact solve requires at least two mesh points"))
    h = r[2] - r[1]
    physical_u_norm(r, h, ones(length(r)))
    rmax = h * (length(r) + 1)
    rebuilt_r, _ = radial_grid(length(r), rmax)
    isapprox(r, rebuilt_r; rtol = 1e-10, atol = 1e-12) || throw(ArgumentError(
        "contact solve requires the interior Dirichlet grid r[i] = i*h",
    ))
    return contact_hyperfine_nonperturbative_states(
        params, masses, L, multiplicity, nlevels;
        solver = with_mesh(solver, length(r), rmax),
    )
end

function contact_hyperfine_nonperturbative_states(
    params::GIParameters,
    masses::ConstituentMasses,
    L,
    multiplicity::Integer,
    nlevels::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    ContactHyperfine(params, masses, multiplicity) # validate before solving
    orbital = L_SYMBOLS[String(L)]
    # Contact is independent of J; choose one allowed representative channel.
    J = multiplicity == 1 ? orbital : orbital + 1
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

"""
First-order local contact approximation, for any `L`. This explicit diagnostic
uses the historical `(1 + epsilon_c)` coefficient; production callers use
`contact_hyperfine_shift_active` to select the configured prescription.

For the reduced radial wavefunction `u(r)`, with `∫|u|² dr = 1`, a radial
smeared contact kernel `K(r)` has expectation `∫ |u(r)|² K(r) dr` for every L.
There is no extra `4π` factor. Solving and reporting use the same native matrix.

Uses only [`FineStructureMultiplet`](@ref).`L_label` and `.multiplicity`; `.J` is unused (same multiplet object as fine-structure).
"""
function contact_hyperfine_shift(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    wave::RadialWave,
)
    return contact_expectation(
        ContactHyperfine(params, masses, multiplet.multiplicity; sandwich = false),
        L_SYMBOLS[multiplet.L_label], wave,
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
    wave::RadialWave,
)
    return contact_expectation(
        ContactHyperfine(params, masses, multiplet.multiplicity; sandwich = true),
        L_SYMBOLS[multiplet.L_label], wave,
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
    wave::RadialWave,
)
    return contact_expectation(ContactHyperfine(params, masses, multiplet.multiplicity),
                               L_SYMBOLS[multiplet.L_label], wave)
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
    return contact_expectation(ContactHyperfine(params, ConstituentMasses(m1, m2), multiplicity),
                               L_SYMBOLS[L], MeshWave(vector, r))
end
