# Public API (exported from GIModel.jl):
#   spin_dot, contact_smearing_sigma

"""
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
    side_exponent = gi_spin_dependent_side_exponent(params.factors.epsilon_c)
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
        -> (values, waves, r)

Diagonalize the central Hamiltonian **plus** a spin-dependent operator `V`,
resumming `V` rather than treating it at first order. This is one of the two
places in the model where numerics actually happen (the other is the central
solve, [`channel_solution`](@ref)); everything downstream consumes the returned
`u(r)`.

Two methods, one per solver, both returning the same physical quantity:

  - [`FiniteDifferenceSolver`](@ref) — build `H` on the uniform mesh, add `V`,
    diagonalize (this method).
  - [`OscillatorSolver`](@ref) — project `V` into the finite oscillator space and
    diagonalize there (the paper's own method; see `harmonic_oscillator_basis.jl`).

`V` lives on the mesh either way, so the solver's `ngrid`/`rmax` must match it.

`V` is a dense operator on the same uniform mesh, e.g.
[`contact_hyperfine_operator`](@ref GIModel.contact_hyperfine_operator).

**Normalization is physical for both methods: `∫u² dr = 1`.** The finite-
difference path used to return Euclidean eigenvectors (`Σu² = 1`) instead,
differing from the oscillator path by exactly `√h`; consumers had to know which
convention they held, and anything quadratic in `u` given the wrong one was off
by `h`. Both now agree.
"""
function resummed_channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    V::AbstractMatrix;
    solver::RadialSolver = FiniteDifferenceSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
    ngrid::Integer = solver.ngrid,
    rmax::Real = solver.rmax,
)
    return _resummed_channel_solution(
        with_mesh(solver, ngrid, rmax), params, masses, L, V, nlevels)
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
        Symmetric(Matrix(hamiltonian) + Matrix(V)), nlevels)
    waves = physically_normalized_waves(Matrix(vectors), h)
    return collect(values), waves, collect(Float64, r)
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
    if !params.factors.contact_momentum_sandwich || L != "S" || !(multiplicity in (1, 3)) || length(r) < 2
        return Float64[]
    end
    levels, _waves, _r =
        _resummed_contact_solve(params, masses, L, multiplicity, r, nlevels, solver)
    return levels
end

"""
    contact_hyperfine_nonperturbative_states(params, masses, L, multiplicity, r, nlevels)
        -> (levels::Vector{Float64}, vectors::Matrix{Float64}, r::Vector{Float64})

Like [`contact_hyperfine_nonperturbative_levels`](@ref) but also returns the
eigenvectors (reduced radial waves `u(r)`, one per column) of the S-wave
Hamiltonian with the contact-hyperfine operator added non-perturbatively, plus
the rebuilt grid. The singlet/triplet split of these waves is what makes the
`^1S_0` (e.g. `pi`) more compact than the `^3S_1` (e.g. `rho`) and drives the
Eq. (20)/(21) realistic-factor ratios. Returns `(Float64[], zeros(0,0),
Float64[])` when the non-perturbative contact path is inactive.
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
    if !params.factors.contact_momentum_sandwich || L != "S" || !(multiplicity in (1, 3)) || length(r) < 2
        return Float64[], zeros(Float64, 0, 0), Float64[]
    end
    return _resummed_contact_solve(params, masses, L, multiplicity, r, nlevels, solver)
end

# Shared body: build the contact operator for this (L, multiplicity) and hand it
# to the solver-dispatched resummed solve. Callers hold `r` rather than a mesh
# specification, so the solver is re-pointed at exactly that mesh -- the operator
# `V` is built on it and `resummed_channel_solution` rejects a mismatch.
function _resummed_contact_solve(params, masses, L, multiplicity, r, nlevels, solver)
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

# Both entry points above are solver-generic. The guard (S-wave, multiplicity 1
# or 3, momentum sandwich on) is physics; building the contact operator is
# method-free (`p2_operator` returns the same mesh operator either way); and the
# solve itself is `resummed_channel_solution`, which dispatches on the solver.
#
# Until U2 unified that solve, this function inlined the finite-difference one,
# so it had to be FD-only -- and its catch-all returned empty for any other
# method, which callers read as "unavailable" and replaced with first-order PT
# (the light 1S0 came out 0.2842 GeV instead of 0.0950). U1 turned that into a
# loud failure; with the paths unified there is nothing left to fail about.

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
