# One radial solve per fixed spectroscopic sector. The solver chooses only the
# representation; both methods assemble the same central + contact + symmetric
# spin-orbit + diagonal tensor Hamiltonian before diagonalization.

function _zero_operator(n::Integer)
    return Symmetric(zeros(Float64, n, n))
end

function _fixed_channel_matrices(
    solver::OscillatorSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    beta::Real,
    nbasis::Integer,
    terms::SpinTerms,
)
    L = L_SYMBOLS[multiplet.L_label]
    central = oscillator_central_matrix(params, masses, L, beta, nbasis)
    contact = terms.contact_hyperfine ? ho_contact_hyperfine_matrix(
        params, masses, L, multiplet.multiplicity, beta, nbasis,
    ) : _zero_operator(nbasis)
    fine = terms.fine_structure ? ho_fine_structure_matrices(
        params,
        masses,
        L,
        multiplet.multiplicity,
        multiplet.J,
        beta,
        nbasis,
    ) : (
        spin_orbit_vector = _zero_operator(nbasis),
        spin_orbit_thomas = _zero_operator(nbasis),
        spin_orbit = _zero_operator(nbasis),
        tensor = _zero_operator(nbasis),
        total = _zero_operator(nbasis),
    )
    return (
        central = central,
        contact = contact,
        spin_orbit_vector = fine.spin_orbit_vector,
        spin_orbit_thomas = fine.spin_orbit_thomas,
        spin_orbit = fine.spin_orbit,
        tensor = fine.tensor,
        fine_structure = fine.total,
        total = Symmetric(Matrix(central) + Matrix(contact) + Matrix(fine.total)),
    )
end

function _fixed_channel_matrices(
    solver::FiniteDifferenceSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    terms::SpinTerms,
)
    L = L_SYMBOLS[multiplet.L_label]
    central, r = relativistic_hamiltonian(params, masses, L; solver = solver)
    h = r[2] - r[1]
    n = length(r)
    contact = if terms.contact_hyperfine && L == 0 &&
                 multiplet.multiplicity in (1, 3)
        if params.factors.contact_momentum_sandwich
            contact_hyperfine_operator(
                params, masses, multiplet.L_label, multiplet.multiplicity, r,
            )
        else
            sigma = contact_smearing_sigma(params, masses)
            strength = (1 + params.factors.epsilon_c) *
                       (32pi / (9 * masses.m1_GeV * masses.m2_GeV)) *
                       spin_dot(multiplet.multiplicity)
            Diagonal([
                strength * alpha_s_r(ri) * delta_sigma_3d(ri, sigma) for ri in r
            ])
        end
    else
        _zero_operator(n)
    end
    fine_total = if terms.fine_structure && params.fine_structure.enabled &&
                    L > 0 && multiplet.multiplicity == 3
        fine_structure_grid_operator(params, masses, multiplet.J, r, h; L = L)
    else
        _zero_operator(n)
    end
    # The FD operator builder currently exposes the summed fine-structure
    # matrix. Component reporting remains available through wave expectations.
    return (
        central = central,
        contact = contact,
        spin_orbit_vector = _zero_operator(n),
        spin_orbit_thomas = _zero_operator(n),
        spin_orbit = _zero_operator(n),
        tensor = _zero_operator(n),
        fine_structure = fine_total,
        total = Symmetric(Matrix(central) + Matrix(contact) + Matrix(fine_total)),
        r = r,
    )
end

"""
    fixed_channel_solution(params, masses, multiplet; solver, terms, nlevels)

Diagonalize the complete radial Hamiltonian for one fixed `(L,S,J)` sector:
central + contact hyperfine + symmetric spin-orbit + diagonal tensor. HO and FD
return the same [`ChannelRadialSolution`](@ref) contract with native waves.
"""
function _fixed_channel_solution(
    solver::OscillatorSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet;
    terms::SpinTerms = SpinTerms(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    nbasis = max(solver.nbasis, nlevels + 4)
    L = L_SYMBOLS[multiplet.L_label]
    best = nothing
    for beta in solver.beta_grid
        matrices = _fixed_channel_matrices(
            solver, params, masses, multiplet, beta, nbasis, terms,
        )
        values, vectors = lowest_eigenpairs(Matrix(matrices.total), nlevels, solver)
        if isnothing(best) || values[end] < best.values[end]
            best = (beta = beta, values = values, vectors = vectors)
        end
    end
    _warn_if_beta_railed(best.beta, solver, masses, L)
    waves = [
        OscillatorWave(L, best.beta, view(best.vectors, :, n)) for n in 1:nlevels
    ]
    return ChannelRadialSolution(best.values, waves)
end

function _fixed_channel_solution(
    solver::FiniteDifferenceSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet;
    terms::SpinTerms = SpinTerms(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    matrices = _fixed_channel_matrices(solver, params, masses, multiplet, terms)
    values, vectors = lowest_eigenpairs(Matrix(matrices.total), nlevels, solver)
    waves = physically_normalized_waves(Matrix(vectors), matrices.r[2] - matrices.r[1])
    return ChannelRadialSolution(values, waves, matrices.r)
end

function fixed_channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet;
    solver::RadialSolver = FiniteDifferenceSolver(),
    terms::SpinTerms = SpinTerms(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    return _fixed_channel_solution(
        solver, params, masses, multiplet; terms = terms, nlevels = nlevels,
    )
end
