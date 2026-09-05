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
    terms::SpinTerms,
    nbasis::Integer,
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
    solver::FiniteDifferenceSolver{:relativistic},
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    ::Nothing,
    terms::SpinTerms,
    nlevels::Integer,
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
    fine = if terms.fine_structure && params.fine_structure.enabled &&
              L > 0 && multiplet.multiplicity == 3
        fine_structure_grid_matrices(params, masses, multiplet.J, r, h; L = L)
    else
        zero_matrix = _zero_operator(n)
        (
            spin_orbit_vector = zero_matrix,
            spin_orbit_thomas = zero_matrix,
            spin_orbit = zero_matrix,
            tensor = zero_matrix,
            total = zero_matrix,
        )
    end
    return (
        central = central,
        contact = contact,
        spin_orbit_vector = fine.spin_orbit_vector,
        spin_orbit_thomas = fine.spin_orbit_thomas,
        spin_orbit = fine.spin_orbit,
        tensor = fine.tensor,
        fine_structure = fine.total,
        total = Symmetric(Matrix(central) + Matrix(contact) + Matrix(fine.total)),
        r = r,
    )
end

function _fixed_channel_matrices(
    solver::FiniteDifferenceSolver{:nonrelativistic},
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet,
    ::Nothing,
    terms::SpinTerms,
    nlevels::Integer,
)
    throw(ArgumentError(
        "the complete GI fixed-sector Hamiltonian is relativistic; " *
        "use channel_solution with this nonrelativistic solver as a central-only comparator",
    ))
end

function _fixed_channel_waves(
    ::FiniteDifferenceSolver,
    multiplet::FineStructureMultiplet,
    ::Nothing,
    vectors::AbstractMatrix,
    matrices,
)
    normalized = physically_normalized_waves(
        Matrix(vectors), matrices.r[2] - matrices.r[1],
    )
    return [MeshWave(view(normalized, :, n), matrices.r) for n in axes(normalized, 2)]
end

"""
    fixed_channel_solution(params, masses, multiplet; solver, terms, nlevels)

Diagonalize the complete radial Hamiltonian for one fixed `(L,S,J)` sector:
central + contact hyperfine + symmetric spin-orbit + diagonal tensor. HO and FD
return the same [`ChannelRadialSolution`](@ref) contract with native waves.
"""
function fixed_channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    multiplet::FineStructureMultiplet;
    solver::RadialSolver = FiniteDifferenceSolver(),
    terms::SpinTerms = SpinTerms(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    evaluate = function (beta, nbasis)
        matrices = _fixed_channel_matrices(
            solver, params, masses, multiplet, beta, terms, nbasis,
        )
        values, vectors = lowest_eigenpairs(Matrix(matrices.total), nlevels, solver)
        return (values = values, vectors = vectors, matrices = matrices)
    end
    result, waves, convergence = _fixed_channel_search(
        solver, multiplet, nlevels, evaluate,
    )
    return ChannelRadialSolution(
        result.values, waves; convergence = convergence,
    )
end

function _fixed_channel_search(
    solver::OscillatorSolver,
    multiplet::FineStructureMultiplet,
    nlevels::Integer,
    evaluate,
)
    best, waves, convergence = _oscillator_solution_search(
        solver, L_SYMBOLS[multiplet.L_label], nlevels, evaluate,
    )
    return best.result, waves, convergence
end

function _fixed_channel_search(
    solver::FiniteDifferenceSolver,
    multiplet::FineStructureMultiplet,
    nlevels::Integer,
    evaluate,
)
    result = evaluate(nothing, nlevels)
    waves = _fixed_channel_waves(
        solver, multiplet, nothing, result.vectors, result.matrices,
    )
    return result, waves, nothing
end
