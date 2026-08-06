# Public API (exported from GIModel.jl):
#   channel_solution

# Minimum grid points across the ground state's RMS radius. Degradation is
# smooth, not a cliff, so this is calibrated against measured error rather than
# guessed. Ground-state energy error vs a converged (ngrid = 2400) reference:
#
#   pts across   19.4   15.4   10.8   10.3    6.9    5.7
#   1S error/MeV -0.53  -0.98  -2.49  -1.91  -7.54  -8.11
#
# Spin-dependent quantities degrade much faster than eigenvalues: at 6.9 points
# the quarkonium hyperfine splitting collapses (0.0025 GeV, breaking an
# otherwise monotone trend), while at 10.8 it still tracks. The break therefore
# sits between ~10 and ~7 points. Nine keeps the genuinely broken cases loud and
# stays silent both for the default grid up to m_Q ~ 15 GeV and for the coarser
# ngrid = 220 grid the Table III mixing audit deliberately uses for bottomonium
# (10.3 points, ~2 MeV — imprecise, not broken).
const MIN_POINTS_ACROSS_STATE = 9

# The FD grid is fixed (ngrid, rmax), so a sufficiently compact state simply
# falls between grid points and is silently under-resolved — the eigenvalues
# still come back, they are just wrong. Warn instead.
function _warn_if_underresolved(values, vectors, r, h, masses, L)
    isempty(values) && return nothing
    u = @view vectors[:, 1]
    nrm = sum(abs2, u) * h
    nrm > 0 || return nothing
    rms = sqrt(sum(abs2.(u) .* r .^ 2) * h / nrm)
    pts = rms / h
    pts < MIN_POINTS_ACROSS_STATE || return nothing
    @warn """
    Radial grid too coarse for this state: its RMS radius spans only \
    $(round(pts; digits = 1)) grid points (want >= $MIN_POINTS_ACROSS_STATE). \
    Expect several MeV of error on the eigenvalues and much worse on \
    spin-dependent quantities, which degrade faster. Increase `ngrid` (or \
    reduce `rmax`) for constituent masses well above bottomonium.""" m1 =
        masses.m1_GeV m2 = masses.m2_GeV L rms_radius = rms grid_spacing = h maxlog = 1
    return nothing
end

"""
    channel_solution(params, masses, L; solver, nlevels) -> (values, waves, r)

Solve the radial equation for one orbital channel. Together with
[`resummed_channel_solution`](@ref) this is where all the numerics in the model
live; everything downstream consumes the returned `u(r)`.

The `solver` picks the method — [`FiniteDifferenceSolver`](@ref) or
[`OscillatorSolver`](@ref) — and both return the same physical quantity in the
same convention: eigenvalues in GeV, and waves normalized `∫u² dr = 1` on the
returned mesh `r`.
"""
function channel_solution(
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer;
    solver::RadialSolver = FiniteDifferenceSolver(),
    nlevels::Integer = solver.nlevels_per_channel,
)
    # `nlevels` is a per-call quantity rather than a solver setting -- the
    # annihilation wave cache legitimately asks for fewer levels than the solver
    # would otherwise keep -- so it stays a keyword here.
    return _channel_solution(solver, params, masses, L, nlevels)
end

function _channel_solution(
    solver::FiniteDifferenceSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    nlevels::Integer,
)
    kinetic = solver.kinetic
    hamiltonian, r = if kinetic == :relativistic
        relativistic_hamiltonian(params, masses, L; solver = solver)
    elseif kinetic == :nonrelativistic
        nonrelativistic_hamiltonian(params, masses, L; solver = solver)
    else
        error("unknown kinetic mode: $kinetic")
    end
    values, vectors = lowest_eigenpairs(hamiltonian, nlevels; solver = solver)
    _warn_if_underresolved(values, vectors, r, length(r) > 1 ? r[2] - r[1] : 0.0, masses, L)
    waves = physically_normalized_waves(Matrix(vectors), length(r) > 1 ? r[2] - r[1] : 1.0)
    if kinetic == :relativistic
        values, waves, r
    else
        values .+ (masses.m1_GeV + masses.m2_GeV), waves, r
    end
end

_channel_solution(
    solver::OscillatorSolver,
    params::GIParameters,
    masses::ConstituentMasses,
    L::Integer,
    nlevels::Integer,
) = oscillator_channel_solution(params, masses, L; nlevels = nlevels, solver = solver)

function solve_channel(args...; kwargs...)
    values, _vectors, _r = channel_solution(args...; kwargs...)
    values
end
