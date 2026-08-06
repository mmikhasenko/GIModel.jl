using GIModel
using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PARAMETER_FILE = joinpath(ROOT, "data", "parameters.provisional.toml")

"""Rebuild today's concrete-Float64 parameter object with selected scalar overrides."""
function replace_continuous(
    params::GIParameters;
    b = params.potential.b,
    c = params.potential.c,
    sigma0 = params.smearing.sigma0,
    s = params.smearing.s,
)
    return GIParameters(
        params;
        potential = ConfinementPotential(params.potential; b = b, c = c),
        smearing = RelativisticSmearing(params.smearing; sigma0 = sigma0, s = s),
    )
end

function matrix_for(params, masses, solver)
    H, _ = GIModel.relativistic_hamiltonian(params, masses, 0; solver = solver)
    return Matrix(H)
end

function check_direction(name, x0, build; relative_steps = (1e-3, 3e-4, 1e-4), nlevels = 3)
    H0 = build(x0)
    F0 = eigen(Symmetric(H0))
    worst_difference = 0.0
    println("\n", name, " (x0 = ", x0, ")")
    println("  level    relative step       eigenvalue FD       v' (dH) v       abs difference")
    for relative_step in relative_steps
        step = relative_step * max(abs(x0), 1.0)
        Hplus = build(x0 + step)
        Hminus = build(x0 - step)
        dH = (Hplus - Hminus) / (2step)
        Eplus = eigvals(Symmetric(Hplus))
        Eminus = eigvals(Symmetric(Hminus))
        dE_fd = (Eplus - Eminus) / (2step)
        for level in 1:nlevels
            v = view(F0.vectors, :, level)
            dE_hf = dot(v, dH, v)
            worst_difference = max(worst_difference, abs(dE_fd[level] - dE_hf))
            @printf("  %5d    %13.1e    %17.10f    %17.10f    %14.3e\n",
                level, relative_step, dE_fd[level], dE_hf, abs(dE_fd[level] - dE_hf))
        end
    end
    @assert worst_difference < 2e-5 "$name: spectral derivative check failed"
    return nothing
end

function main()
    params, mass_table = load_parameters_and_quark_masses(PARAMETER_FILE)
    charm = mass_table["c"]
    masses = ConstituentMasses(charm, charm)
    # Small enough for a quick research probe, not a production-accuracy setting.
    solver = FiniteDifferenceSolver(ngrid = 160, rmax = 24.0, nlevels_per_channel = 3)

    H0 = matrix_for(params, masses, solver)
    F0 = eigen(Symmetric(H0))
    @printf("FD central charmonium S-wave probe: ngrid=%d, rmax=%.1f\n", solver.ngrid, solver.rmax)
    @printf("lowest central masses [GeV]: %s\n", join((@sprintf("%.8f", x) for x in F0.values[1:3]), ", "))
    @printf("smallest adjacent gap among those levels [GeV]: %.8f\n", minimum(diff(F0.values[1:3])))

    # On a fixed grid p2 is independent of every continuous physics parameter.
    r, h = GIModel.radial_grid(solver.ngrid, solver.rmax)
    p2a = Matrix(GIModel.p2_operator(params, masses.m1_GeV, 0, r, h))
    p2b = Matrix(GIModel.p2_operator(params, 2masses.m1_GeV, 0, r, h))
    @printf("max |p2(m)-p2(2m)|: %.3e\n", maximum(abs, p2a - p2b))
    @assert p2a == p2b

    check_direction("b [GeV^2]", params.potential.b,
        x -> matrix_for(replace_continuous(params; b = x), masses, solver))
    check_direction("c [GeV]", params.potential.c,
        x -> matrix_for(replace_continuous(params; c = x), masses, solver))
    check_direction("sigma0 [GeV]", params.smearing.sigma0,
        x -> matrix_for(replace_continuous(params; sigma0 = x), masses, solver))
    check_direction("s", params.smearing.s,
        x -> matrix_for(replace_continuous(params; s = x), masses, solver))
    check_direction("common charm mass [GeV]", charm,
        x -> matrix_for(params, ConstituentMasses(x, x), solver))

    step = 1e-5
    dHdc = (
        matrix_for(replace_continuous(params; c = params.potential.c + step), masses, solver) -
        matrix_for(replace_continuous(params; c = params.potential.c - step), masses, solver)
    ) / (2step)
    @printf("\nmax |dH/dc - I|: %.3e\n", maximum(abs, dHdc - I))
    @assert maximum(abs, dHdc - I) < 1e-7

    dHdb = (
        matrix_for(replace_continuous(params; b = params.potential.b + step), masses, solver) -
        matrix_for(replace_continuous(params; b = params.potential.b - step), masses, solver)
    ) / (2step)
    unit_b = replace_continuous(params; b = 1.0, c = 0.0)
    dHdb_exact = Diagonal([
        GIModel.smeared_confinement_S_closed(
            unit_b, masses.m1_GeV, masses.m2_GeV, ri,
        ) for ri in r
    ])
    @printf("max |dH/db - analytic|: %.3e\n", maximum(abs, dHdb - dHdb_exact))
    @assert maximum(abs, dHdb - dHdb_exact) < 1e-7
    return nothing
end

main()
