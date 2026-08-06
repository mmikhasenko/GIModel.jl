using GIModel
using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PARAMETER_FILE = joinpath(ROOT, "data", "parameters.provisional.toml")

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

function matrix_for(params, masses, r, h, beta, nbasis)
    H, _ = GIModel.oscillator_hamiltonian_for_beta(
        params, masses, 0, r, h, beta; nbasis = nbasis,
    )
    return Matrix(H)
end

function check_direction(name, x0, build; relative_steps = (3e-4, 1e-4), nlevels = 3)
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
            difference = abs(dE_fd[level] - dE_hf)
            worst_difference = max(worst_difference, difference)
            @printf("  %5d    %13.1e    %17.10f    %17.10f    %14.3e\n",
                level, relative_step, dE_fd[level], dE_hf, difference)
        end
    end
    @assert worst_difference < 2e-5 "$name: fixed-beta spectral derivative check failed"
    return nothing
end

function main()
    params, mass_table = load_parameters_and_quark_masses(PARAMETER_FILE)
    charm = mass_table["c"]
    masses = ConstituentMasses(charm, charm)

    # This intentionally checks one smooth HO branch. It does not differentiate
    # the discrete beta-grid selection used by oscillator_channel_solution.
    beta = 0.65
    nbasis = 8
    r, h = GIModel.radial_grid(160, 24.0)
    build(params_here, masses_here) =
        matrix_for(params_here, masses_here, r, h, beta, nbasis)

    H0 = build(params, masses)
    F0 = eigen(Symmetric(H0))
    @printf("fixed-beta HO charmonium S-wave probe: beta=%.2f, nbasis=%d\n", beta, nbasis)
    @printf("lowest central masses [GeV]: %s\n", join((@sprintf("%.8f", x) for x in F0.values[1:3]), ", "))
    @printf("smallest adjacent gap among those levels [GeV]: %.8f\n", minimum(diff(F0.values[1:3])))

    check_direction("b [GeV^2]", params.potential.b,
        x -> build(replace_continuous(params; b = x), masses))
    check_direction("c [GeV]", params.potential.c,
        x -> build(replace_continuous(params; c = x), masses))
    check_direction("sigma0 [GeV]", params.smearing.sigma0,
        x -> build(replace_continuous(params; sigma0 = x), masses))
    check_direction("s", params.smearing.s,
        x -> build(replace_continuous(params; s = x), masses))
    check_direction("common charm mass [GeV]", charm,
        x -> build(params, ConstituentMasses(x, x)))

    step = 1e-5
    dHdc = (
        build(replace_continuous(params; c = params.potential.c + step), masses) -
        build(replace_continuous(params; c = params.potential.c - step), masses)
    ) / (2step)
    @printf("\nmax |dH/dc - I|: %.3e\n", maximum(abs, dHdc - I))
    @assert maximum(abs, dHdc - I) < 1e-7
    return nothing
end

main()
