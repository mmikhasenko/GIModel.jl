@testset "parameter loading" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    @test mq["c"] ≈ 1.628
    @test mq["b"] ≈ 4.977
    @test params.potential.b ≈ 0.18
    @test params.central isa AppendixAMomentumSandwich
    @test params.factors.epsilon_c ≈ -0.168
    @test params.factors.epsilon_t ≈ 0.025
    @test params.factors.epsilon_so_vector ≈ -0.035
    @test params.factors.epsilon_so_scalar ≈ 0.055
    @test params.fine_structure.enabled == true
    @test isconcretetype(typeof(params))
    @test fieldtype(typeof(params), :central) === AppendixAMomentumSandwich

    varied = GIParameters(params;
        potential = ConfinementPotential(params.potential; b = 0.19),
        smearing = RelativisticSmearing(params.smearing; s = 1.6),
    )
    @test varied.potential.b == 0.19
    @test varied.potential.c == params.potential.c
    @test varied.smearing.s == 1.6
    @test varied.central === params.central
    varied_factors = RelativisticFactors(params.factors; epsilon_c = -0.2)
    @test varied_factors.epsilon_c == -0.2
    @test varied_factors.epsilon_t == params.factors.epsilon_t
    varied_fine = FineStructure(params.fine_structure; enabled = false)
    @test !varied_fine.enabled
    varied_annihilation = AnnihilationAmplitudes(params.annihilation; s1_A = 2.6)
    @test varied_annihilation.s1_A == 2.6
    @test varied_annihilation.p1_A_np == params.annihilation.p1_A_np
end

@testset "central numerical boundaries are inferred" begin
    params, mq = load_parameters_and_quark_masses(
        joinpath(root, "data", "parameters.provisional.toml"),
    )
    masses = ConstituentMasses(mq["c"], mq["c"])
    r = collect(range(0.1, 4.0; length = 24))
    @test (@inferred GIModel.potential_diagonal(
        params, masses.m1_GeV, masses.m2_GeV, r,
    )) isa Vector{Float64}

    fd = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0, nlevels_per_channel = 2)
    Hfd, _ = @inferred GIModel.relativistic_hamiltonian(params, masses, 0; solver = fd)
    @test Hfd isa Symmetric{Float64,Matrix{Float64}}
    fd_solution = @inferred GIModel._channel_solution(
        fd, params, masses, 0, 2,
    )
    @test fd_solution isa ChannelRadialSolution{MeshWave}

    ho = OscillatorSolver(
        nbasis = 8,
        beta_grid = [0.65],
        converge = false,
        nlevels_per_channel = 2,
    )
    ho_solution = @inferred GIModel._channel_solution(
        ho, params, masses, 0, 2,
    )
    @test ho_solution isa ChannelRadialSolution{OscillatorWave}

    # Solver settings compare by value, so spectra solved with separately
    # constructed but identical solvers can be combined.
    @test OscillatorSolver() == OscillatorSolver()
    @test hash(OscillatorSolver()) == hash(OscillatorSolver())
    @test OscillatorSolver() != OscillatorSolver(nbasis = 32)
    @test OscillatorSolver() != OscillatorSolver(beta_grid = [0.7])

    matrix = Symmetric([2.0 1.0; 1.0 2.0])
    for eigensolver in (:full, :krylov)
        solver = FiniteDifferenceSolver(eigensolver = eigensolver)
        @test typeof(solver) === FiniteDifferenceSolver{:relativistic,eigensolver}
        values, vectors = @inferred GIModel.lowest_eigenpairs(matrix, 1, solver)
        @test values isa Vector{Float64}
        @test vectors isa Matrix{Float64}
    end
end

@testset "continuous parameter leaves preserve numeric types" begin
    potential32 = ConfinementPotential(b = 0.18f0, c = -0.253f0)
    smearing32 = RelativisticSmearing(sigma0 = 1.8f0, s = 1.55f0)
    masses32 = ConstituentMasses(1.628f0, 1.628f0)
    @test potential32 isa ConfinementPotential{Float32}
    @test smearing32 isa RelativisticSmearing{Float32}
    @test masses32 isa ConstituentMasses{Float32}
    @test ConfinementPotential(b = 1, c = 0.5) isa ConfinementPotential{Float64}

    params, mq = load_parameters_and_quark_masses(
        joinpath(root, "data", "parameters.provisional.toml"),
    )
    big_params = GIParameters(
        params;
        potential = ConfinementPotential(
            b = BigFloat(params.potential.b),
            c = BigFloat(params.potential.c),
        ),
        smearing = RelativisticSmearing(
            sigma0 = BigFloat(params.smearing.sigma0),
            s = BigFloat(params.smearing.s),
        ),
    )
    big_masses = ConstituentMasses(BigFloat(mq["c"]), BigFloat(mq["c"]))
    solver = FiniteDifferenceSolver(ngrid = 24, rmax = 8.0)
    Hfd, _ = GIModel.relativistic_hamiltonian(big_params, big_masses, 0; solver = solver)
    @test eltype(Hfd) === BigFloat

    r, h = GIModel.radial_grid(24, 8.0)
    Hho, _ = GIModel.oscillator_hamiltonian_for_beta(
        big_params, big_masses, 0, r, h, 0.65; nbasis = 6,
    )
    @test eltype(Hho) === BigFloat
end

@testset "baseline solver shape" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    cc = solve_sector(
        params, mq["c"];
        maxn = 4,
        solver = FiniteDifferenceSolver(ngrid = 250, rmax = 20.0),
    )
    bb = solve_sector(
        params, mq["b"];
        maxn = 4,
        solver = FiniteDifferenceSolver(ngrid = 250, rmax = 16.0),
    )

    @test cc[(1, "S")] < cc[(2, "S")] < cc[(3, "S")]
    @test bb[(1, "S")] < bb[(2, "S")] < bb[(3, "S")]
    @test cc[(1, "S")] < cc[(1, "P")] < cc[(1, "D")]
    @test bb[(1, "S")] < bb[(1, "P")] < bb[(1, "D")]
end

@testset "Krylov eigensolver matches full eigensolver" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    for kinetic in (:relativistic, :nonrelativistic)
        full = channel_solution(
            params, ConstituentMasses(m, m), 0;
            nlevels = 3,
            solver = FiniteDifferenceSolver(ngrid = 120, rmax = 16.0, kinetic = kinetic, eigensolver = :full),
        )
        krylov = channel_solution(
            params, ConstituentMasses(m, m), 0;
            nlevels = 3,
            solver = FiniteDifferenceSolver(ngrid = 120, rmax = 16.0, kinetic = kinetic, eigensolver = :krylov),
        )
        @test krylov.eigenvalues_GeV ≈ full.eigenvalues_GeV rtol = 1e-10 atol = 1e-10
    end
end

@testset "central: Coulomb+confinement = pointwise V" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    for r0 in (0.15, 0.4, 1.2, 3.0)
        a = GIModel.central_potential(r0, params)
        b = GIModel.static_coulomb_G(r0, params) + GIModel.static_confinement_S(r0, params)
        @test a ≈ b rtol = 1e-12
    end
end

@testset "Coulomb central derivative matches finite difference" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    Vg(r) = -4.0 * GIModel.alpha_s_r(r) / (3.0 * r)
    for r0 in (0.08, 0.15, 0.4, 1.2, 3.0)
        δ = 1e-6 * max(1.0, r0)
        num = (Vg(r0 + δ) - Vg(r0 - δ)) / (2δ)
        ana = GIModel.dV_coul_central_dr(r0, params)
        @test ana ≈ num rtol = 2e-7 atol = 1e-9
    end
end

@testset "Gaussian contact regulator is 3D-normalized" begin
    # δ_σ(r) is implemented as a 3D density. Its defining normalization is
    #   ∫ d³r δ_σ(r) = 1  ⇔  4π ∫ r² δ_σ(r) dr = 1.
    # We check it on a finite mesh that captures essentially all support.
    for σ in (0.2, 0.5, 1.0, 2.0, 5.0)
        rmax = 12.0 / σ
        n = 4000
        h = rmax / n
        s = 0.0
        for i = 1:n
            r = (i - 0.5) * h
            s += 4π * r^2 * GIModel.delta_sigma_3d(r, σ) * h
        end
        @test s ≈ 1.0 rtol = 2e-6 atol = 2e-6
    end
end

@testset "central_potential_path (default = GI momentum sandwich)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    @test params isa GIParameters
    p = GIModel.central_potential_path(params)
    @test p.name == "appendix_a_momentum_sandwich"
    @test params.central isa AppendixAMomentumSandwich
    @test params.central === central_potential_method("appendix_a_momentum_sandwich")
    @test params.annihilation.p1_A_np ≈ 0.50
    @test params.annihilation.p1_m_eta ≈ 0.548
    @test params.annihilation.p2_A_np ≈ 0.55
    @test params.annihilation.p2_M0 ≈ 1.17
    @test params.annihilation.s1_A ≈ 2.5
    @test params.annihilation.a_3p2 ≈ -0.8
end

@testset "The mesh p² is method-free" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r, h = GIModel.radial_grid(12, 3.0)
    # The convenience form taking parameters is the bare one. It used to dispatch
    # on a basis type parameter carried by GIParameters, with two methods whose
    # bodies were identical.
    @test GIModel.p2_operator(params, mq["c"], 1, r, h) ≈
          GIModel.p2_operator(mq["c"], 1, r, h)
end

@testset "OscillatorSolver channel solve returns native waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    sol = channel_solution(
        params,
        ConstituentMasses(mq["c"], mq["c"]),
        0;
        solver = OscillatorSolver(beta_grid = [0.7], converge = false),
        nlevels = 3,
    )
    @test sol isa ChannelRadialSolution{OscillatorWave}
    @test all(w -> w isa OscillatorWave && wave_norm(w) ≈ 1.0, sol.waves)

    # Sampling is an explicit view, not the stored representation.
    r, _ = GIModel.radial_grid(120, 14.0)
    vals, vecs, r = sampled_arrays(sol; grid = r)
    @test length(vals) == 3
    @test size(vecs) == (length(r), 3)
    @test vals[1] < vals[2] < vals[3]
    h = r[2] - r[1]
    for col in axes(vecs, 2)
        @test sum(abs2, vecs[:, col]) * h ≈ 1.0 rtol = 1e-10
    end

    @test_throws MethodError OscillatorSolver(ngrid = 75)
    @test_throws MethodError OscillatorSolver(rmax = 9.0)
end

@testset "reduced_mass" begin
    @test reduced_mass(ConstituentMasses(1.5, 0.3)) ≈ (1.5 * 0.3) / (1.5 + 0.3)
    @test reduced_mass(ConstituentMasses(2.0, 2.0)) ≈ 1.0
end

@testset "public package metadata and orbital labels" begin
    @test isfile(default_parameters_path())
    @test load_parameters(default_parameters_path()) isa GIParameters
    @test orbital_angular_momentum("S") == 0
    @test orbital_angular_momentum("G") == 4
    @test orbital_label(0) == "S"
    @test orbital_label(4) == "G"
    @test_throws ArgumentError orbital_angular_momentum("H")
    @test_throws ArgumentError orbital_label(5)
end

# Rewrite the active `central` method in the parameters TOML and reload.
function load_params_with_central(dir, central_name)
    p = joinpath(dir, "p.toml")
    s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
    @assert occursin("central = \"appendix_a_momentum_sandwich\"", s)
    write(
        p,
        replace(
            s,
            "central = \"appendix_a_momentum_sandwich\"" => "central = \"$central_name\"",
        ),
    )
    return load_parameters_and_quark_masses(p)
end

@testset "Strict runtime parameter files" begin
    original = GIModel.TOML.parsefile(default_parameters_path())
    @test GIModel.gi_parameters_from_raw(original) isa GIParameters
    for (section, key) in (("potential", "central"), ("relativistic_factors", "epsilon_c"), ("relativistic_factors", "contact_momentum_sandwich"), ("annihilation", "p1_A_np"))
        raw = deepcopy(original)
        delete!(raw[section], key)
        @test_throws ArgumentError GIModel.gi_parameters_from_raw(raw)
    end
    for (section, key, value) in (("potential", "Lambda_MeV", 200), ("potential", "b_GeV2", Inf), ("relativistic_factors", "epsilon_c", true), ("fine_structure", "enabled", 1))
        raw = deepcopy(original); raw[section][key] = value
        @test_throws ArgumentError GIModel.gi_parameters_from_raw(raw)
    end
    raw = deepcopy(original); delete!(raw, "annihilation")
    @test GIModel.gi_parameters_from_raw(raw).annihilation == AnnihilationAmplitudes()
    raw = deepcopy(original); raw["masses"]["m_c_MeV"] = -1
    @test_throws ArgumentError GIModel.quark_masses_from_raw(raw)
end

@testset "Parameter construction never supplies implicit physics" begin
    p = load_parameters(default_parameters_path())
    @test_throws UndefKeywordError GIParameters(; potential = p.potential, smearing = p.smearing)
    rebuilt = GIParameters(; (name => getfield(p, name) for name in fieldnames(typeof(p)))...)
    @test rebuilt == p
    @test GIParameters(p; potential = ConfinementPotential(p.potential; b = 0.2)).factors == p.factors
end
