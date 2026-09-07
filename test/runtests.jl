# Invoked by `Pkg.test()` with the package environment already active.
# To run manually: `julia --project=. test/runtests.jl` from the repo root.

using Test
using FiniteDifferences
using LinearAlgebra
using QuadGK
using SpecialFunctions: erf
using GIModel

root = dirname(@__DIR__)

# Explicit diagnostic/plot adapter. Production code consumes `radial_wave`
# directly; tests that validate sampled arrays request that view deliberately.
function sampled_arrays(sol::ChannelRadialSolution; grid = nothing)
    sampled = if all(w -> w isa MeshWave, sol.waves)
        sol.waves
    else
        isnothing(grid) && throw(ArgumentError("sampling an oscillator solution requires `grid`"))
        [sample_wave(w, grid) for w in sol.waves]
    end
    r = isnothing(grid) ? first(sampled).r : collect(Float64, grid)
    return sol.eigenvalues_GeV, hcat((w.u for w in sampled)...), r
end

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

@testset "central method `$name` code path (finite S-wave energy)" for (
    name,
    MethodType,
    ngrid,
    rmax,
) in [
    ("appendix_a_smearing", AppendixASmearing3D, 120, 12.0),
    ("coulomb_1d_smear", Coulomb1DSmearing, 120, 12.0),
    ("appendix_a_derivative_g", AppendixADerivativeG, 120, 12.0),
    ("appendix_a_closed_form", AppendixAClosedForm, 120, 12.0),
    ("appendix_a_momentum_sandwich", AppendixAMomentumSandwich, 80, 10.0),
]
    mktempdir() do d
        params, mq = load_params_with_central(d, name)
        @test params.central isa MethodType
        mc = mq["c"]
        sol = GIModel.channel_solution(
            params, ConstituentMasses(mc, mc), 0;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax, kinetic = :relativistic),
        )
        @test all(isfinite, sol.eigenvalues_GeV)
        @test sol.eigenvalues_GeV[1] < sol.eigenvalues_GeV[2]
    end
end

@testset "1D Coulomb smear: constant vector unchanged on uniform grid" begin
    h = 0.05
    n = 200
    r = collect(h:h:(h*n))
    c = GIModel.convolve_1d_gaussian_same_length(ones(n), r, h, 0.3)
    @test maximum(abs.(c .- 1.0)) < 1e-12
end

@testset "3D isotropic Gaussian smear: constant vector unchanged (with tail coverage)" begin
    # For a normalized 3D Gaussian ρ(|R-r|), the convolution should preserve constants.
    # On a finite radial mesh this is only true once the mesh extends far enough past the
    # region of interest to include the Gaussian tail. We mimic the production extension
    # used by `smeared_central_values`: Δr = 8/σ gives exp(-64) tail suppression.
    σ = 1.0
    h = 0.05
    rmax0 = 8.0
    r = collect(h:h:rmax0)
    n_tail = Int(ceil(8 / (σ * h)))
    r_ext = vcat(r, collect(range(rmax0 + h, rmax0 + n_tail * h; step = h)))
    out = GIModel.smear_3d_radial(ones(length(r_ext)), r_ext, σ)
    @test maximum(abs.(out[1:length(r)] .- 1.0)) < 2e-3

    # Also exercise the R→0 branch by explicitly including r=0 in the mesh.
    r0 = collect(0.0:h:rmax0)
    r0_ext = vcat(r0, collect(range(rmax0 + h, rmax0 + n_tail * h; step = h)))
    out0 = GIModel.smear_3d_radial(ones(length(r0_ext)), r0_ext, σ)
    @test abs(out0[1] - 1.0) < 2e-3
end

@testset "central_potential_values method dispatch" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mc = mq["c"]
    r, _h = GIModel.radial_grid(80, 8.0)
    pointwise =
        GIModel.central_potential_values(params, mc, mc, r; method = PointwiseCentral())
    @test pointwise ≈ [GIModel.central_potential(ri, params) for ri in r]
    @test GIModel.central_potential_values(params, mc, mc, r; method = Coulomb1DSmearing()) ≈
          GIModel.coulomb_1d_smeared_central_values(params, mc, mc, r)
    derivative = GIModel.central_potential_values(
        params,
        mc,
        mc,
        r;
        method = AppendixADerivativeG(),
    )
    @test all(isfinite, derivative)
    @test maximum(abs.(derivative .- pointwise)) > 0.0
    closed =
        GIModel.central_potential_values(params, mc, mc, r; method = AppendixAClosedForm())
    @test all(isfinite, closed)
    @test maximum(abs.(closed .- pointwise)) > 0.0
    sandwich_view = GIModel.central_potential_values(
        params,
        mc,
        mc,
        r;
        method = AppendixAMomentumSandwich(),
    )
    @test sandwich_view ≈ closed
    @test GIModel.smeared_coulomb_G_closed(params, mc, mc, 0.0) ≈
          GIModel.smeared_coulomb_G_closed(params, mc, mc, 1.0e-10)
    @test GIModel.smeared_confinement_S_closed(params, mc, mc, 0.0) ≈
          GIModel.smeared_confinement_S_closed(params, mc, mc, 1.0e-10)
    rcheck = 1.4
    δ = 1.0e-4
    fd_g_prime =
        (
            GIModel.smeared_coulomb_G_closed(params, mc, mc, rcheck + δ) -
            GIModel.smeared_coulomb_G_closed(params, mc, mc, rcheck - δ)
        ) / (2δ)
    fd_g_second =
        (
            GIModel.smeared_coulomb_G_closed(params, mc, mc, rcheck + δ) -
            2 * GIModel.smeared_coulomb_G_closed(params, mc, mc, rcheck) +
            GIModel.smeared_coulomb_G_closed(params, mc, mc, rcheck - δ)
        ) / δ^2
    fd_s_prime =
        (
            GIModel.smeared_confinement_S_closed(params, mc, mc, rcheck + δ) -
            GIModel.smeared_confinement_S_closed(params, mc, mc, rcheck - δ)
        ) / (2δ)
    @test GIModel.smeared_coulomb_G_prime_closed(params, mc, mc, rcheck) ≈ fd_g_prime rtol =
        1e-6
    @test GIModel.smeared_coulomb_G_second_closed(params, mc, mc, rcheck) ≈ fd_g_second rtol =
        1e-5
    @test GIModel.smeared_confinement_S_prime_closed(params, mc, mc, rcheck) ≈ fd_s_prime rtol =
        1e-7
    @test GIModel.tensor_kernel_smeared_coulomb(params, mc, mc, rcheck) ≈
          GIModel.smeared_coulomb_G_prime_closed(params, mc, mc, rcheck) / rcheck -
          GIModel.smeared_coulomb_G_second_closed(params, mc, mc, rcheck)
end

@testset "central_potential_method resolves names and rejects unknown ones" begin
    @test central_potential_method("pointwise") isa PointwiseCentral
    @test central_potential_method("coulomb_1d_smear") isa Coulomb1DSmearing
    @test central_potential_method("appendix_a_smearing") isa AppendixASmearing3D
    @test central_potential_method("appendix_a_derivative_g") isa AppendixADerivativeG
    @test central_potential_method("appendix_a_closed_form") isa AppendixAClosedForm
    @test central_potential_method("appendix_a_momentum_sandwich") isa
          AppendixAMomentumSandwich
    @test_throws ArgumentError central_potential_method("appendix_a_typo")
    @test GIModel.central_potential_path(AppendixASmearing3D()).name ==
          "experimental_3d_convl_a7a8"
    @test GIModel.central_potential_path(AppendixADerivativeG()).name ==
          "appendix_a_derivative_g"
    @test GIModel.central_potential_path(PointwiseCentral()).name == "pointwise_fd"
end

@testset "radial_laplacian_values constant is zero" begin
    h = 0.05
    n = 200
    r = collect(h:h:(h*n))
    lap = GIModel.radial_laplacian_values(ones(n), r)
    @test maximum(abs.(lap)) < 1e-12
end

@testset "Appendix A 3D smearing (constant preserves norm)" begin
    h = 0.02
    n = 2000
    r = collect(h:h:(h*n))
    σ = 1.0
    n_tail = σ > 0 ? max(0, Int(ceil(8 / (σ * h)))) : 0
    r_ext = n_tail > 0 ? vcat(r, collect((r[end]+h):h:(r[end]+n_tail*h))) : r
    w = GIModel.smear_3d_radial(ones(length(r_ext)), r_ext, σ)
    w = w[1:n]
    @test maximum(abs.(w .- 1.0)) < 0.01
    @test w[1] ≈ 1.0 atol = 0.01
    @test w[div(n, 2)] ≈ 1.0 atol = 0.01
end

@testset "closed-form Coulomb smearing matches QuadGK convolution" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    σ = GIModel.contact_smearing_sigma(params, m, m)

    function quadgk_smeared_coulomb(R)
        total = 0.0
        for (α, γ) in zip(GIModel.ALPHA_COEFFS, GIModel.ALPHA_GAMMAS)
            integrand(rp) = begin
                g =
                    rp == 0 ? -8 * α * γ / (3 * sqrt(π)) :
                    -4 * α * erf(γ * rp) / (3 * rp)
                pre = σ / (sqrt(π) * R)
                pre * rp * (exp(-(σ * (R - rp))^2) - exp(-(σ * (R + rp))^2)) * g
            end
            val, _err = quadgk(integrand, 0.0, Inf; rtol = 1e-9)
            total += val
        end
        total
    end

    for R in (0.2, 0.8, 2.0)
        @test GIModel.smeared_coulomb_G_closed(params, m, m, R) ≈ quadgk_smeared_coulomb(R) rtol =
            1e-10 atol = 1e-10
    end
end

@testset "L·S and spin_dot algebra" begin
    @test GIModel.LdotS(1, 1, 0) ≈ -2.0
    @test GIModel.LdotS(1, 1, 1) ≈ -1.0
    @test GIModel.LdotS(1, 1, 2) ≈ 1.0
    @test GIModel.spin_dot(1) ≈ -0.75
    @test GIModel.spin_dot(3) ≈ 0.25
end

@testset "triplet fine-structure angular factors" begin
    for L = 1:4
        js = collect((L-1):(L+1))
        weights = [2J + 1 for J in js]
        ldot = [GIModel.LdotS(L, 1, J) for J in js]
        tensor = [GIModel.tensor_triplet_LJ(L, J, 1) for J in js]
        @test sum(weights .* ldot) ≈ 0.0 atol = 1e-12
        @test sum(weights .* tensor) ≈ 0.0 atol = 1e-12
    end
    @test GIModel.tensor_triplet_LJ(1, 0, 1) ≈ -4.0
    @test GIModel.tensor_triplet_LJ(1, 1, 1) ≈ 2.0
    @test GIModel.tensor_triplet_LJ(1, 2, 1) ≈ -0.4
end

@testset "fine_structure_split: S-wave and P-wave triplet" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    sol_s = GIModel.channel_solution(
        params, ConstituentMasses(m, m), 0;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
    )
    wave_s = radial_wave(sol_s, 1)
    @test GIModel.fine_structure_split(
        params,
        ConstituentMasses(m, m),
        FineStructureMultiplet("S", 3, 1),
        wave_s,
    ) == 0.0
    @test GIModel.fine_structure_split(
        params,
        ConstituentMasses(m, m),
        FineStructureMultiplet("S", 1, 0),
        wave_s,
    ) == 0.0
    sol_p = GIModel.channel_solution(
        params, ConstituentMasses(m, m), 1;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
    )
    wave_p = radial_wave(sol_p, 1)
    δ0 = GIModel.fine_structure_split(
        params,
        ConstituentMasses(m, m),
        FineStructureMultiplet("P", 3, 0),
        wave_p,
    )
    δ1 = GIModel.fine_structure_split(
        params,
        ConstituentMasses(m, m),
        FineStructureMultiplet("P", 3, 1),
        wave_p,
    )
    δ2 = GIModel.fine_structure_split(
        params,
        ConstituentMasses(m, m),
        FineStructureMultiplet("P", 3, 2),
        wave_p,
    )
    @test isfinite(δ0) && isfinite(δ1) && isfinite(δ2)
    @test δ0 != δ1 || δ1 != δ2
    @test GIModel.fine_structure_split(
        params,
        ConstituentMasses(m, m),
        FineStructureMultiplet("P", 1, 0),
        wave_p,
    ) == 0.0
end

@testset "fine_structure_components: decomposition sums correctly" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    sol_p = GIModel.channel_solution(
        params, ConstituentMasses(m, m), 1;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
    )
    wave_p = radial_wave(sol_p, 1)
    for J in (0, 1, 2)
        comp = GIModel.fine_structure_components(
            params,
            ConstituentMasses(m, m),
            FineStructureMultiplet("P", 3, J),
            wave_p,
        )
        @test comp.spin_orbit ≈ comp.spin_orbit_vector + comp.spin_orbit_thomas atol = 1e-12
        @test comp.total ≈ comp.spin_orbit + comp.tensor atol = 1e-12
        @test comp.total ≈ GIModel.fine_structure_split(
            params,
            ConstituentMasses(m, m),
            FineStructureMultiplet("P", 3, J),
            wave_p,
        ) atol = 1e-12
    end
end

@testset "same-J spin-orbit mixing diagnostics" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mc = mq["c"]
    mb = mq["b"]

    generic_block = GIModel.MixingBlock(
        "three-state smoke test",
        [
            GIModel.BasisState(1, "S", 3, 1),
            GIModel.BasisState(1, "D", 3, 1),
            GIModel.BasisState(2, "S", 3, 1),
        ],
        [3.0 0.01 0.0; 0.01 3.2 0.02; 0.0 0.02 3.6];
        mechanism = "test_mixing",
    )
    generic = GIModel.diagonalize_mixing_block(generic_block)
    @test length(generic.masses) == 3
    @test issorted(generic.masses)
    @test generic.vectors' * generic.vectors ≈ [i == j ? 1.0 : 0.0 for i = 1:3, j = 1:3] atol =
        1e-12
    @test_throws ArgumentError GIModel.BasisState(0, "S", 3, 1)
    @test_throws ArgumentError GIModel.BasisState(1, "S", 2, 1)
    @test_throws ArgumentError GIModel.BasisState(1, "S", 3, -1)
    @test_throws ArgumentError GIModel.MixingBlock(
        "empty block", GIModel.BasisState[], zeros(0, 0),
    )
    @test_throws ArgumentError GIModel.MixingBlock(
        "bad asymmetric block",
        [GIModel.BasisState(1, "S", 3, 1), GIModel.BasisState(1, "D", 3, 1)],
        [1.0 0.2; 0.1 2.0],
    )
    @test_throws ArgumentError GIModel.MixingBlock(
        "duplicate identities",
        [
            GIModel.BasisState(1, "S", 3, 1; label = "first"),
            GIModel.BasisState(1, "S", 3, 1; label = "second"),
        ],
        Matrix{Float64}(I, 2, 2),
    )
    @test_throws ArgumentError GIModel.MixingResult(
        generic_block, reverse(generic.masses), generic.vectors,
    )
    @test_throws ArgumentError GIModel.MixingResult(
        generic_block, generic.masses, 2 .* generic.vectors,
    )

    sol_cc = GIModel.channel_solution(
        params, ConstituentMasses(mc, mc), 1;
        nlevels = 1,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
    )
    radial_cc = radial_wave(sol_cc, 1)
    off_cc = GIModel.spin_orbit_mixing_components(
        params,
        ConstituentMasses(mc, mc),
        "P",
        radial_cc,
    )
    @test off_cc.total == 0.0
    mix_cc = GIModel.same_j_mixing(3.5, 3.6, off_cc.total)
    @test mix_cc.masses ≈ [3.5, 3.6]
    @test mix_cc.theta_deg ≈ 0.0 atol = 1e-12
    @test mix_cc.block isa GIModel.MixingBlock
    @test mix_cc.block.mechanism == "antisymmetric_spin_orbit"

    sol_bc = GIModel.channel_solution(
        params, ConstituentMasses(mb, mc), 1;
        nlevels = 1,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 24.0),
    )
    radial_bc = radial_wave(sol_bc, 1)
    off_bc = GIModel.spin_orbit_mixing_components(
        params,
        ConstituentMasses(mb, mc),
        "P",
        radial_bc,
    )
    @test isfinite(off_bc.total)
    @test off_bc.total != 0.0

    triplet_shift = GIModel.fine_structure_split(
        params,
        ConstituentMasses(mb, mc),
        FineStructureMultiplet("P", 3, 1),
        radial_bc,
    )
    central_bc = sol_bc.eigenvalues_GeV[1]
    mix_bc = GIModel.same_j_mixing(central_bc, central_bc + triplet_shift, off_bc.total)
    @test isfinite(mix_bc.theta_deg)
    @test minimum(mix_bc.masses) < central_bc < maximum(mix_bc.masses)
end

@testset "legacy fine-structure ε factors are scalar (1+ε) multipliers" begin
    mktempdir() do d
        p0 = joinpath(d, "p0.toml")
        pt = joinpath(d, "pt.toml")
        pv = joinpath(d, "pv.toml")
        ps = joinpath(d, "ps.toml")

        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s0 = replace(s, r"^epsilon_c\s*=.*$"m => "epsilon_c = 0.0")
        s0 = replace(s0, r"^epsilon_t\s*=.*$"m => "epsilon_t = 0.0")
        s0 = replace(s0, r"^epsilon_so_vector\s*=.*$"m => "epsilon_so_vector = 0.0")
        s0 = replace(s0, r"^epsilon_so_scalar\s*=.*$"m => "epsilon_so_scalar = 0.0")
        s0 = replace(
            s0,
            "fine_structure_momentum_sandwich = true" => "fine_structure_momentum_sandwich = false",
        )
        write(p0, s0)

        write(pt, replace(s0, r"^epsilon_t\s*=.*$"m => "epsilon_t = 0.5"))
        write(pv, replace(s0, r"^epsilon_so_vector\s*=.*$"m => "epsilon_so_vector = 0.5"))
        write(ps, replace(s0, r"^epsilon_so_scalar\s*=.*$"m => "epsilon_so_scalar = 0.5"))

        params0, mq0 = load_parameters_and_quark_masses(p0)
        mc = mq0["c"]
        _v_p, umat_p, r_p = sampled_arrays(GIModel.channel_solution(
            params0, ConstituentMasses(mc, mc), 1;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
        ))
        h_p = r_p[2] - r_p[1]
        u1p = collect(umat_p[:, 1])

        # Use J=2 to avoid any accidental L·S or tensor zeros.
        comp0 = GIModel.fine_structure_components(
            params0,
            mc,
            mc,
            "P",
            3,
            2,
            u1p,
            r_p,
            h_p,
        )
        @test comp0.tensor != 0.0
        @test comp0.spin_orbit_vector != 0.0
        @test comp0.spin_orbit_thomas != 0.0

        compt = GIModel.fine_structure_components(
            load_parameters(pt),
            mc,
            mc,
            "P",
            3,
            2,
            u1p,
            r_p,
            h_p,
        )
        compv = GIModel.fine_structure_components(
            load_parameters(pv),
            mc,
            mc,
            "P",
            3,
            2,
            u1p,
            r_p,
            h_p,
        )
        comps = GIModel.fine_structure_components(
            load_parameters(ps),
            mc,
            mc,
            "P",
            3,
            2,
            u1p,
            r_p,
            h_p,
        )

        @test compt.tensor / comp0.tensor ≈ 1.5 rtol = 1e-12 atol = 0.0
        @test compv.spin_orbit_vector / comp0.spin_orbit_vector ≈ 1.5 rtol = 1e-12 atol =
            0.0
        @test comps.spin_orbit_thomas / comp0.spin_orbit_thomas ≈ 1.5 rtol = 1e-12 atol =
            0.0
    end
end

@testset "tensor proxy keeps Coulomb color factor" begin
    mktempdir() do d
        p0 = joinpath(d, "p0.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s0 = replace(s, r"^epsilon_t\s*=.*$"m => "epsilon_t = 0.0")
        s0 = replace(
            s0,
            "fine_structure_momentum_sandwich = true" => "fine_structure_momentum_sandwich = false",
        )
        s0 = replace(
            s0,
            "fine_structure_smeared_kernels = true" => "fine_structure_smeared_kernels = false",
        )
        write(p0, s0)

        params0, mq0 = load_parameters_and_quark_masses(p0)
        mc = mq0["c"]
        _v_p, umat_p, r_p = sampled_arrays(GIModel.channel_solution(
            params0, ConstituentMasses(mc, mc), 1;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
        ))
        h_p = r_p[2] - r_p[1]
        u1p = collect(umat_p[:, 1])

        Itk = GIModel.radial_expect_udr(
            u1p,
            r_p,
            h_p,
            (ri, i) -> GIModel.tensor_kernel_coulomb_running(ri),
        )
        comp = GIModel.fine_structure_components(
            params0,
            mc,
            mc,
            "P",
            3,
            2,
            u1p,
            r_p,
            h_p,
        )
        @test comp.I_tk ≈ Itk rtol = 1e-12 atol = 0.0
        expected = Itk * GIModel.tensor_triplet_LJ(1, 2, 1) / (12.0 * mc * mc)
        @test comp.tensor ≈ expected rtol = 1e-12 atol = 0.0
    end
end

@testset "tensor off-diagonal same-J angular factor" begin
    @test GIModel.tensor_triplet_offdiag_sameJ(0, 1) ≈ 0.0
    @test GIModel.tensor_triplet_offdiag_sameJ(1, 1) ≈ sqrt(8.0)
    @test GIModel.tensor_triplet_offdiag_sameJ(2, 1) ≈ 6sqrt(6.0) / 5
    @test GIModel.tensor_triplet_offdiag_sameJ(1, 0) ≈ 0.0
end

@testset "contact_hyperfine_shift: normalization + spin algebra" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    _vals, umat, r = sampled_arrays(GIModel.channel_solution(
        params, ConstituentMasses(m, m), 0;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
    ))
    u1s = collect(umat[:, 1])

    δ_triplet = GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u1s, r)
    δ_singlet = GIModel.contact_hyperfine_shift(params, m, m, "S", 1, u1s, r)
    @test isfinite(δ_triplet) && isfinite(δ_singlet)
    # spin_dot(3) / spin_dot(1) = 0.25 / (-0.75) = -1/3
    @test δ_triplet ≈ (-1 / 3) * δ_singlet rtol = 1e-12

    # contact_hyperfine_shift is defined as an expectation over u(r) with an
    # internal physical normalization ∫|u|^2 dr = 1, so scaling u must not
    # change the shift.
    δ_scaled = GIModel.contact_hyperfine_shift(params, m, m, "S", 3, 7.0 .* u1s, r)
    @test δ_scaled ≈ δ_triplet rtol = 1e-12
end

@testset "fine structure uses u(r) normalization (scale invariant)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    _v_p, umat_p, r_p = sampled_arrays(GIModel.channel_solution(
        params, ConstituentMasses(m, m), 1;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = 200, rmax = 20.0),
    ))
    h_p = r_p[2] - r_p[1]
    u1p = collect(umat_p[:, 1])
    for J in (0, 1, 2)
        base = GIModel.fine_structure_components(
            params,
            m,
            m,
            "P",
            3,
            J,
            u1p,
            r_p,
            h_p,
        )
        scaled_hi = GIModel.fine_structure_components(
            params,
            m,
            m,
            "P",
            3,
            J,
            3.0 .* u1p,
            r_p,
            h_p,
        )
        scaled_lo = GIModel.fine_structure_components(
            params,
            m,
            m,
            "P",
            3,
            J,
            0.2 .* u1p,
            r_p,
            h_p,
        )
        @test scaled_hi.total ≈ base.total rtol = 1e-12 atol = 0.0
        @test scaled_lo.total ≈ base.total rtol = 1e-12 atol = 0.0
        @test scaled_hi.spin_orbit_vector ≈ base.spin_orbit_vector rtol = 1e-12 atol = 0.0
        @test scaled_lo.spin_orbit_vector ≈ base.spin_orbit_vector rtol = 1e-12 atol = 0.0
        @test scaled_hi.spin_orbit_thomas ≈ base.spin_orbit_thomas rtol = 1e-12 atol = 0.0
        @test scaled_lo.spin_orbit_thomas ≈ base.spin_orbit_thomas rtol = 1e-12 atol = 0.0
        @test scaled_hi.spin_orbit ≈ base.spin_orbit rtol = 1e-12 atol = 0.0
        @test scaled_lo.spin_orbit ≈ base.spin_orbit rtol = 1e-12 atol = 0.0
        @test scaled_hi.tensor ≈ base.tensor rtol = 1e-12 atol = 0.0
        @test scaled_lo.tensor ≈ base.tensor rtol = 1e-12 atol = 0.0
    end
end

@testset "GI erf profile basic symmetries" begin
    @test erf(0.0) ≈ 0.0 atol = 1e-7
    for x in (0.05, 0.3, 0.8, 1.5)
        @test erf(x) + erf(-x) ≈ 0.0 atol = 1e-12
    end
end

@testset "GI momentum-space alpha_s profile" begin
    @test GIModel.alpha_s_q(0.0) ≈ sum(GIModel.ALPHA_COEFFS)
    @test GIModel.alpha_s_q(1.0) ≈
          sum(a * exp(-1.0^2 / (4g^2)) for (a, g) in zip(GIModel.ALPHA_COEFFS, GIModel.ALPHA_GAMMAS))
    @test GIModel.alpha_s_q(3.0) < GIModel.alpha_s_q(1.0) < GIModel.alpha_s_q(0.0)
end

@testset "GI erf profile second derivative consistency" begin
    for x in (0.1, 0.4, 1.1, 1.8)
        δ = 1e-6 * max(1.0, x)
        num = (GIModel.erf_prime(x + δ) - GIModel.erf_prime(x - δ)) / (2δ)
        ana = GIModel.erf_second(x)
        @test ana ≈ num rtol = 2e-5 atol = 1e-9
        @test GIModel.erf_second(-x) ≈ -ana rtol = 1e-12 atol = 1e-12
    end
end

@testset "Coulomb derivative consistency (GI erf profile)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    for r0 in (0.05, 0.2, 1.3)
        δ = 1e-6 * max(1.0, r0)
        num = (GIModel.alpha_s_r(r0 + δ) - GIModel.alpha_s_r(r0 - δ)) / (2δ)
        @test GIModel.alpha_s_prime_r(r0) ≈ num rtol = 1e-6 atol = 1e-10
    end

    for r0 in (0.05, 0.2, 1.3)
        δ = 1e-6 * max(1.0, r0)
        num = (GIModel.alpha_s_prime_r(r0 + δ) - GIModel.alpha_s_prime_r(r0 - δ)) / (2δ)
        @test GIModel.alpha_s_second_r(r0) ≈ num rtol = 2e-5 atol = 1e-10
    end

    V(r) = -(4 / 3) * GIModel.alpha_s_r(r) / r
    for r0 in (0.2, 1.0)
        δ = 1e-6 * max(1.0, r0)
        num = (V(r0 + δ) - V(r0 - δ)) / (2δ)
        @test GIModel.dV_coul_central_dr(r0, params) ≈ num rtol = 1e-6 atol = 1e-10
    end
end

@testset "analytic radial derivatives match FiniteDifferences" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    d1 = central_fdm(5, 1)
    d2 = central_fdm(5, 2)

    for r0 in (0.08, 0.4, 1.4)
        @test GIModel.alpha_s_prime_r(r0) ≈ d1(GIModel.alpha_s_r, r0) rtol = 1e-10 atol =
            1e-10
        @test GIModel.alpha_s_second_r(r0) ≈ d2(GIModel.alpha_s_r, r0) rtol = 1e-9 atol =
            1e-9
    end

    Gs(r) = GIModel.smeared_coulomb_G_closed(params, m, m, r)
    for r0 in (0.2, 0.8, 2.0)
        @test GIModel.smeared_coulomb_G_prime_closed(params, m, m, r0) ≈ d1(Gs, r0) rtol =
            1e-9 atol = 1e-9
        @test GIModel.smeared_coulomb_G_second_closed(params, m, m, r0) ≈ d2(Gs, r0) rtol =
            1e-7 atol = 1e-8
    end
end

@testset "spin–orbit convention: Eq. (6) uses α_s/r^3 (no α_s')" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s0 = replace(
            s,
            "fine_structure_momentum_sandwich = true" => "fine_structure_momentum_sandwich = false",
        )
        s0 = replace(
            s0,
            "fine_structure_smeared_kernels = true" => "fine_structure_smeared_kernels = false",
        )
        write(p, s0)
        params, mq = load_parameters_and_quark_masses(p)
        h = 0.02
        r = collect(0.10:h:4.00)
        u = exp.(-r)

        I_cm = GIModel.radial_expect_udr(
            u,
            r,
            h,
            (ri, i) -> begin
                r0 = max(ri, 1.0e-8)
                (4.0 / 3.0) * GIModel.alpha_s_r(r0) / r0^3
            end,
        )
        I_deriv = GIModel.radial_expect_udr(
            u,
            r,
            h,
            (ri, i) -> begin
                r0 = max(ri, 1.0e-8)
                (1.0 / r0) * GIModel.dV_coul_central_dr(r0, params)
            end,
        )
        # With running α_s(r)=∑ α_k erf(γ_k r), α_s'(r)>0 so these must differ:
        #   (4/3) α_s(r)/r^3  - (1/r) d/dr[-4α_s(r)/(3r)] = (4/3) α_s'(r)/r^2  > 0.
        @test I_cm > I_deriv
        @test abs(I_cm - I_deriv) > 1e-6

        # Guardrail: fine_structure_components must use I_cm (Eq. (6)), not I_deriv.
        m = mq["c"]
        comp = GIModel.fine_structure_components(
            params,
            m,
            m,
            "P",
            3,
            2,
            collect(u),
            collect(r),
            h,
        )
        scaled_cm = (1 + params.factors.epsilon_so_vector) * I_cm
        @test comp.I_vector_11 ≈ scaled_cm rtol = 1e-12 atol = 0.0
        @test comp.I_vector_22 ≈ scaled_cm rtol = 1e-12 atol = 0.0
        @test comp.I_vector_12 ≈ scaled_cm rtol = 1e-12 atol = 0.0
        ls = GIModel.LdotS(1, 1, 2)
        expected_vec = ls * (3 / (2m^2)) * scaled_cm
        @test comp.spin_orbit_vector ≈ expected_vec rtol = 1e-12 atol = 0.0

        I_scalar = GIModel.radial_expect_udr(
            u,
            r,
            h,
            (ri, i) -> begin
                r0 = max(ri, 1.0e-8)
                params.potential.b / r0
            end,
        )
        scaled_scalar = (1 + params.factors.epsilon_so_scalar) * I_scalar
        @test comp.I_scalar_11 ≈ scaled_scalar rtol = 1e-12 atol = 0.0
        @test comp.I_scalar_22 ≈ scaled_scalar rtol = 1e-12 atol = 0.0
        expected_tp = -ls * scaled_scalar / (2m^2)
        @test comp.spin_orbit_thomas ≈ expected_tp rtol = 1e-12 atol = 0.0
    end
end

@testset "OGE tensor kernel matches finite-difference derivatives of G(r)" begin
    # Tensor kernel is built from the Coulomb piece G(r) = -4 α_s(r)/(3 r) as
    #   K(r) = (1/r) dG/dr - d²G/dr².
    # This regression test guards signs/factors in the running-α_s implementation.
    G(r) = GIModel.coulomb_G_running(r)
    for r0 in (0.08, 0.15, 0.4, 1.0, 2.0)
        δ = 1e-5 * max(1.0, r0)
        d1 = (G(r0 + δ) - G(r0 - δ)) / (2δ)
        d2 = (G(r0 + δ) - 2G(r0) + G(r0 - δ)) / (δ^2)
        num = d1 / r0 - d2
        ana = GIModel.tensor_kernel_coulomb_running(r0)
        @test ana ≈ num rtol = 2e-4 atol = 1e-7
    end
end

@testset "reduced radial u(r) expectation normalization" begin
    h = 0.1
    r = collect(h:h:(3h))
    u = [1.0, 2.0, 3.0]
    @test GIModel.radial_expect_udr(u, r, h, (ri, i) -> 1.0) ≈ 1.0 atol = 1e-12

    manual_weights = abs2.(u) ./ sum(abs2.(u))
    manual_mean_r = sum(manual_weights .* r)
    @test GIModel.radial_expect_udr(u, r, h, (ri, i) -> ri) ≈ manual_mean_r atol = 1e-12
end

@testset "reduced radial u(r) expectations validate mesh alignment" begin
    r = [0.1, 0.2, 0.3]
    u = [1.0, 2.0]
    @test_throws ArgumentError GIModel.physical_u_norm(r, 0.1, u)
    @test_throws ArgumentError GIModel.radial_expect_udr(u, r, 0.1, (ri, i) -> 1.0)
    @test_throws ArgumentError GIModel.physical_u_norm(r[1:2], -0.1, u)
    @test_throws ArgumentError GIModel.physical_u_norm(r[1:2], NaN, u)

    # Guardrail: expectation values assume a uniform r mesh, and the caller-supplied
    # spacing `h` must match the actual r[i]-r[i-1].
    u3 = [1.0, 1.0, 1.0]
    r_nonuniform = [0.1, 0.2, 0.31]
    @test_throws ArgumentError GIModel.physical_u_norm(r_nonuniform, 0.1, u3)
    @test_throws ArgumentError GIModel.radial_expect_udr(
        u3,
        r_nonuniform,
        0.1,
        (ri, i) -> 1.0,
    )
    @test_throws ArgumentError GIModel.physical_u_norm(r, 0.11, u3)
end

@testset "contact hyperfine: only S-waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    r = collect(0.05:0.05:1.0)
    u = exp.(-2.0 .* r)
    @test GIModel.contact_hyperfine_shift(params, m, m, "P", 3, u, r) == 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "D", 3, u, r) == 0.0
end

@testset "contact smearing σ implements Appendix A (A9)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m1, m2 = mq["c"], mq["b"]

    σ = GIModel.contact_smearing_sigma(params, m1, m2)
    σ_swapped = GIModel.contact_smearing_sigma(params, m2, m1)
    @test σ ≈ σ_swapped rtol = 0.0 atol = 0.0

    mass_factor = 4 * m1 * m2 / (m1 + m2)^2
    reduced_twice = 2 * m1 * m2 / (m1 + m2)
    σ_manual = sqrt(
        params.smearing.sigma0^2 * (0.5 + 0.5 * mass_factor^4) +
        params.smearing.s^2 * reduced_twice^2,
    )
    @test σ ≈ σ_manual rtol = 0.0 atol = 0.0

    m = mq["c"]
    σ_equal = GIModel.contact_smearing_sigma(params, m, m)
    σ_equal_manual = sqrt(params.smearing.sigma0^2 + params.smearing.s^2 * m^2)
    @test σ_equal ≈ σ_equal_manual rtol = 0.0 atol = 0.0
end

@testset "contact hyperfine shift uses u(r) normalization" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    r = collect(0.05:0.05:1.0)
    u = exp.(-2.0 .* r)
    base = GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u, r)
    @test base != 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, 3.0 .* u, r) ≈ base rtol =
        1e-12 atol = 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, 0.2 .* u, r) ≈ base rtol =
        1e-12 atol = 0.0
end

@testset "contact hyperfine matches radial_expect_udr convention" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    r = collect(0.05:0.05:1.0)
    h = r[2] - r[1]
    u = exp.(-2.0 .* r)
    masses = ConstituentMasses(m, m)
    expectation = GIModel.radial_expect_udr(
        u,
        r,
        h,
        (ri, i) -> GIModel.smeared_contact_kernel(params, masses, ri),
    )
    manual =
        (1.0 + params.factors.epsilon_c) *
        (32 * pi / (9 * m * m)) *
        expectation *
        GIModel.spin_dot(3)
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u, r) ≈ manual rtol = 1e-12 atol =
        0.0
end

@testset "A15 contact is the Laplacian of the smeared Coulomb kernel" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["q"], mq["s"])
    sigma = GIModel.contact_smearing_sigma(params, masses)

    for r in (0.0, 0.15, 0.8, 2.0)
        expected = sum(zip(GIModel.ALPHA_COEFFS, GIModel.ALPHA_GAMMAS)) do (alpha, gamma)
            tau = inv(sqrt(inv(sigma^2) + inv(gamma^2)))
            alpha * GIModel.delta_sigma_3d(r, tau)
        end
        @test GIModel.smeared_contact_kernel(params, masses, r) ≈ expected rtol = 2e-15
    end

    @test fieldnames(FineStructure) == (:enabled,)
    @test_throws MethodError FineStructure(; k_tensor = 0.42)
    parameter_text = read(joinpath(root, "data", "parameters.provisional.toml"), String)
    @test !occursin("k_spin_orbit", parameter_text)
    @test !occursin("k_tensor", parameter_text)
end

@testset "post-A14 spin-dependent momentum exponent gives energy denominators" begin
    @test GIModel.gi_spin_dependent_side_exponent(0.0) ≈ 0.5
    @test GIModel.gi_spin_dependent_side_exponent(-0.168) ≈ 0.332

    p2_fact = eigen(Diagonal([0.0, 3.0]))
    B = GIModel.momentum_relativization_matrix(0.22, 0.22, GIModel.gi_spin_dependent_side_exponent(0.0), p2_fact)
    λ = 3.0
    E = sqrt(λ + 0.22^2)
    @test B[2, 2]^2 / 0.22^2 ≈ 1 / E^2 rtol = 1e-12
end

@testset "RadialChannelKey collapses nearly-equal masses" begin
    a = RadialChannelKey(1.628, 1.628, "S")
    nearby = 1.6280000000004
    b = RadialChannelKey(nearby, nearby, "S")
    @test a == b
    @test hash(a) == hash(b)
    @test RadialChannelKey(1.6, 1.6, "S") != RadialChannelKey(1.6, 1.6, "P")
    masses = GIModel.ConstituentMasses(1.6, 1.6)
    @test RadialChannelKey(masses, "S") == RadialChannelKey(1.6, 1.6, "S")
    exact = ConstituentMasses(nearby, nearby)
    @test exact.m1_GeV == nearby
    @test exact.m1_GeV != 1.628
end

@testset "ChannelRadialSolution retains its native radial wave" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    sol = GIModel.channel_solution(
            params, ConstituentMasses(m, m), 1;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = 120, rmax = 16.0),
        )
    wave = GIModel.radial_wave(sol, 1)
    @test wave === sol.waves[1]
    @test wave isa MeshWave
    @test wave.h ≈ wave.r[2] - wave.r[1]
    mult = GIModel.FineStructureMultiplet("P", 3, 2)
    masses = GIModel.ConstituentMasses(m, m)
    c1 = GIModel.fine_structure_components(
        params,
        masses,
        mult,
        wave,
    )
    split = GIModel.fine_structure_split(
        params, masses, mult, wave,
    )
    @test c1.total ≈ split rtol = 1e-12 atol = 0.0
end

@testset "mixing mechanisms identify model mass-matrix stages" begin
    @test GIModel.AntisymmetricSpinOrbit() isa GIModel.MixingMechanism
    @test GIModel.TensorMixing() isa GIModel.MixingMechanism
    @test GIModel.IsoscalarAnnihilation() isa GIModel.MixingMechanism
end

@testset "Table V strong-decay model (light 1S+1P)" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)
    # the two fit rows are exact by construction
    @test strong_decay_amplitude(model, sqrt(4 / 3), :A, 1, q_rho) ≈ 12.4 atol = 1e-9
    @test strong_decay_amplitude(model, -sqrt(2 / 9), :S, 0, q_B) ≈ -11.0 atol = 1e-9
    # predictions against the paper's numeric column (2 significant figures)
    q_kstar = decay_momentum(0.8921, 0.4957, 0.138)
    @test isapprox(strong_decay_amplitude(model, 1.0, :A, 1, q_kstar), 7.9; rtol = 0.05)
    q_a2 = decay_momentum(1.318, 0.5488, 0.138)
    @test isapprox(strong_decay_amplitude(model, sqrt(1 / 30), :A, 2, q_a2), 4.5; rtol = 0.05)
    q_f = decay_momentum(1.273, 0.138, 0.138)
    @test isapprox(strong_decay_amplitude(model, -sqrt(1 / 10), :A, 2, q_f), -11.0; rtol = 0.05)
    # below threshold: zero momentum and zero amplitude
    @test decay_momentum(1.0, 0.6, 0.6) == 0.0
    @test strong_decay_amplitude(model, 1.0, :A, 1, 0.0) == 0.0
    # reduced-amplitude classes
    @test reduced_decay_amplitude(model, :A0, 1.0) == model.A
    @test reduced_decay_amplitude(model, :Aprime, 2.0) == model.A
    @test reduced_decay_amplitude(model, :S, 1.0) ≈ model.S0 - 0.5 * model.A
    @test reduced_decay_amplitude(model, :D, 1.0) ≈ model.S0 - 0.3 * model.A
    @test reduced_decay_amplitude(model, :P, 1.0) ≈ model.S0 - 0.75 * model.A
    @test_throws ArgumentError reduced_decay_amplitude(model, :bogus, 1.0)
end

@testset "Table V charmed decays (A_c/S_c, footnote d)" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)

    # A_c is structure-independent, so it is the light A unchanged.
    @test reduced_decay_amplitude(model, :A_c, 3.0) == model.A

    # The charm mass ratio comes from the parameters TOML, not from constants
    # frozen inside src/ (which used to shadow it).
    _, mq_charm = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m_c, m_d = mq_charm["c"], mq_charm["d"]
    r = m_c / (m_c + m_d)
    @test r ≈ 1.628 / (1.628 + 0.220)   # unfreezing changed no value

    # Table IV (paper/vision_ocr/pages/page-013.md:39,55) prints
    #   S   = [3h - (1/2)      (g + h/4) q^2/beta^2  ] beta
    #   S_c = [3h - m_c/(m_d+m_c) A_c/beta q^2/beta_c^2] beta_c
    # so the S-family polynomial coefficient is the heavy fraction r, and the
    # light 1/2 is simply r at equal constituent masses. S_c must therefore be
    # MORE suppressed than S (r = 0.881 > 1/2); it is not the light formula.
    @test reduced_decay_amplitude(model, :S_c, 1.0; heavy_fraction = 0.5) ≈
          reduced_decay_amplitude(model, :S, 1.0)
    @test reduced_decay_amplitude(model, :S_c, 1.0; heavy_fraction = r) ≈
          model.S0 - r * model.A
    @test reduced_decay_amplitude(model, :S_c, 1.0; heavy_fraction = r) <
          reduced_decay_amplitude(model, :S, 1.0)

    # The charmed form factor exp[-(1/4)(m_c/(m_c+m_d))^2 q^2/beta^2] replaces
    # the light exp(-q^2/16 beta^2); at fixed q, coefficient, and class the
    # charm amplitude equals the light one rescaled by the Gaussian ratio.
    let q = 0.4, c = 0.9, L = 2, beta = model.beta_GeV
        light = strong_decay_amplitude(model, c, :A, L, q)
        # rebuild the light amplitude under the :A_c class = same reduced value,
        # so the only difference is the form factor and (optionally) recoil.
        charm = strong_decay_amplitude(model, c, :A_c, L, q; heavy_fraction = r)
        ratio = exp(-0.25 * r^2 * q^2 / beta^2) / exp(-q^2 / (16 * beta^2))
        @test isapprox(charm / light, ratio; rtol = 1e-10)
    end

    # The recoil multiplier m_c*beta/((m_c+m_d)*beta_c) = m_c/(m_c+m_d) with
    # beta_c=beta; A_c P-wave rows carry it, S-wave/1^3S_1 rows do not.
    let q = 0.5, c = -sqrt(1 / 5), L = 2
        no_rec = strong_decay_amplitude(model, c, :A_c, L, q; heavy_fraction = r, recoil = false)
        with_rec = strong_decay_amplitude(model, c, :A_c, L, q; heavy_fraction = r, recoil = true)
        # The multiplier itself is exactly r when beta_c = beta, but this ratio
        # is a quotient of two products, so it carries its own rounding and is
        # not bit-exact (it lands 1 ulp off r).
        @test isapprox(with_rec / no_rec, r; rtol = 1e-10)
    end

    # below threshold -> zero
    @test strong_decay_amplitude(model, 1.0, :A_c, 1, 0.0; heavy_fraction = r) == 0.0

    # Clean 1^3S_1 D* -> D pi rows reproduce the paper's column with NO refit
    # (masses: D*+=2.010, D0=1.865, pi+=0.1396; D*0=2.007, pi0=0.135).
    q_dstarp = decay_momentum(2.010, 1.865, 0.1396)
    @test isapprox(strong_decay_amplitude(model, -sqrt(2 / 3), :A_c, 1, q_dstarp; heavy_fraction = r),
        -0.34; atol = 0.05)
    q_dstar0 = decay_momentum(2.007, 1.865, 0.135)
    @test isapprox(strong_decay_amplitude(model, -sqrt(1 / 3), :A_c, 1, q_dstar0; heavy_fraction = r),
        -0.27; atol = 0.05)
end

@testset "Row-oriented decay API (DecayChannel decomposition)" begin
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)          # :table_iv default
    model_lead = calibrate_strong_decay_model(q_rho, q_B; convention = :leading)

    # spatial_overlap is the pure [DERIVED] SHO momentum factor.
    @test spatial_overlap(0.0, 1, 0.40) == 0.0
    @test spatial_overlap(0.359, 0, 0.40) ≈ GIModel._suppressed_factor(0.359, 0.40)

    # The 3-way decomposition multiplies back to the scalar amplitude.
    ch_rho = DecayChannel("rho", "pi", "pi", sqrt(4 / 3), :A, 1)
    a = decay_amplitude(model, ch_rho, q_rho; convention = :table_iv)
    @test a.total ≈ a.coefficient * a.reduced * a.spatial_overlap
    @test a.total ≈ strong_decay_amplitude(model, sqrt(4 / 3), :A, 1, q_rho)
    @test a.total ≈ 12.4 atol = 1e-9
    @test matrix_element(a) ≈ a.coefficient * a.reduced
    @test decay_width(a) ≈ a.total^2

    # leading vs table_iv: identical for structure-independent A, differ for D.
    ch_D = DecayChannel("rho2", "omega", "pi", -sqrt(1 / 36), :D, 1)
    q = 0.66
    @test decay_amplitude(model, ch_rho, q; convention = :leading).total ≈
          decay_amplitude(model, ch_rho, q; convention = :table_iv).total
    @test decay_amplitude(model_lead, ch_D, q; convention = :leading).reduced ≈ model_lead.S0
    @test decay_amplitude(model, ch_D, q; convention = :table_iv).reduced <
          decay_amplitude(model, ch_D, q; convention = :leading).reduced

    # MesonMasses resolves the momentum internally.
    masses = MesonMasses(Dict("rho" => 0.769, "pi" => 0.138))
    @test meson_mass(masses, "rho") == 0.769
    @test masses["pi"] == 0.138
    @test decay_amplitude(model, ch_rho, masses; convention = :table_iv).total ≈ 12.4 atol = 1e-9
    @test_throws Exception meson_mass(masses, "unregistered")

    # A charmed channel now carries its OWN mass ratio r = m_c/(m_c+m_d) rather
    # than having it inferred from the class; the A_c P-wave recoil multiplier is
    # still keyed on the class.
    _, mq_row = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r_c = mq_row["c"] / (mq_row["c"] + mq_row["d"])
    ch_c = DecayChannel("Kstar_c", "D", "pi", -sqrt(1 / 5), :A_c, 2; heavy_fraction = r_c)
    q_c = 0.4
    got = decay_amplitude(model, ch_c, q_c; convention = :leading).total
    want = strong_decay_amplitude(model, -sqrt(1 / 5), :A_c, 2, q_c;
        heavy_fraction = r_c, recoil = true, convention = :leading)
    @test got ≈ want

    # Omitting it on an unequal-mass class is a construction error, never a
    # silent fallback to the light (r = 1/2) form factor.
    @test_throws ArgumentError DecayChannel("Kstar_c", "D", "pi", 1.0, :A_c, 2)
    @test_throws ArgumentError DecayChannel("x", "y", "z", 1.0, :A, 1; heavy_fraction = 1.5)
    @test DecayChannel("rho", "pi", "pi", 1.0, :A, 1).heavy_fraction == 0.5
end

@testset "Quark types (mass is dynamics, charge is the only discrete datum)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    c = HeavyQuark{:up}(mq["c"], :c)
    b = HeavyQuark{:down}(mq["b"], :b)
    s = StrangeQuark(mq["s"])
    l = LightQuark(mq["q"])

    # Charge is the discrete flavor datum, carried by the weak-isospin tag.
    @test charge(c) == 2 // 3
    @test charge(b) == -1 // 3
    @test charge(s) == -1 // 3

    # LightQuark has NO charge: the model sets m_u = m_d, so it cannot resolve
    # up from down. Charge lives on the meson's flavor wavefunction instead
    # (pi+ = u dbar is charged, pi0 = (u ubar - d dbar)/sqrt(2) is not, at the
    # same constituent masses). This must raise, never guess.
    @test_throws MethodError charge(l)

    # The tag is weak-isospin class, validated at construction.
    @test_throws ArgumentError HeavyQuark{:sideways}(1.5, :x)
    @test_throws ArgumentError HeavyQuark{:up}(-1.0, :x)
    @test_throws ArgumentError LightQuark(0.0)
    @test_throws ArgumentError StrangeQuark(-0.4)

    @test flavor_symbol(l) === :q
    @test flavor_symbol(s) === :s
    @test flavor_symbol(c) === :c
    @test mass_GeV(c) == mq["c"]

    # Quark-built mesons agree with the flavor-table form.
    @test Meson(c, c) == Meson(mq, :c, :c)
    @test Meson(l, s) == Meson(mq, :q, :s)

    # Mass is the only dynamical input: same mass + different charge class ->
    # identical constituent masses, different charge.
    up_at_mb = HeavyQuark{:up}(mq["b"], :hypothetical)
    @test Meson(up_at_mb, up_at_mb).constituent_masses ==
          Meson(b, b).constituent_masses
    @test charge(up_at_mb) != charge(b)
end

@testset "The oscillator basis has a fixed phase" begin
    # LAPACK's QR picks its own column signs, and they varied with nbasis, beta
    # and the mesh (flips at n = 10 and n = 18 for beta = 0.65, nbasis = 24).
    # Invisible numerically -- a sign flip is unitary -- but fatal once a closed
    # form enters, since analytic elements are written in ho_reduced_radial's
    # convention. Wiring exact p^2 into an unfixed basis put charmonium 1S
    # 9.4 MeV BELOW the finite-difference answer, which a variational
    # calculation in a finite basis cannot do.
    r, h = GIModel.radial_grid(450, 24.0)
    for β in (0.25, 0.65, 1.35, 2.35), nb in (6, 12, 24)
        B = GIModel.ho_basis_matrix(0, β, r, nb)
        U = GIModel.orthonormalize_physical_basis(B, h)
        # every column keeps the sign of the basis function it came from
        @test all(sum(view(U, :, j) .* view(B, :, j)) > 0 for j in 1:nb)
        # and it is still orthonormal under the physical inner product
        @test maximum(abs, h * (transpose(U) * U) - I) < 1e-10
    end

    # With the phase fixed, the projected p^2 agrees with the closed form up to
    # mesh error, and that error is a property of the mesh vs beta: fine at the
    # variational optimum, poor where the basis is too compact or too diffuse
    # for the grid (the resolution wall and the beta railing, as numbers).
    ana(β, nb) = Matrix(ho_p2_matrix(0, β, nb))
    proj(β, nb) = Matrix(GIModel.projected_matrix(
        GIModel.orthonormalize_physical_basis(GIModel.ho_basis_matrix(0, β, r, nb), h),
        h, GIModel.p2_operator(1.0, 0, r, h)))
    rel(β, nb) = maximum(abs, proj(β, nb) - ana(β, nb)) / maximum(abs, ana(β, nb))
    @test rel(0.65, 12) < 0.01          # near the variational optimum
    @test rel(0.65, 24) < 0.01
    @test rel(2.35, 24) > 0.05          # too compact for h: the resolution wall
    @test rel(0.25, 24) > 1.0           # too diffuse for rmax: the beta railing
end

@testset "Oscillator matrix elements validate themselves" begin
    # These are a self-contained mathematics problem: matrix elements of
    # operators in the 3D oscillator basis. They can be — and here are —
    # validated with no mesh, no reference data and no reference to the GI
    # model at all. Doing that BEFORE any table comparison means a later
    # disagreement with the paper is about the physics, not about whether the
    # algebra is right.

    # 1. The strongest check available: p^2 and r^2 must RECONSTRUCT the
    #    oscillator Hamiltonian. H = p^2/2mu + (1/2) mu w^2 r^2 is diagonal with
    #    eigenvalue (2n+L+3/2)w, and beta^2 = mu*w, so
    #        p^2/2mu + (beta^4/2mu) r^2
    #    must be EXACTLY diagonal with exactly those eigenvalues. An error in
    #    either operator's magnitude, sign, or power of beta destroys the
    #    cancellation, so this tests both operators and their relative
    #    normalization simultaneously.
    for L in (0, 1, 2), β in (0.35, 0.65, 1.35, 2.35)
        nb, μ = 24, 0.7                      # any mu; beta^2 = mu*w fixes w
        ω = β^2 / μ
        H = Matrix(ho_p2_matrix(L, β, nb)) ./ (2μ) .+
            (β^4 / (2μ)) .* Matrix(ho_r2_matrix(L, β, nb))
        @test maximum(abs, H - Diagonal(diag(H))) < 1e-12          # exactly diagonal
        @test maximum(abs, diag(H) .- [(2n + L + 1.5) * ω for n = 0:(nb-1)]) < 1e-12
    end

    # 2. Virial theorem for the oscillator: <T> = <V> in every eigenstate.
    for L in (0, 2), β in (0.45, 1.85)
        nb, μ = 16, 1.3
        T = Matrix(ho_p2_matrix(L, β, nb)) ./ (2μ)
        V = (β^4 / (2μ)) .* Matrix(ho_r2_matrix(L, β, nb))
        for n in 1:nb
            @test isapprox(T[n, n], V[n, n]; rtol = 1e-12)
        end
    end

    # 3. Dimensional scaling: p^2 ~ beta^2, r^2 ~ 1/beta^2, so their product is
    #    beta-independent — a check no single operator can provide alone.
    for L in (0, 1)
        A = Matrix(ho_p2_matrix(L, 0.5, 10)) * Matrix(ho_r2_matrix(L, 0.5, 10))
        B = Matrix(ho_p2_matrix(L, 1.9, 10)) * Matrix(ho_r2_matrix(L, 1.9, 10))
        @test isapprox(A, B; rtol = 1e-10)
    end

    # 4. Structure: symmetric, tridiagonal, and r^2 positive definite.
    for op in (ho_p2_matrix(1, 0.8, 9), ho_r2_matrix(1, 0.8, 9))
        M = Matrix(op)
        @test M ≈ transpose(M)
        @test all(iszero, [M[i, j] for i in 1:9, j in 1:9 if abs(i - j) > 1])
        @test all(>(0), eigvals(Symmetric(M)))
    end
    @test_throws ArgumentError ho_r2_matrix(0, -1.0, 4)
    @test_throws ArgumentError ho_r2_matrix(0, 0.5, 0)
end

@testset "A17 position side: Gauss-Laguerre in Golub-Welsch form" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    # The Jacobi matrix for weight x^(L+1/2) e^{-x}, x = (beta r)^2, IS
    # beta^2 * ho_r2_matrix. Diagonalizing it gives nodes as eigenvalues and
    # sqrt(w_i) p_n(x_i) as eigenvector entries -- so the 1e19 polynomial value
    # and the 1e-95 weight never exist separately. That is the whole trick:
    # forming them apart loses the digits before any summation happens, which no
    # compensated summation can undo.
    for L in (0, 1, 2), β in (0.25, 0.65, 2.35)
        # g = 1 must be the identity, because Z is orthogonal
        @test maximum(abs, Matrix(ho_operator_matrix(L, β, 24, r -> 1.0)) - I) < 1e-13
        # g = r^2 must reconstruct the very Jacobi matrix that generated the rule
        @test maximum(abs, Matrix(ho_operator_matrix(L, β, 24, r -> r^2)) -
                           Matrix(ho_r2_matrix(L, β, 24))) < 1e-10
    end

    # Higher pure moments too, against the Gaussian-moment closed form
    # <n|r^0|n> = 1 and the r^2 case above; r^4 is checked for consistency
    # between two independent quadrature sizes rather than a closed form.
    for L in (0, 1), β in (0.45, 1.35)
        a = Matrix(ho_operator_matrix(L, β, 12, r -> r^4; nq = 64))
        b = Matrix(ho_operator_matrix(L, β, 12, r -> r^4; nq = 256))
        @test maximum(abs, a - b) < 1e-9
    end

    # The real operator: the Appendix-A smeared Coulomb + confinement, checked
    # against independent adaptive quadrature (QuadGK), not against itself.
    mc = mq["c"]
    g(r) = GIModel.smeared_coulomb_G_closed(params, mc, mc, r) +
           GIModel.smeared_confinement_S_closed(params, mc, mc, r)
    L, β, nb = 0, 0.65, 24
    M = Matrix(ho_operator_matrix(L, β, nb, g))
    for (a, b) in ((1, 1), (1, 2), (3, 3), (5, 8), (12, 12))
        ref, _ = quadgk(r -> GIModel.ho_reduced_radial(a - 1, L, β, r) * g(r) *
                             GIModel.ho_reduced_radial(b - 1, L, β, r),
                        0, Inf; rtol = 1e-13, order = 21)
        @test isapprox(M[a, b], ref; atol = 1e-11)
    end

    # Symmetry, and independence of the quadrature size once converged.
    @test M ≈ transpose(M)
    @test maximum(abs, M - Matrix(ho_operator_matrix(L, β, nb, g; nq = 512))) < 1e-10
    @test_throws ArgumentError ho_operator_matrix(0, -1.0, 4, r -> 1.0)
    @test_throws ArgumentError ho_operator_matrix(0, 0.5, 0, r -> 1.0)

    # The rule itself carries no beta. In x = (beta r)^2 the Jacobi matrix is
    # beta^2 * ho_r2_matrix(L, beta, n) = ho_r2_matrix(L, 1, n), every beta
    # cancelling, so nodes and DVR weights are functions of (L, nq) alone and
    # beta enters only as r_i = sqrt(x_i)/beta. `gauss_laguerre_dvr` memoizes on
    # that fact; if it ever stopped holding, the cache would silently hand one
    # beta's rule to another.
    for L in (0, 1, 2), β in (0.25, 0.65, 2.35)
        @test Matrix(β^2 * ho_r2_matrix(L, β, 40)) ≈ Matrix(ho_r2_matrix(L, 1, 40))
    end
    # The memo is capped by BYTES, not entry count: entries span 24x64 to
    # 24x8192, so 256 small ones cost less than one large one. Flushing can only
    # cost time -- a rebuilt entry is bit-identical, which is what makes the
    # whole memo safe to discard at any moment.
    let c = GIModel._GAUSS_LAGUERRE_DVR
        a = GIModel.gauss_laguerre_dvr(0, 24, 1024)
        empty!(c)
        b = GIModel.gauss_laguerre_dvr(0, 24, 1024)
        @test a[1] == b[1] && a[2] == b[2]
        @test GIModel._dvr_cache_bytes() == sizeof(b[1]) + sizeof(b[2])
    end

    sqrt_x, Z = GIModel.gauss_laguerre_dvr(1, 12, 64)
    @test Z * transpose(Z) ≈ I            # orthogonality of the leading rows
    @test issorted(sqrt_x)                # nodes come out ordered
    @test length(sqrt_x) == 64 && size(Z) == (12, 64)
end

@testset "Exact oscillator p^2 (A17 momentum side)" begin
    # p^2 = 2*mu*H_osc - beta^4 r^2 is tridiagonal in the oscillator basis:
    #   <n|p^2|n>   = beta^2 (2n + L + 3/2)
    #   <n|p^2|n+1> = beta^2 sqrt((n+1)(n + L + 3/2))
    # Checked against the mesh projection it is meant to replace -- and the
    # check is CONVERGENCE, not a tolerance: the difference must shrink as the
    # mesh refines, which is what shows the formula is exact and the mesh is
    # the approximation.
    for (L, β, nb) in ((0, 0.55, 6), (1, 0.75, 6), (2, 0.45, 5))
        ana = Matrix(ho_p2_matrix(L, β, nb))
        errs = Float64[]
        for (ng, rm) in ((450, 24.0), (2000, 40.0), (8000, 56.0))
            r, h = GIModel.radial_grid(ng, rm)
            U = GIModel.orthonormalize_physical_basis(GIModel.ho_basis_matrix(L, β, r, nb), h)
            num = Matrix(GIModel.projected_matrix(U, h, GIModel.p2_operator(1.0, L, r, h)))
            push!(errs, maximum(abs, num .- ana))
        end
        @test errs[1] > errs[2] > errs[3]        # monotone convergence to the formula
        @test errs[3] < 1e-3
    end

    # Structure: symmetric, tridiagonal, positive definite (p^2 is).
    M = Matrix(ho_p2_matrix(0, 0.6, 8))
    @test M ≈ transpose(M)
    @test all(iszero, [M[i, j] for i in 1:8, j in 1:8 if abs(i - j) > 1])
    @test all(>(0), eigvals(Symmetric(M)))
    # beta scaling: p^2 has dimensions of beta^2.
    @test Matrix(ho_p2_matrix(0, 1.2, 5)) ≈ 4 .* Matrix(ho_p2_matrix(0, 0.6, 5))
    @test_throws ArgumentError ho_p2_matrix(0, -1.0, 4)
    @test_throws ArgumentError ho_p2_matrix(0, 0.5, 0)
end

@testset "RadialWave interface: invariants, not fixed numbers" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r, h = GIModel.radial_grid(450, 24.0)

    @test MeshWave <: RadialWave

    # --- closed-form cross-checks: the operations must reproduce integrals we
    # --- can do by hand, independently of any solver.
    # u = r e^{-r}  =>  <r^n> = Gamma(3+n)/(2^n Gamma(3)) exactly.
    u = r .* exp.(-r)
    w = MeshWave(u ./ sqrt(sum(abs2, u) * h), r)
    @test isapprox(wave_norm(w), 1.0; rtol = 1e-10)
    @test isapprox(radial_expect(w, x -> 1.0), 1.0; rtol = 1e-10)
    @test isapprox(radial_expect(w, x -> x), 1.5; rtol = 1e-5)          # 3/2
    @test isapprox(radial_expect(w, x -> x^2), 3.0; rtol = 1e-5)        # 3
    # A normalized wave overlapped with itself is 1, and the overlap is
    # invariant under rescaling either wave (it normalizes internally).
    @test isapprox(radial_overlap(w, w, x -> 1.0), 1.0; rtol = 1e-10)
    w_scaled = MeshWave(7.3 .* w.u, r)
    @test isapprox(radial_overlap(w, w_scaled, x -> 1.0), 1.0; rtol = 1e-10)
    @test isapprox(radial_expect(w_scaled, x -> x), radial_expect(w, x -> x); rtol = 1e-12)

    # Two different waves: overlap is symmetric and bounded by 1 (Cauchy-Schwarz).
    v = r .* exp.(-0.7 .* r)
    wv = MeshWave(v ./ sqrt(sum(abs2, v) * h), r)
    ov = radial_overlap(w, wv, x -> 1.0)
    @test isapprox(ov, radial_overlap(wv, w, x -> 1.0); rtol = 1e-12)
    @test abs(ov) <= 1.0 + 1e-10

    # --- Parseval: the transform must conserve probability. This ties the
    # --- position-space and momentum-space halves of the interface together.
    for key in ("q", "c", "b")
        masses = ConstituentMasses(mq[key], mq[key])
        sol = channel_solution(params, masses, 0; nlevels = 2)
        for n in 1:2
            wn = radial_wave(sol, n)
            mw = momentum_wave(wn, 0)
            @test isapprox(wave_norm(wn), 1.0; rtol = 1e-10)
            @test isapprox(momentum_expect(mw, p -> 1.0), 1.0; rtol = 1e-6)
            # <E> >= m for a relativistic quark energy, with equality only at p=0.
            @test momentum_expect(mw, p -> sqrt(mq[key]^2 + p^2)) > mq[key]
        end
    end

    # --- orthogonality: distinct radial levels of one channel are orthogonal.
    # --- Nothing in the solve enforces this; it is a property of the operator,
    # --- so it is a real check on the solve rather than on the interface.
    masses = ConstituentMasses(mq["c"], mq["c"])
    sol = channel_solution(params, masses, 0; nlevels = 3)
    for i in 1:3, j in 1:3
        w_i = radial_wave(sol, i)
        w_j = radial_wave(sol, j)
        target = i == j ? 1.0 : 0.0
        @test isapprox(abs(radial_overlap(w_i, w_j, x -> 1.0)), target; atol = 1e-8)
    end

    # --- node counting: the n-th radial level has n-1 interior nodes. Count only
    # --- where the wave is physically present: the exponential tail flips sign
    # --- in floating point (measured: crossings at r = 17.6 and 23.8 with |u|
    # --- at 2e-12 and 2e-16 of the peak), so a naive sign count reports phantom
    # --- nodes in the ground state.
    for n in 1:3
        u_n = radial_wave(sol, n).u
        cut = 1e-8 * maximum(abs, u_n)
        big = [i for i in eachindex(u_n) if abs(u_n[i]) > cut]
        signs = sign.(u_n[big])
        nodes = count(i -> signs[i] != signs[i+1], 1:(length(signs)-1))
        @test nodes == n - 1
    end

    # --- momentum space is its own noun, with a linear and a quadratic operation.
    # --- momentum_expect is quadratic in Phi, momentum_functional is linear;
    # --- conflating them is a whole class of factor-of-Phi bug.
    @test MeshMomentumWave <: MomentumWave
    mc = mq["c"]
    solq = channel_solution(params, ConstituentMasses(mc, mc), 0; nlevels = 1)
    wq = radial_wave(solq, 1)
    mwq = momentum_wave(wq, 0)
    # quadratic with g = 1 is the norm; linear with K = 1 is NOT (different object)
    @test isapprox(momentum_expect(mwq, p -> 1.0), 1.0; rtol = 1e-6)
    @test !isapprox(momentum_functional(mwq, p -> 1.0), 1.0; rtol = 1e-3)
    # both are linear in their kernel argument
    @test isapprox(momentum_expect(mwq, p -> 3.0), 3 * momentum_expect(mwq, p -> 1.0); rtol = 1e-12)
    @test isapprox(momentum_functional(mwq, p -> 3.0), 3 * momentum_functional(mwq, p -> 1.0); rtol = 1e-12)
    # the functional scales linearly in Phi, the expectation quadratically
    scaled = MeshMomentumWave(mwq.p, 2 .* mwq.phi)
    @test isapprox(momentum_functional(scaled, p -> 1.0), 2 * momentum_functional(mwq, p -> 1.0); rtol = 1e-12)
    @test isapprox(momentum_expect(scaled, p -> 1.0), 4 * momentum_expect(mwq, p -> 1.0); rtol = 1e-12)

    # --- scale invariance: every consumer of a wave must give the same answer
    # --- for a rescaled wave, because the amplitude of u carries no physics.
    # --- This is what catches a deleted normalization: migrating
    # --- charge_radius_squared onto `radial_expect` removed the normalization
    # --- that `_rel_momentum_average` was silently relying on, and only a
    # --- rescaled input revealed it (the value moved by a factor of 13.7^2).
    sol3 = channel_solution(params, ConstituentMasses(mq["q"], mq["s"]), 0; nlevels = 1)
    w_one = radial_wave(sol3, 1)
    w_big = MeshWave(13.7 .* w_one.u, w_one.r)
    for probe in (
        w -> charge_radius_squared(w, mq["q"], 2 // 3, mq["s"], 1 // 3),
        w -> radial_expect(w, x -> x^2),
        w -> wave_norm(w) / wave_norm(w),
        w -> momentum_expect(momentum_wave(w, 0), p -> 1.0),
    )
        @test isapprox(probe(w_one), probe(w_big); rtol = 1e-10)
    end

    # Mismatched meshes are an error, not a silently wrong overlap.
    other, _ = GIModel.radial_grid(200, 24.0)
    @test_throws ArgumentError radial_overlap(w, MeshWave(other .* 0 .+ 1.0, other), x -> 1.0)
end

@testset "Every solve returns u with the same normalization" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    # The central solve, with both solvers and every orbital: integral u^2 dr = 1.
    # FD used to return Euclidean eigenvectors here (sum u^2 = 1) while HO
    # returned physical ones, differing by sqrt(h).
    for key in ("q", "s", "c", "b"), L in 0:2
        masses = ConstituentMasses(mq[key], mq[key])
        for slv in (FiniteDifferenceSolver(), OscillatorSolver())
            sol = channel_solution(params, masses, L; solver = slv, nlevels = 2)
            for wave in sol.waves
                @test isapprox(wave_norm(wave), 1.0; rtol = 1e-10)
            end
            if slv isa OscillatorSolver
                @test sol.convergence.status == :converged
                @test sol.convergence.energy_delta_GeV <= slv.energy_tolerance_GeV
                @test sol.convergence.refinements >=
                      GIModel.HO_REQUIRED_CONVERGED_REFINEMENTS
                @test sol.convergence.max_wave_overlap_defect >= 0
            else
                @test isnothing(sol.convergence)
            end
        end
    end

    # The convention has one definition, and it is idempotent.
    raw = reshape(collect(1.0:10.0), 10, 1)
    once = GIModel.physically_normalized_waves(raw, 0.25)
    @test isapprox(sum(abs2, once) * 0.25, 1.0; rtol = 1e-12)
    @test GIModel.physically_normalized_waves(once, 0.25) ≈ once
    # A zero column is left alone rather than producing NaN.
    @test all(iszero, GIModel.physically_normalized_waves(zeros(5, 1), 0.25))
end

@testset "Native spin matrices equal wave-interface expectations" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["c"], mq["c"])
    beta, nbasis = 0.75, 12

    cS = collect(range(1.0, 0.1; length = nbasis))
    wS = OscillatorWave(0, beta, cS)
    contact = GIModel.ho_contact_hyperfine_matrix(
        params, masses, 0, 1, beta, nbasis,
    )
    @test Matrix(contact) ≈ Matrix(contact)' atol = 1e-13
    @test dot(wS.coefficients, contact * wS.coefficients) ≈
          GIModel.contact_hyperfine_shift_momentum_sandwich(
              params, masses, FineStructureMultiplet("S", 1, 0), wS,
          ) rtol = 1e-5 atol = 1e-9

    cP = [sin(i) + 0.2cos(2i) for i in 1:nbasis]
    wP = OscillatorWave(1, beta, cP)
    matrices = GIModel.ho_fine_structure_matrices(
        params, masses, 1, 3, 2, beta, nbasis,
    )
    components = fine_structure_components(
        params, masses, FineStructureMultiplet("P", 3, 2), wP,
    )
    @test all(
        M -> isapprox(Matrix(M), Matrix(M)'; atol = 1e-13),
        (matrices.spin_orbit_vector, matrices.spin_orbit_thomas,
         matrices.spin_orbit, matrices.tensor, matrices.total),
    )
    @test dot(wP.coefficients, matrices.spin_orbit_vector * wP.coefficients) ≈
          components.spin_orbit_vector rtol = 1e-5 atol = 1e-9
    @test dot(wP.coefficients, matrices.spin_orbit_thomas * wP.coefficients) ≈
          components.spin_orbit_thomas rtol = 1e-5 atol = 1e-9
    @test dot(wP.coefficients, matrices.tensor * wP.coefficients) ≈
          components.tensor rtol = 1e-5 atol = 1e-9
    @test Matrix(matrices.total) ≈
          Matrix(matrices.spin_orbit) + Matrix(matrices.tensor) atol = 1e-13

    r, h = GIModel.radial_grid(120, 12.0)
    wfd = MeshWave(r .^ 2 .* exp.(-r), r)
    fdmat = GIModel.fine_structure_grid_matrices(params, masses, 2, r, h; L = 1)
    fdcomp = fine_structure_components(
        params, masses, FineStructureMultiplet("P", 3, 2), wfd,
    )
    expect(M) = h * dot(wfd.u, M * wfd.u) / wave_norm(wfd)
    @test expect(fdmat.spin_orbit_vector) ≈ fdcomp.spin_orbit_vector atol = 1e-10
    @test expect(fdmat.spin_orbit_thomas) ≈ fdcomp.spin_orbit_thomas atol = 1e-10
    @test expect(fdmat.tensor) ≈ fdcomp.tensor atol = 1e-10
    @test Matrix(fdmat.total) ≈
          Matrix(fdmat.spin_orbit) + Matrix(fdmat.tensor) atol = 1e-13
end

@testset "Both solvers solve the same native fixed-(L,S,J) Hamiltonian" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r, h = GIModel.radial_grid(450, 24.0)

    for key in ("q", "c", "b")
        masses = ConstituentMasses(mq[key], mq[key])
        multiplet = FineStructureMultiplet("S", 1, 0)
        terms = SpinTerms(fine_structure = false)
        sol_fd = fixed_channel_solution(
            params, masses, multiplet;
            solver = FiniteDifferenceSolver(), terms = terms, nlevels = 2)
        sol_ho = fixed_channel_solution(
            params, masses, multiplet;
            solver = OscillatorSolver(), terms = terms, nlevels = 2)

        # Same normalization, independent of representation.
        # FD used to return Euclidean eigenvectors (sum u^2 = 1), differing from
        # HO by exactly sqrt(h); anything quadratic in u was then off by h.
        for wave in vcat(sol_fd.waves, sol_ho.waves)
            @test isapprox(wave_norm(wave), 1.0; rtol = 1e-10)
        end
        # Same quantity, to the combined accuracy of two now-INDEPENDENT methods.
        # This was 1e-3 while the oscillator path projected the finite-difference
        # p^2 through the mesh: it was then a Galerkin restriction of the FD
        # problem, inherited FD's discretization error, and so agreed artificially
        # well. With Eq. (A17) assembled from exact matrix elements the two share
        # no error, and the gap measures both: +0.18 MeV (charm) and +0.55 MeV
        # (bottom) with the oscillator answer correctly ABOVE — variational in a
        # finite basis — and ~1.6 MeV for light quarks, where the finite-difference
        # mesh is itself least converged and sits above the true value.
        @test abs(sol_ho.eigenvalues_GeV[1] - sol_fd.eigenvalues_GeV[1]) < 3e-3
    end

    # A mesh operator is not an HO input. Unsupported routes fail explicitly.
    masses = ConstituentMasses(mq["c"], mq["c"])
    V = Matrix(GIModel.contact_hyperfine_operator(params, masses, "S", 1, r))
    ho = OscillatorSolver()
    @test_throws ArgumentError resummed_channel_solution(
        params, masses, 0, V; solver = ho, nlevels = 2)

    # V must live on the solver's mesh, whichever solver that is.
    bad = zeros(10, 10)
    @test_throws Exception resummed_channel_solution(
        params, masses, 0, bad; solver = FiniteDifferenceSolver(), nlevels = 1)
    @test_throws Exception resummed_channel_solution(
        params, masses, 0, bad; solver = ho, nlevels = 1)
end

@testset "An unimplemented solve fails instead of degrading" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["q"], mq["q"])
    r, _ = GIModel.radial_grid(450, 24.0)

    # FD has the non-perturbative contact solve, and returns it.
    solution = contact_hyperfine_nonperturbative_states(params, masses, "S", 1, r, 2)
    @test length(solution.eigenvalues_GeV) == 2 && length(solution.waves) == 2

    # Empty is still the right answer where the path is genuinely inactive:
    # not an S wave, or a multiplicity the contact term does not touch.
    @test contact_hyperfine_nonperturbative_states(params, masses, "P", 1, r, 2) === nothing
    @test GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 2, r, 2) == Float64[]

    # The oscillator basis used to answer "empty" here, which
    # the old spectrum correction stage read as "fall back to first-order PT": the light
    # 1S0 came out 0.2842 GeV against the correctly resummed 0.149 GeV. U1 made that a loud
    # failure, U2 unified the solve, and B1 made this wrapper basis-generic — so
    # the oscillator path now resums in its own space and lands on the same
    # answer. That agreement is what the throw was standing in for.
    ho = OscillatorSolver()
    solution_ho = contact_hyperfine_nonperturbative_states(
        params, masses, "S", 1, 2; solver = ho)
    @test length(solution_ho.eigenvalues_GeV) == 2 && length(solution_ho.waves) == 2
    @test abs(solution_ho.eigenvalues_GeV[1] - solution.eigenvalues_GeV[1]) < 1e-3
    @test 0.12 < solution_ho.eigenvalues_GeV[1] < 0.18
    @test GIModel.contact_hyperfine_nonperturbative_levels(
        params, masses, "S", 1, 2; solver = ho) ≈ solution_ho.eigenvalues_GeV
    @test_throws ArgumentError contact_hyperfine_nonperturbative_states(
        params, masses, "S", 1, r, 2; solver = ho)
    ho_pi = spectrum_state(
        compute_spectrum(params, Meson(mq, :q, :q); solver = ho, levels = spectrum_levels(1)),
        "1^1S_0")
    @test 0.12 < ho_pi.mass_GeV < 0.18

    # The FD front door is untouched and still resums.
    fd_pi = spectrum_state(
        compute_spectrum(params, Meson(mq, :q, :q); levels = spectrum_levels(1)), "1^1S_0")
    @test 0.12 < fd_pi.mass_GeV < 0.18 # resummed; first-order PT lands near 0.28
end

@testset "RadialSolver is numerics, SpinTerms is physics" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    meson = Meson(mq, :q, :s)          # unequal mass: every mixing mechanism is live
    levels = spectrum_levels(2)
    base = compute_spectrum(params, meson; levels = levels)
    masses(spec) = [s.mass_GeV for s in spec.states]

    # The struct defaults ARE the old inline defaults: passing them explicitly
    # must reproduce the default call bit for bit.
    @test masses(compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 450, rmax = 24.0, kinetic = :relativistic,
            eigensolver = :full, nlevels_per_channel = 6),
        terms = SpinTerms(contact_hyperfine = true, fine_structure = true,
            same_j_spin_orbit = true, tensor = true))) == masses(base)
    @test masses(compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(), terms = SpinTerms())) == masses(base)

    # RadialSolver is "how well": refining the mesh must not move a mass more
    # than the discretization error it removes (sub-MeV here).
    fine = compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 900, rmax = 32.0))
    @test maximum(abs.(masses(fine) .- masses(base))) < 0.003   # < 3 MeV

    # SpinTerms is "what physics": every switch must move a mass, and turning
    # them all off must land exactly on the central eigenvalues.
    for sw in (:contact_hyperfine, :fine_structure, :same_j_spin_orbit, :tensor)
        off = compute_spectrum(params, meson; levels = levels,
            terms = SpinTerms(; sw => false))
        @test maximum(abs.(masses(off) .- masses(base))) > 1e-4   # > 0.1 MeV
    end
    none = compute_spectrum(params, meson; levels = levels,
        terms = SpinTerms(contact_hyperfine = false, fine_structure = false,
            same_j_spin_orbit = false, tensor = false))
    @test all(s.mass_GeV === s.central_GeV for s in none.states)

    # The solver threads down to the low tier, and its copy constructor keeps
    # the untouched fields.
    sol = channel_solution(params, meson.constituent_masses, 0;
        solver = FiniteDifferenceSolver(nlevels_per_channel = 3))
    @test length(sol.eigenvalues_GeV) == 3
    tuned = FiniteDifferenceSolver(FiniteDifferenceSolver(ngrid = 900); rmax = 32.0)
    @test tuned.ngrid == 900 && tuned.rmax == 32.0 && tuned.kinetic === :relativistic

    # Every path is silent. There is nothing left to deprecate: the loose
    # keywords are gone, so the shim that used to fold them in -- and warn that
    # they silently overrode the objects -- is gone with them. A settings name
    # that no longer exists is now a MethodError at the call, not a mass that
    # quietly came from the wrong grid.
    @test_logs compute_spectrum(params, meson; levels = levels)
    @test_logs compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(), terms = SpinTerms())
    @test_logs central_spectrum(params, meson; levels = levels)
    @test_logs channel_solution(params, meson.constituent_masses, 0;
        solver = FiniteDifferenceSolver())
    # Every retired keyword, by its old name, on the entry point that used to
    # accept it. NOTE: these calls are deliberately written the old way -- that
    # they no longer parse into a method IS the assertion -- so do not rewrite
    # them into the solver/terms form.
    corrected = fixed_spectrum(params, meson; levels = levels)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels, ngrid = 450)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels, rmax = 24.0)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        eigensolver = :krylov)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        nlevels_per_channel = 3)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        tensor_mixing = true)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        contact_hyperfine = true)
    @test_throws MethodError compute_spectrum(params, meson; levels = levels,
        same_j_spin_orbit_mixing = true)
    @test_throws MethodError central_spectrum(params, meson; levels = levels, rmax = 24.0)
    @test_throws MethodError channel_solution(params, meson.constituent_masses, 0;
        kinetic = :relativistic)
    @test_throws MethodError solve_sector(params, 1.628; ngrid = 100)
    @test_throws MethodError fixed_spectrum(
        params, meson; levels = levels, use_fine_structure = false,
    )
    @test_throws MethodError add_intra_meson_mixing(corrected; tensor_mixing = false)

    # Invalid settings are construction errors, not silent fallbacks.
    @test_throws ArgumentError RadialSolver(kinetic = :newtonian)
    @test_throws ArgumentError RadialSolver(eigensolver = :lanczos)
    @test_throws ArgumentError RadialSolver(ngrid = 1)
    @test_throws ArgumentError RadialSolver(rmax = 0)
    @test_throws ArgumentError RadialSolver(nlevels_per_channel = 0)
end

@testset "The solver, not the parameters, chooses the radial method" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    meson = Meson(mq, :c, :c)
    levels = spectrum_levels(2)
    masses(spec) = [s.mass_GeV for s in spec.states]

    # `GIParameters` is the model and nothing else: no basis marker, and one
    # object serves both methods. It used to carry a phantom type parameter that
    # no field used, so switching numerical method meant rebuilding the physics.
    @test isconcretetype(typeof(params))
    @test !any(f -> occursin("asis", String(f)), fieldnames(GIParameters))

    # `RadialSolver` is the interface; calling it builds the default
    # implementation, which is what every pre-split call site meant.
    @test RadialSolver isa Type && !isconcretetype(RadialSolver)
    @test RadialSolver(ngrid = 900) === FiniteDifferenceSolver(ngrid = 900)
    @test FiniteDifferenceSolver() isa RadialSolver
    @test OscillatorSolver() isa RadialSolver

    # Same parameters, two methods, through the front door. Both must answer,
    # and agree to a few MeV -- the residual is each one's own truncation error.
    fd = compute_spectrum(params, meson; levels = levels, solver = FiniteDifferenceSolver())
    ho = compute_spectrum(params, meson; levels = levels, solver = OscillatorSolver())
    @test isconcretetype(typeof(meson))
    @test isconcretetype(typeof(fd.computation))
    @test isconcretetype(typeof(fd))
    @test maximum(abs.(masses(fd) .- masses(ho))) < 0.005

    # The oscillator path is variational in a finite basis, so it can only sit
    # ABOVE the converged answer. Against a fine finite-difference mesh, every
    # charmonium level must therefore come out no lower.
    fine = compute_spectrum(params, meson; levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 1800, rmax = 24.0))
    @test all(masses(ho) .>= masses(fine) .- 1e-9)

    # Raising nbasis can only lower an oscillator eigenvalue, for the same reason.
    e24 = channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(nbasis = 24, converge = false),
        nlevels = 3).eigenvalues_GeV
    e40 = channel_solution(params, meson.constituent_masses, 0;
        solver = OscillatorSolver(nbasis = 40, converge = false),
        nlevels = 3).eigenvalues_GeV
    @test all(e40 .<= e24 .+ 1e-9)

    # The stage-1 solver is recorded and stages 2-3 reuse it, so one spectrum is
    # one calculation. The non-perturbative contact solve in stage 2 resolves its
    # own eigenproblem and would otherwise silently revert to finite differences.
    @test compute_spectrum(params, meson; levels = levels,
        solver = OscillatorSolver()).computation.solver isa OscillatorSolver

    # Settings the oscillator path does not have are absent, not ignored. The
    # `:nonrelativistic` comparator is a finite-difference thing, and asking an
    # OscillatorSolver for it fails at construction rather than silently doing
    # the relativistic calculation you did not ask for.
    @test !hasfield(OscillatorSolver, :kinetic)
    @test !hasfield(OscillatorSolver, :eigensolver)
    @test_throws MethodError OscillatorSolver(kinetic = :nonrelativistic)
    @test_throws MethodError OscillatorSolver(eigensolver = :krylov)
    @test_throws ArgumentError OscillatorSolver(nbasis = 0)
    @test_throws ArgumentError OscillatorSolver(max_nbasis = 1)
    @test_throws ArgumentError OscillatorSolver(basis_step = 0)
    @test_throws ArgumentError OscillatorSolver(energy_tolerance_GeV = 0)
    @test_throws ArgumentError OscillatorSolver(beta_tolerance_GeV = 0)
    @test_throws ArgumentError OscillatorSolver(beta_grid = Float64[])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [0.5, 0.6])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [0.5, -0.5])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [1.5, 0.5])
    @test_throws ArgumentError OscillatorSolver(beta_grid = [0.5, 0.5, 0.6])

    # The paper's complete fixed-sector Hamiltonian is relativistic. The
    # nonrelativistic FD option remains a central-only comparator and must not
    # be silently upgraded to a relativistic production calculation.
    @test_throws ArgumentError compute_spectrum(
        params, meson; levels = levels,
        solver = FiniteDifferenceSolver(kinetic = :nonrelativistic),
    )

    # beta_grid is a solver field, not a module constant to edit in source. A
    # bracket that excludes the minimum now fails rather than returning a
    # warned-but-usable under-resolved result.
    narrow = OscillatorSolver(beta_grid = 0.9:0.1:1.2)
    @test narrow.beta_grid == [0.9, 1.0, 1.1, 1.2]
    @test_throws ErrorException channel_solution(
        params, meson.constituent_masses, 0; solver = narrow, nlevels = 2)

    @testset "PA-12 adaptive HO convergence is certified or fails loudly" begin
        # Explicit fixed-size mode is available only as an unchecked convergence
        # study. Production mode records the continuously refined beta and final
        # basis certificate on the existing solution object.
        unchecked = channel_solution(
            params, meson.constituent_masses, 0;
            solver = OscillatorSolver(nbasis = 24, converge = false),
            nlevels = 2,
        )
        @test unchecked.convergence.status == :unchecked
        @test isnothing(unchecked.convergence.energy_delta_GeV)
        certified = channel_solution(
            params, meson.constituent_masses, 0;
            solver = OscillatorSolver(), nlevels = 2,
        )
        @test certified.convergence.status == :converged
        @test certified.convergence.beta_GeV == radial_wave(certified, 1).beta
        @test certified.convergence.beta_GeV ∉ OscillatorSolver().beta_grid
        @test certified.convergence.energy_delta_GeV <=
              certified.convergence.tolerance_GeV
        high_order_left = OscillatorWave(0, 0.55, [sin(i) for i in 1:64])
        high_order_right = OscillatorWave(0, 0.65, [cos(i / 2) for i in 1:72])
        @test isfinite(radial_overlap(high_order_left, high_order_right, _ -> 1.0))
        high_order_momentum = momentum_wave(high_order_left)
        @test isfinite(momentum_functional(
            high_order_momentum, p -> inv(sqrt(1 + p^2)),
        ))
        @test isfinite(momentum_overlap(
            high_order_momentum, momentum_wave(high_order_right), _ -> 1.0,
        ))
        too_small = OscillatorSolver(
            nbasis = 8,
            max_nbasis = 16,
            basis_step = 4,
            energy_tolerance_GeV = 1e-14,
            beta_grid = [0.8],
            nlevels_per_channel = 1,
        )
        @test_throws ErrorException channel_solution(
            params, meson.constituent_masses, 0; solver = too_small, nlevels = 1)
        @test_throws ArgumentError channel_solution(
            params,
            meson.constituent_masses,
            0;
            solver = OscillatorSolver(nbasis = 1, converge = false),
            nlevels = 2,
        )
    end
end

@testset "Resolution walls are detected, not silent" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    # Grid points spanned by the ground state's RMS radius. The FD grid is fixed,
    # so a compact enough state falls between points and comes back silently
    # under-resolved (at 30 GeV the hyperfine splitting collapses to 0.0025).
    function points_across(m; ngrid = 450, rmax = 24.0)
        sol = channel_solution(
            params, ConstituentMasses(m, m), 0;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax),
        )
        wave = radial_wave(sol, 1)
        rms = sqrt(radial_expect(wave, x -> x^2))
        return rms / wave.h
    end
    # Every sector the paper actually uses must sit clear of the threshold.
    for m in (mq["q"], mq["s"], mq["c"], mq["b"])
        @test points_across(m) > GIModel.MIN_POINTS_ACROSS_STATE
    end
    # ...including bottomonium on the deliberately coarser grid the Table III
    # mixing audit uses (ngrid = 220, rmax = 22): 10.3 points is ~2 MeV of
    # error, imprecise but not broken, so it must NOT warn.
    @test points_across(mq["b"]; ngrid = 220, rmax = 22.0) > GIModel.MIN_POINTS_ACROSS_STATE
    # The genuinely broken case must be caught: at 6.9 points the hyperfine
    # splitting collapses.
    @test points_across(30.0) < GIModel.MIN_POINTS_ACROSS_STATE

    # Same story on the oscillator path: the default beta_grid rails only well
    # above bottomonium. The adaptive production solve rejects that bracket.
    default_grid = OscillatorSolver().beta_grid
    function best_beta(m)
        mm = ConstituentMasses(m, m)
        r, h = GIModel.radial_grid(450, 24.0)
        best = nothing
        for b in default_grid
            H, _ = GIModel.oscillator_hamiltonian_for_beta(params, mm, 0, r, h, b; nbasis = 24)
            v, _ = GIModel.lowest_eigenpairs(Matrix(H), 2, Val(:full))
            (isnothing(best) || v[end] < best[2]) && (best = (b, v[end]))
        end
        return best[1]
    end
    @test best_beta(mq["b"]) < last(default_grid)     # bottomonium is fine
    @test best_beta(30.0) == last(default_grid)       # railed

    @test_throws ErrorException channel_solution(
        params, ConstituentMasses(30.0, 30.0), 0;
        solver = OscillatorSolver(), nlevels = 2)
end

@testset "Isoscalar coherence factor is stated, not string-matched" begin
    r = collect(range(0.05, 6.0; length = 64))
    wave = MeshWave(exp.(-r), r)
    mk(label, coh) = pseudoscalar_annihilation_basis_input(
        BasisState(1, "S", 1, 0; label = label), 0.22, 0.9, wave;
        isoscalar_coherent = coh)

    # The sqrt(2) follows the flag...
    @test GIModel._flavor_coherence_factor(mk("1 ns", true)) ≈ sqrt(2)
    @test GIModel._flavor_coherence_factor(mk("1 ss", false)) == 1.0

    # ...and NOT the label. Previously the factor was chosen by
    # occursin("ns", lowercase(label)), so a channel whose name merely contained
    # those letters silently gained a 41% amplitude factor. Both directions must
    # now be decided by the flag alone.
    @test GIModel._flavor_coherence_factor(mk("1 snsn", false)) == 1.0
    @test GIModel._flavor_coherence_factor(mk("1 cc", true)) ≈ sqrt(2)

    @test_throws UndefKeywordError pseudoscalar_annihilation_basis_input(
        BasisState(1, "S", 1, 0; label = "1 ns"), 0.22, 0.9, wave)
end

@testset "Meson construction and flavor resolution" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    cc = Meson(mq, :c, :c)
    @test cc.constituent_masses.m1_GeV ≈ mq["c"]
    @test is_equal_flavor(cc)
    @test reduced_mass(cc) ≈ mq["c"] / 2
    bu = Meson(mq, :b, :u)
    @test !is_equal_flavor(bu)
    @test bu.constituent_masses.m1_GeV ≈ mq["b"]
    @test bu.constituent_masses.m2_GeV ≈ mq["u"]
    # :n aliases the light average :q, matching the paper's n nbar notation
    nn = Meson(mq, :n, :n)
    qq = Meson(mq, :q, :q)
    @test nn == qq
    @test nn.constituent_masses.m1_GeV ≈ 0.5 * (mq["u"] + mq["d"])
    # unknown flavor fails loudly - no silent fallback masses
    @test_throws ArgumentError Meson(mq, :t, :t)
    @test_throws ArgumentError Meson(mq, :cbar, :c)
    # explicit-mass constructor for parameter scans
    scan = Meson(:c, :c, ConstituentMasses(1.5, 1.5))
    @test scan.constituent_masses.m1_GeV ≈ 1.5
end

@testset "spectrum_levels enumerates n^(2S+1)L_J multiplets" begin
    levels = spectrum_levels(2)
    # per n: S gives {1S0, 3S1}, P and D give {singlet + 3 triplets}
    @test length(levels) == 2 * (2 + 4 + 4)
    s_triplets = [l for l in levels if l.L_label == "S" && l.multiplicity == 3]
    @test all(l.J == 1 for l in s_triplets)
    p_triplets = [l for l in levels if l.L_label == "P" && l.multiplicity == 3 && l.n == 1]
    @test sort([l.J for l in p_triplets]) == [0, 1, 2]
    p_singlets = [l for l in levels if l.L_label == "P" && l.multiplicity == 1]
    @test all(l.J == 1 for l in p_singlets)
    @test_throws ArgumentError spectrum_levels(0)
    @test_throws ArgumentError spectrum_levels(2; L_labels = ("X",))
    only_s = spectrum_levels(3; L_labels = ("S",))
    @test length(only_s) == 6
end

@testset "compute_spectrum breakdown and tensor mixing (no reference data)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    meson = Meson(mq, :c, :c)
    levels = [
        BasisState(1, "S", 3, 1),
        BasisState(2, "S", 3, 1),
        BasisState(1, "D", 3, 1),
        BasisState(1, "P", 3, 2),
    ]
    spec = compute_spectrum(
        params, meson;
        levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
    )
    @test length(spec.states) == 4
    for s in spec.states
        if isempty(s.mixings)
            @test s.mass_GeV ≈ s.central_GeV + s.contact_shift_GeV + s.fine_structure_shift_GeV atol = 1e-12
        end
        @test s.fine_structure_shift_GeV ≈ s.spin_orbit_shift_GeV + s.tensor_shift_GeV atol = 1e-12
    end
    # Every requested radial state in the compatible ^3S_1/^3D_1 sectors enters
    # one tensor block; trace and the shared eigensystem are conserved.
    mixed = [s for s in spec.states if !isempty(s.mixings)]
    @test length(mixed) == 3
    @test all(m.mechanism == "tensor_mixing" for s in mixed for m in s.mixings)
    @test sum(s.mass_GeV for s in mixed) ≈
          sum(s.mixings[end].unmixed_GeV for s in mixed) atol = 1e-10
    @test all(
        s.mixings[end].partner_masses_GeV == mixed[1].mixings[end].partner_masses_GeV for
        s in mixed
    )
    @test all(s.mixings[end].result === mixed[1].mixings[end].result for s in mixed)
    @test fieldnames(StateMixing) == (:result, :eigenstate, :unmixed_GeV)
    @test all(length(physical_components(spec, s)) == 3 for s in mixed)
    components = physical_components(spec, mixed[1])
    coherent_r2 = sum(
        left.coefficient * right.coefficient *
        radial_overlap(left.wave, right.wave, r -> r^2) for
        left in components for right in components if
        left.basis.L_label == right.basis.L_label &&
        left.basis.multiplicity == right.basis.multiplicity &&
        left.basis.J == right.basis.J
    )
    @test radial_expect(spec, mixed[1], r -> r^2) ≈ coherent_r2 atol = 1e-12
    @test_throws ArgumentError radial_wave(spec, mixed[1])
    # A consumer that is specific to one spectroscopic channel receives the
    # coherent radial-n projection of the physical state, not one arbitrarily
    # selected precursor wave.  The remaining D component does not contribute
    # to an S-wave annihilation kernel.
    annihilation_input = GIModel.annihilation_basis_input(
        spec, BasisState(1, "S", 3, 1),
    )
    @test length(annihilation_input.radial_components) == 2
    @test all(isfinite(first(component)) for component in annihilation_input.radial_components)
    @test isfinite(GIModel._sL_smearing_factor(
        MomentumIntegralSmearing(120), annihilation_input, 0,
    ))
    # lookup by quantum numbers
    s = spectrum_state(spec, 1, "P", 3, 2)
    @test s.label == "1^3P_2"
    @test_throws ArgumentError spectrum_state(spec, 3, "S", 1, 0)
    # level beyond the per-channel budget fails loudly
    @test_throws ArgumentError compute_spectrum(
        params, meson;
        levels = [BasisState(7, "S", 1, 0)],
        solver = FiniteDifferenceSolver(ngrid = 80, rmax = 8.0),
    )
end

@testset "compute_spectrum same-J mixing gated by flavor content" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    pair = [BasisState(1, "P", 1, 1), BasisState(1, "P", 3, 1)]
    us = compute_spectrum(
        params, Meson(mq, :u, :s);
        levels = pair,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
    )
    so_mixed = [s for s in us.states if any(m.mechanism == "antisymmetric_spin_orbit" for m in s.mixings)]
    @test length(so_mixed) == 2
    @test all(s.fine_structure_mass_convention == "unequal_mass_same_j_mixed" for s in so_mixed)
    @test sum(s.mass_GeV for s in so_mixed) ≈
          sum(s.mixings[end].unmixed_GeV for s in so_mixed) atol = 1e-10
    @test any(abs(s.mixings[end].offdiag_GeV) > 0 for s in so_mixed)
    # equal flavor: the antisymmetric matrix element vanishes, no block forms
    cc = compute_spectrum(
        params, Meson(mq, :c, :c);
        levels = pair,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
    )
    @test all(isempty(s.mixings) for s in cc.states)
end

@testset "annihilation blocks built from two spectra (no reference data)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    ps_levels = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
    nn = compute_spectrum(
        params, Meson(mq, :q, :q);
        levels = ps_levels,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )
    ss = compute_spectrum(
        params, Meson(mq, :s, :s);
        levels = ps_levels,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )

    # the calibrated rank-one block reproduces its targets by construction
    targets = (0.520, 0.960, 1.440, 1.630)
    cal = pseudoscalar_annihilation_block(CalibratedP1Annihilation(), params, nn, ss;
        targets = targets)
    @test sort(cal.masses) ≈ sort(collect(targets)) atol = 1e-10
    calibrated_diagonal = [
        spectrum_state(nn, ps_levels[1]).mass_GeV,
        spectrum_state(ss, ps_levels[1]).mass_GeV,
        spectrum_state(nn, ps_levels[2]).mass_GeV,
        spectrum_state(ss, ps_levels[2]).mass_GeV,
    ]
    calibrated_update = cal.block.matrix - Diagonal(calibrated_diagonal)
    @test all(diag(calibrated_update) .>= 0.0)
    # explicit targets are mandatory - the digitized values live in GIPaper
    @test_throws ArgumentError pseudoscalar_annihilation_block(
        CalibratedP1Annihilation(), params, nn, ss)

    p1 = pseudoscalar_annihilation_block(PaperP1Annihilation(), params, nn, ss)
    p2 = pseudoscalar_annihilation_block(PaperP2Annihilation(), params, nn, ss)
    @test cal isa MixingResult
    @test p1 isa MixingResult
    @test p2 isa MixingResult
    @test isempty(p1.pole_matrices)
    @test length(p2.pole_matrices) == length(p2.masses)
    @test all(isfinite, p1.masses)
    @test all(isfinite, p2.masses)
    @test p1.masses != p2.masses
    pseudoscalar_basis = [
        annihilation_basis_input(nn, ps_levels[1]),
        annihilation_basis_input(ss, ps_levels[1]),
        annihilation_basis_input(nn, ps_levels[2]),
        annihilation_basis_input(ss, ps_levels[2]),
    ]
    @test_throws ArgumentError GIModel.isoscalar_pseudoscalar_annihilation_solution(
        PaperP2Annihilation(), params, pseudoscalar_basis; maxiter = 0,
    )
    @test_throws ArgumentError GIModel.isoscalar_pseudoscalar_annihilation_solution(
        PaperP2Annihilation(), params, pseudoscalar_basis; tol = 0.0,
    )
    # light nn̄ inputs carry the sqrt(2) flavor-coherence label
    input = annihilation_basis_input(nn, ps_levels[1])
    @test input.label == "1 ns"
    @test annihilation_basis_input(ss, ps_levels[2]).label == "2 ss"

    # general Eq. (16) block for one channel across the two flavors
    s1 = BasisState(1, "S", 3, 1)
    nn3 = compute_spectrum(
        params, Meson(mq, :q, :q);
        levels = [s1],
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )
    ss3 = compute_spectrum(
        params, Meson(mq, :s, :s);
        levels = [s1],
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )
    block = isoscalar_annihilation_block(params, nn3, ss3, s1;
        amplitude_A = params.annihilation.s1_A)
    @test length(block.masses) == 2
    @test all(isfinite, block.masses)
    # trace conservation: eigenvalue sum equals diagonal sum plus block trace
    diag_sum = spectrum_state(nn3, s1).mass_GeV + spectrum_state(ss3, s1).mass_GeV
    annihilation_matrix = block.block.matrix - Diagonal([
        spectrum_state(nn3, s1).mass_GeV,
        spectrum_state(ss3, s1).mass_GeV,
    ])
    @test sum(block.masses) ≈ diag_sum + sum(diag(annihilation_matrix)) atol = 1e-10
end

@testset "PA-16 final isoscalar spectrum owns flavor annihilation" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    solver = FiniteDifferenceSolver(ngrid = 80, rmax = 8.0, nlevels_per_channel = 2)
    levels = [
        BasisState(1, "S", 1, 0),
        BasisState(2, "S", 1, 0),
        BasisState(1, "S", 3, 1),
        BasisState(2, "S", 3, 1),
        BasisState(1, "D", 3, 1),
    ]
    nn = compute_spectrum(params, Meson(mq, :q, :q); levels = levels, solver = solver)
    ss = compute_spectrum(params, Meson(mq, :s, :s); levels = levels, solver = solver)
    final = add_isoscalar_annihilation(
        params, nn, ss;
        amplitudes = Dict(("S", 3, 1) => params.annihilation.s1_A),
    )

    @test final isa MixedSpectrum
    @test [(m.flavor1, m.flavor2) for m in final.channels] == [(:q, :q), (:s, :s)]
    @test length(final.states) == length(nn.states) + length(ss.states)
    @test all(!isnothing(state.basis.flavors) for state in final.states)
    @test_throws ArgumentError spectrum_state(final, 1, "S", 3, 1)

    qS = spectrum_state(final, BasisState(1, "S", 3, 1; flavors = (:q, :q)))
    sS = spectrum_state(final, BasisState(1, "S", 3, 1; flavors = (:s, :s)))
    q_ann = last(qS.mixings)
    s_ann = last(sS.mixings)
    @test q_ann.result === s_ann.result
    @test q_ann.mechanism == "general_eq16_annihilation"
    @test qS.mass_GeV == q_ann.result.masses[q_ann.eigenstate]
    @test sS.mass_GeV == s_ann.result.masses[s_ann.eigenstate]
    @test_throws ArgumentError radial_wave(final, qS)

    # Tensor S/D mixing precedes flavor mixing. The final state must flatten
    # both transformations, preserve flavor identity, and remain normalized.
    components = physical_components(final, qS)
    @test Set(component.basis.flavors for component in components) ==
          Set([(:q, :q), (:s, :s)])
    @test Set(component.basis.L_label for component in components) == Set(["S", "D"])
    @test sum(abs2(component.coefficient) for component in components) ≈ 1.0 atol = 1e-12
    @test radial_expect(final, qS, _ -> 1.0) ≈ 1.0 atol = 1e-12

    # The four-state radial/flavor block uses the same final-spectrum and
    # composition contract: unique native identities and normalized weights.
    p1_final = add_isoscalar_annihilation(
        params, nn, ss; pseudoscalar = PaperP1Annihilation(),
    )
    p1_state = first(sort(
        [
            state for state in p1_final.states if
            any(m -> m.mechanism == "paper_p1_pseudoscalar_annihilation", state.mixings)
        ];
        by = state -> state.mass_GeV,
    ))
    p1_components = physical_components(p1_final, p1_state)
    p1_identities = [
        (component.basis.n, component.basis.L_label, component.basis.multiplicity,
         component.basis.J, component.basis.flavors) for component in p1_components
    ]
    @test length(p1_identities) == length(unique(p1_identities))
    @test sum(abs2(component.coefficient) for component in p1_components) ≈ 1.0 atol = 1e-12
    weights = Dict(
        component.basis => float(i) for (i, component) in enumerate(p1_components)
    )
    @test physical_state_amplitude(p1_final, p1_state) do component
        weights[component.basis]
    end ≈ sum(
        component.coefficient * weights[component.basis] for
        component in p1_components
    )
    @test physical_transition_amplitude(p1_final, p1_state, p1_state) do left, right
        left.basis.flavors == right.basis.flavors ?
        radial_overlap(left.wave, right.wave, _ -> 1.0) : 0.0
    end ≈ 1.0 atol = 1e-10

    # No reference ordering or implicit amplitude exists in the model stage.
    combined_only = add_isoscalar_annihilation(params, nn, ss)
    @test all(
        mixing -> !occursin("annihilation", mixing.mechanism),
        Iterators.flatten(state.mixings for state in combined_only.states),
    )
    mismatched_ss = compute_spectrum(
        params, Meson(mq, :s, :s);
        levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 90, rmax = 8.0, nlevels_per_channel = 2),
    )
    @test_throws ArgumentError add_isoscalar_annihilation(params, nn, mismatched_ss)
    @test_throws ArgumentError add_isoscalar_annihilation(
        params, nn, ss; pseudoscalar_targets = [0.5, 1.0, 1.5, 2.0],
    )
    @test_throws ArgumentError compute_spectrum(
        params, Meson(mq, :q, :q);
        levels = [
            BasisState(1, "S", 3, 1),
            BasisState(1, "S", 3, 1; label = "duplicate"),
        ],
        solver = solver,
    )

    # Architecture boundary: core model code must never learn paper-row types,
    # target masses, or comparison assignment.
    core_bridge = read(joinpath(root, "src", "flavor_mixing.jl"), String)
    @test !occursin("ReferenceState", core_bridge)
    @test !occursin("reference_GeV", core_bridge)
    @test !occursin("compare_reference", core_bridge)

    # The general builder retains the actual requested orbital channel in its
    # shared basis; it must never relabel a P-wave block as an S wave.
    p_level = BasisState(1, "P", 3, 2)
    nnP = compute_spectrum(params, Meson(mq, :q, :q); levels = [p_level], solver = solver)
    ssP = compute_spectrum(params, Meson(mq, :s, :s); levels = [p_level], solver = solver)
    p_block = isoscalar_annihilation_block(
        params, nnP, ssP, p_level; amplitude_A = params.annihilation.a_3p2,
    )
    @test all(basis.L_label == "P" for basis in p_block.block.basis)
    @test all(basis.flavors in ((:q, :q), (:s, :s)) for basis in p_block.block.basis)
end

@testset "Annihilation phase convention, on the spectrum's own waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    lv = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
    solver = FiniteDifferenceSolver(ngrid = 80, rmax = 8.0)
    spec = compute_spectrum(params, Meson(mq, :q, :q); levels = lv, solver = solver)

    # A spectrum holds exactly one wave per level, from its own solver. The
    # second, oscillator-basis cache that used to sit beside it existed to work
    # around a normalization bug (FD Euclidean vs HO physical, a factor
    # 1/sqrt(h)); with that fixed the two solvers agree and the cache is gone.
    @test fieldnames(GIModel.SectorComputation) == (:params, :solver, :channel_cache)

    # The phase convention is NOT part of that fossil and must survive: the
    # eigensolver returns arbitrary column signs, and Table III amplitude signs
    # depend on Phi(0) proportional to the integral of r*u(r) being positive.
    for n in 1:2
        input = GIModel.annihilation_basis_input(spec, BasisState(n, "S", 1, 0))
        coefficient, w = only(input.radial_components)
        @test coefficient == 1.0
        @test sum(w.r .* w.u) * w.h > 0
    end

    # Applied without mutating the cached wave. Channel solutions already have
    # a deterministic outer-lobe phase; the annihilation convention may choose
    # the same or opposite sign for a particular radial excitation.
    key = RadialChannelKey(only(spec.channels).constituent_masses, "S", 1, 0)
    cached = spec.computation.channel_cache[key]
    before = [copy(w.u) for w in cached.waves]
    raw = radial_wave(spec, "1^1S_0")
    fixed = only(GIModel.annihilation_basis_input(
        spec, BasisState(1, "S", 1, 0),
    ).radial_components)[2]
    @test abs.(raw.u) ≈ abs.(fixed.u)
    @test sum(fixed.r .* fixed.u) > 0
    @test [w.u for w in cached.waves] == before

    # Both solvers give the same smeared origin factor -- the quantity the second
    # cache was introduced to correct. Agreement here is what makes it removable.
    m = mq["q"]
    sfd = wavefunction_origin_smearing(
        radial_wave(compute_spectrum(params, Meson(mq, :q, :q); levels = lv,
            solver = FiniteDifferenceSolver()), "1^1S_0"), m; L = 0)
    sho = wavefunction_origin_smearing(
        radial_wave(compute_spectrum(params, Meson(mq, :q, :q); levels = lv,
            solver = OscillatorSolver()), "1^1S_0"), m; L = 0)
    @test isapprox(abs(sfd), abs(sho); rtol = 0.01)
end

@testset "central diagnostic is independent of fixed -> mixed production" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    us = Meson(mq, :u, :s)
    levels = spectrum_levels(2; L_labels = ("S", "P"))
    kwargs = (solver = FiniteDifferenceSolver(ngrid = 250, rmax = 16.0),)

    central = central_spectrum(params, us; levels = levels, kwargs...)
    @test central isa CentralSpectrum
    @test eltype(central.states) === CentralState
    @test parameters(central) === central.computation.params
    @test all(isfinite(s.central_GeV) for s in central.states)

    corrected = fixed_spectrum(params, us; levels = levels, kwargs...)
    @test corrected isa CorrectedSpectrum
    @test corrected.computation !== central.computation
    @test all(k.multiplicity == 0 for k in keys(central.computation.channel_cache))
    @test all(k.multiplicity != 0 for k in keys(corrected.computation.channel_cache))
    for s in corrected.states
        @test s.basis isa BasisState
        @test !hasfield(CorrectedState, :central)
        @test s.mass_GeV ≈ s.central_GeV + s.contact_shift_GeV + s.fine_structure_shift_GeV
    end

    mixed = add_intra_meson_mixing(corrected)
    @test mixed isa MixedSpectrum
    @test mixed.computation === corrected.computation

    # compute_spectrum composes only the production stages
    direct = compute_spectrum(params, us; levels = levels, kwargs...)
    @test [s.mass_GeV for s in direct.states] == [s.mass_GeV for s in mixed.states]
    @test [s.fine_structure_mass_convention for s in direct.states] ==
          [s.fine_structure_mass_convention for s in mixed.states]

    # mixing overrides the stage view, never the corrected provenance
    so_mixed = [
        s for s in mixed.states if
        any(m.mechanism == "antisymmetric_spin_orbit" for m in s.mixings)
    ]
    @test !isempty(so_mixed)
    for s in so_mixed
        @test s.fine_structure_mass_convention == "unequal_mass_same_j_mixed"
        @test s.corrected.fine_structure_mass_convention == "unequal_mass_equal_share_LdotS"
        @test s.mixings[end].unmixed_GeV == s.corrected.mass_GeV
    end

    # stage skipping: no corrections means no fine structure and no mixing blocks
    bare = add_intra_meson_mixing(
        fixed_spectrum(
            params,
            us;
            levels = levels,
            kwargs...,
            terms = SpinTerms(contact_hyperfine = false, fine_structure = false),
        ),
    )
    @test all(s.mass_GeV == s.central_GeV for s in bare.states)
    @test all(isempty(s.mixings) for s in bare.states)
    @test all(s.fine_structure_mass_convention == "disabled" for s in bare.states)

    # stage order is enforced by dispatch
    @test !isdefined(GIModel, :add_spin_corrections)
    @test !hasmethod(add_intra_meson_mixing, Tuple{CentralSpectrum})
end

@testset "nonperturbative contact states expose hyperfine-distinct S waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = Meson(mq, :q, :q).constituent_masses
    central = central_spectrum(params, Meson(mq, :q, :q); levels = spectrum_levels(1))
    central_wave = radial_wave(central.computation.channel_cache[RadialChannelKey(masses, "S")], 1)
    r = central_wave.r

    sol1 = contact_hyperfine_nonperturbative_states(params, masses, "S", 1, r, 2)
    sol3 = contact_hyperfine_nonperturbative_states(params, masses, "S", 3, r, 2)
    @test radial_wave(sol1, 1).r == radial_wave(sol3, 1).r == r
    @test !isempty(sol1.waves) && !isempty(sol3.waves)
    # levels agree with the energy-only accessor
    @test sol1.eigenvalues_GeV ≈ GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 1, r, 2)
    @test sol3.eigenvalues_GeV ≈ GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 3, r, 2)

    # ^1S_0 (pi) is more compact than ^3S_1 (rho): smaller <r^2>, lower energy
    @test radial_expect(radial_wave(sol1, 1), x -> x^2) <
          radial_expect(radial_wave(sol3, 1), x -> x^2)
    @test sol1.eigenvalues_GeV[1] < sol3.eigenvalues_GeV[1]

    # The oscillator basis resums the same operator in its own space and lands
    # on the same answer, so this wrapper is basis-generic. "Empty" is reserved
    # for the genuinely inactive cases (non-S wave, sandwich off); it is never
    # "this basis has no implementation", which is what used to make callers
    # substitute first-order PT and report 0.2842 GeV for the light 1S0.
    sol_ho = contact_hyperfine_nonperturbative_states(
        params, masses, "S", 1, 2; solver = OscillatorSolver())
    @test abs(sol_ho.eigenvalues_GeV[1] - sol1.eigenvalues_GeV[1]) < 1e-3
end

@testset "Table VII gluonic annihilation (Eq. 17 S_L)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mb = mq["b"]
    masses = Meson(mq, :b, :b).constituent_masses
    tomev(a) = abs(a) * sqrt(1000)   # GeV^1/2 -> MeV^1/2

    # bottomonium (most paper-faithful) S-wave: eta_b / Upsilon central wave
    sol = channel_solution(
        params, masses, 0;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = 900, rmax = 24.0),
    )
    wave0 = radial_wave(sol, 1)
    S0 = wavefunction_origin_smearing(wave0, mb; L = 0)
    a0 = GIModel.alpha_s_q(sol.eigenvalues_GeV[1])
    # zero-parameter amplitudes vs paper (eta_b -> 2g = 2.5, Upsilon -> 3g = 0.21)
    @test 0.85 < tomev(gluonic_annihilation_amplitude(:S0_2g, S0, a0, mb)) / 2.5  < 1.15
    @test 0.85 < tomev(gluonic_annihilation_amplitude(:S1_3g, S0, a0, mb)) / 0.21 < 1.25

    # width is amplitude squared, for every channel
    for ch in GLUONIC_CHANNELS
        @test gluonic_annihilation_width(ch, S0, a0, mb) ≈
              gluonic_annihilation_amplitude(ch, S0, a0, mb)^2
    end

    # P-wave chi_2b via S1 (paper chi_2b -> 2g = 0.35)
    solP = channel_solution(
        params, masses, 1;
        nlevels = 1,
        solver = FiniteDifferenceSolver(ngrid = 900, rmax = 24.0),
    )
    S1 = wavefunction_origin_smearing(radial_wave(solP, 1), mb; L = 1)
    @test 0.8 < tomev(gluonic_annihilation_amplitude(:P2_2g, S1, GIModel.alpha_s_q(solP.eigenvalues_GeV[1]), mb)) / 0.35 < 1.2

    # S_L is normalization-invariant (the wave is renormalized internally)
    S0_scaled = wavefunction_origin_smearing(MeshWave(3.0 .* wave0.u, wave0.r), mb; L = 0)
    @test S0_scaled ≈ S0
    # unknown channel is rejected
    @test_throws ArgumentError gluonic_annihilation_amplitude(:bogus, S0, a0, mb)
end

@testset "Table VII leptonic decay constants (mock-meson factors)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    outer!(u) = (pk = maximum(abs, u); i = findlast(x -> abs(x) > 0.2pk, u);
                 (i !== nothing && u[i] < 0) && (u .*= -1); u)

    # triplet ³S₁ wave for a QQ̄; return (M, unit-phase wave)
    function triplet_swave(m, n)
        masses = ConstituentMasses(m, m)
        sol = contact_hyperfine_nonperturbative_states(
            params, masses, "S", 3, 2;
            solver = FiniteDifferenceSolver(ngrid = 1000, rmax = 24.0),
        )
        wave = radial_wave(sol, n)
        return sol.eigenvalues_GeV[n], MeshWave(outer!(copy(wave.u)), wave.r)
    end

    # ψ -> e+e-: f = (16/3)^(1/2) V_ψ, paper 0.12
    Mψ, wψ = triplet_swave(mq["c"], 1)
    fψ = sqrt(16 / 3) * leptonic_decay_factor(:V_V, wψ, mq["c"], mq["c"], Mψ)
    @test 0.85 < fψ / 0.12 < 1.15

    # ρ -> e+e-: f = √6 V_ρ, paper 0.20 (use light triplet 1S)
    Mρ, wρ = triplet_swave(mq["q"], 1)
    fρ = sqrt(6) * leptonic_decay_factor(:V_V, wρ, mq["q"], mq["q"], Mρ)
    @test 0.85 < fρ / 0.20 < 1.15

    # mock mass is the free-pair energy: M̃ = <E1+E2> ≥ 2m, and finite
    Mt = mock_meson_mass(wψ, mq["c"], mq["c"]; L = 0)
    @test Mt > 2 * mq["c"]
    @test Mt < 2 * mq["c"] + 2.0

    # kind guards: D-/P-wave factors require equal masses; unknown kind rejected
    @test_throws ArgumentError leptonic_decay_factor(:Vp_V, wψ, mq["c"], mq["b"], Mψ)
    @test_throws ArgumentError leptonic_decay_factor(:bogus, wψ, mq["c"], mq["c"], Mψ)
end

@testset "Decay widths from decay constants (Eqs. D7-D9)" begin
    hbar = 6.582119e-25   # GeV·s

    # D7: Γ(P→ℓν) reproduces the measured π→μν width when fed the EXPERIMENTAL
    # f_π/M_π and physical masses (the ~4% excess is the Cabibbo cos²θ_C the
    # paper's reduced G² omits).
    f_pi = 0.1307 / 0.13957
    Γ_pi = leptonic_pseudoscalar_width(f_pi, 0.13957, 0.10566)
    Γ_pi_pdg = hbar / 2.6033e-8
    @test 0.98 < Γ_pi / Γ_pi_pdg < 1.08
    @test leptonic_pseudoscalar_width(f_pi, 0.10, 0.10566) == 0.0   # m_ℓ ≥ M_P: closed
    @test leptonic_pseudoscalar_width(0.0, 0.13957, 0.10566) == 0.0 # f=0

    # D8: Γ(V→ℓ⁺ℓ⁻) round-trips the measured ψ→ee width; scales as M·f².
    f_V = sqrt(5.55e-6 / ((4π / 3) * GIModel.ALPHA_EM^2 * 3.0969))
    @test isapprox(dilepton_vector_width(f_V, 3.0969), 5.55e-6; rtol = 1e-6)
    @test isapprox(dilepton_vector_width(2f_V, 3.0969), 4 * 5.55e-6; rtol = 1e-6)  # ∝ f²

    # D9: analytic value + kinematic bracket (closes as M_A1 → m_τ).
    @test isapprox(axial_tau_width(0.15, 1.26, 1.777), 5.412e-13; rtol = 5e-3)
    @test axial_tau_width(0.15, 1.777, 1.777) == 0.0          # M_A1 = m_τ: closed
    @test axial_tau_width(0.15, 1.80, 1.777) == 0.0           # M_A1 > m_τ: closed
    @test axial_tau_width(0.15, 1.26, 1.777) > 0
end

@testset "Appendix-D mock-meson overlap kernels (D2-D3)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["q"]
    masses = ConstituentMasses(m, m)
    H, r = GIModel.relativistic_hamiltonian(params, masses, 0; ngrid = 450, rmax = 24.0)
    # hyperfine-distinct nn 1S waves: ¹S₀ (pi), ³S₁ (rho); physical ∫u²dr=1
    op1 = GIModel.contact_hyperfine_operator(params, masses, "S", 1, r)
    op3 = GIModel.contact_hyperfine_operator(params, masses, "S", 3, r)
    _, v1 = GIModel.lowest_eigenpairs(Symmetric(Matrix(H) + Matrix(op1)), 2)
    _, v3 = GIModel.lowest_eigenpairs(Symmetric(Matrix(H) + Matrix(op3)), 2)
    h = r[2] - r[1]
    phase!(u) = (sum(r .* u) < 0 && (u .*= -1); u)
    u_pi = phase!(v1[:, 1] ./ sqrt(h)); u_rho = phase!(v3[:, 1] ./ sqrt(h))
    w_pi = MeshWave(u_pi, r); w_rho = MeshWave(u_rho, r)

    mw_pi = mock_momentum_wave(w_pi, 0); mw_rho = mock_momentum_wave(w_rho, 0)
    # momentum wave normalized ∫p²Φ²dp=1; mock mass ≥ 2m and finite
    @test isapprox(sum(0.5 * (mw_pi.p[2:end] .^ 2 .* mw_pi.phi[2:end] .^ 2 .+
                              mw_pi.p[1:end-1] .^ 2 .* mw_pi.phi[1:end-1] .^ 2) .*
                        diff(mw_pi.p)), 1.0; atol = 1e-3)
    @test mock_mean_energy(mw_pi, m) ≥ m
    Mpi = mock_wave_mass(mw_pi, m, m); Mrho = mock_wave_mass(mw_rho, m, m)
    @test 2m ≤ Mpi ≤ 2m + 2.0

    # I_i drives ρ→πγ, the paper's 0.7-exponent fit row: μ = (1/3) I_ρπ M_N ≈ +0.69
    M_N = 0.93827
    μ_rho = (1 / 3) * mock_meson_overlap(mw_pi, mw_rho, m; Mx = Mpi, My = Mrho) * M_N
    @test 0.60 < μ_rho < 0.72        # audit/report value +0.650, paper +0.69

    # Eₙⁱ radial moment: mesh guard + the n=1 self-moment recovers ⟨r⟩ scaling
    @test_throws ArgumentError mock_meson_radial_moment(
        w_pi, MeshWave(u_pi, r .+ 1.0), 1.0, 1.0, m)
    E1 = mock_meson_radial_moment(w_pi, w_pi, m, m, m; n = 1, exponent = 0.0)
    @test isapprox(E1, sum(@. u_pi^2 * r) * h; rtol = 1e-9)   # exponent 0 ⇒ ∫u²r dr
    @test mock_meson_radial_moment(w_pi, w_pi, m, m, m; n = 1) > 0

    # --- Eq. (22) assembly, promoted out of the Table VI audit ---------------
    # Photon momentum q = (M² - M'²)/2M, and q → 0 as the masses close up.
    @test isapprox(photon_momentum(3.686, 2.980), (3.686^2 - 2.980^2) / (2 * 3.686);
        rtol = 1e-12)
    @test photon_momentum(1.0, 1.0) == 0.0

    # m1_transition_moment is Σ c·I_i·M_N: the ρ→πγ fit row rebuilt through the
    # public assembly must equal the hand-written kernel product above.
    @test isapprox(m1_transition_moment(mw_pi, mw_rho, m, m, [(1 / 3, m)]), μ_rho;
        rtol = 1e-12)
    @test NUCLEON_MASS_GEV == M_N
    # Coefficients are linear and the antiquark charge enters flipped, so the
    # neutral combination (+2/3, -1/3) is the sum of its two single-line terms.
    @test isapprox(
        m1_transition_moment(mw_pi, mw_rho, m, m, [(2 / 3, m), (-1 / 3, m)]),
        m1_transition_moment(mw_pi, mw_rho, m, m, [(1 / 3, m)]); rtol = 1e-12)

    # e1_transition_amplitude = coeff(q)·E₁ⁱ·√(α q_MeV); an explicit q overrides
    # the mass-implied one, and the amplitude is linear in the row coefficient.
    amp = e1_transition_amplitude(w_pi, mw_pi, w_rho, mw_rho, m, q -> 4q / 9, 1.5, 1.0)
    q_implied = photon_momentum(1.5, 1.0)
    E1i = mock_meson_radial_moment(w_pi, w_rho, mock_mean_energy(mw_pi, m),
        mock_mean_energy(mw_rho, m), m; n = 1)
    @test isapprox(amp, (4 * q_implied / 9) * E1i * sqrt(ALPHA_EM * 1000 * q_implied);
        rtol = 1e-12)
    @test isapprox(
        e1_transition_amplitude(w_pi, mw_pi, w_rho, mw_rho, m, q -> 4q / 9, 9.9, 9.9;
            q = q_implied), amp; rtol = 1e-12)
    @test isapprox(
        e1_transition_amplitude(w_pi, mw_pi, w_rho, mw_rho, m, q -> 8q / 9, 1.5, 1.0),
        2amp; rtol = 1e-12)
end

@testset "Table VII two-photon amplitudes (part b)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    outer!(u) = (pk = maximum(abs, u); i = findlast(x -> abs(x) > 0.2pk, u);
                 (i !== nothing && u[i] < 0) && (u .*= -1); u)
    tokeV(a) = a * sqrt(1e6)

    # eta_c -> gamma gamma: ¹S₀ cc̄, q_eff = 4/9, paper 2.6 keV^½
    mc = mq["c"]
    solηc = contact_hyperfine_nonperturbative_states(
        params, ConstituentMasses(mc, mc), "S", 1, 2;
        solver = FiniteDifferenceSolver(ngrid = 1000, rmax = 24.0),
    )
    rawηc = radial_wave(solηc, 1)
    wηc = MeshWave(outer!(copy(rawηc.u)), rawηc.r)
    Aηc = tokeV(two_photon_amplitude(:P, wηc, mc, solηc.eigenvalues_GeV[1], 4 / 9))
    @test 0.85 < Aηc / 2.6 < 1.25

    # A2 -> gamma gamma: ³P₂ light isovector, q_eff = (e_u²−e_d²)/√2, paper −1.2 keV^½
    mqk = mq["q"]
    solP = channel_solution(
        params, ConstituentMasses(mqk, mqk), 1;
        nlevels = 1,
        solver = FiniteDifferenceSolver(ngrid = 1000, rmax = 24.0),
    )
    rawA2 = radial_wave(solP, 1)
    wA2 = MeshWave(outer!(copy(rawA2.u)), rawA2.r)
    AA2 = tokeV(two_photon_amplitude(:P2, wA2, mqk, solP.eigenvalues_GeV[1], (4 / 9 - 1 / 9) / sqrt(2)))
    @test AA2 < 0                       # −√(4/5) prefactor => negative amplitude
    @test 0.8 < abs(AA2) / 1.2 < 1.2

    @test_throws ArgumentError two_photon_amplitude(:bogus, wηc, mc, solηc.eigenvalues_GeV[1], 4 / 9)
end

@testset "Table VII charge radii (part d)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mu, ms = mq["q"], mq["s"]

    function ps_wave(m1, m2)
        sol = contact_hyperfine_nonperturbative_states(
            params, ConstituentMasses(m1, m2), "S", 1, 2;
            solver = FiniteDifferenceSolver(ngrid = 1200, rmax = 26.0),
        )
        return radial_wave(sol, 1)
    end

    # K+ = u s̄ : charges +2/3, +1/3 ; paper r_E² = +(0.59)² fm²  (a prediction)
    rEK = charge_radius_squared(ps_wave(mu, ms), mu, 2 / 3, ms, 1 / 3) * HBARC_FM2
    @test rEK > 0
    @test 0.85 < rEK / 0.59^2 < 1.15

    # K0 = d s̄ : charges −1/3, +1/3 ; paper r_E² = −(0.30)² fm²  (sign is the point)
    rEK0 = charge_radius_squared(ps_wave(mu, ms), mu, -1 / 3, ms, 1 / 3) * HBARC_FM2
    @test rEK0 < 0                       # negative charge on the larger-radius light quark
    @test 0.8 < abs(rEK0) / 0.30^2 < 1.3

    # normalization-invariant (wave renormalized internally)
    w = ps_wave(mu, mu)
    @test charge_radius_squared(MeshWave(3.0 .* w.u, w.r), mu, 2 / 3, mu, 1 / 3) ≈
          charge_radius_squared(w, mu, 2 / 3, mu, 1 / 3)
end

@testset "Table VII mixed eta/eta' two-photon (P1 coherent sum)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mu, ms = mq["q"], mq["s"]
    Qnn = (4 / 9 + 1 / 9) / sqrt(2); Qss = 1 / 9
    psl = [GIModel.BasisState(1, "S", 1, 0), GIModel.BasisState(2, "S", 1, 0)]
    final = compute_isoscalar_spectrum(
        params,
        Meson(LightQuark(mu), LightQuark(mu)),
        Meson(StrangeQuark(ms), StrangeQuark(ms));
        levels = psl,
        solver = FiniteDifferenceSolver(ngrid = 1000, rmax = 24.0),
        pseudoscalar = PaperP1Annihilation(),
    )
    states = sort(final.states; by = state -> state.mass_GeV)
    function ggamp(state, Mphys)
        return physical_state_amplitude(final, state) do component
            mass, charge = component.basis.flavors == (:q, :q) ? (mu, Qnn) : (ms, Qss)
            two_photon_amplitude(:P, component.wave, mass, Mphys, charge)
        end * sqrt(1e6)
    end
    Aη = ggamp(states[1], 0.548)
    Aη′ = ggamp(states[2], 0.958)
    # both positive, and the η<η' ordering that IDEAL mixing gets backwards
    @test Aη > 0 && Aη′ > 0
    @test abs(Aη) < abs(Aη′)
    # magnitudes in the paper's ballpark (η≈0.5, η'≈1.3 keV^½)
    @test 0.6 < abs(Aη) / 0.5 < 1.6
    @test 0.6 < abs(Aη′) / 1.3 < 1.4
end

@testset "Native HO fixed-channel spin-distorted waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    solver_ho = OscillatorSolver()
    mc = mq["c"]
    masses = ConstituentMasses(mc, mc)
    ngrid, rmax = 900, 24.0
    tomev(a) = abs(a) * sqrt(1000)
    sm(w::RadialWave, L) = wavefunction_origin_smearing(w, mc; L = L)

    # 1. basis-fidelity control: HO central S_L matches FD to <2% for charm —
    #    the 15-20% gluonic row residuals were never a basis artifact
    fdsol = channel_solution(
        params, masses, 0;
        nlevels = 2,
        solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax),
    )
    hosol = channel_solution(params, masses, 0; solver = solver_ho, nlevels = 2)
    @test 0.98 < abs(wavefunction_origin_smearing(radial_wave(hosol, 1), mc; L = 0)) /
                 sm(radial_wave(fdsol, 1), 0) < 1.02

    # 2. Paper-order full diagonalization with literal A15-A16 spin blocks
    #    lands the charm gluonic rows on the paper.
    amp(ch, S, M) = tomev(gluonic_annihilation_amplitude(ch, S, GIModel.alpha_s_q(M), mc))
    ratios = Dict{Symbol,Float64}()
    cases = (
        (:eta_c, FineStructureMultiplet("S", 1, 0), :S0_2g, 4.700),
        (:psi, FineStructureMultiplet("S", 3, 1), :S1_3g, 0.420),
        (:chi_0c, FineStructureMultiplet("P", 3, 0), :P0_2g, 2.500),
        (:chi_2c, FineStructureMultiplet("P", 3, 2), :P2_2g, 0.880),
    )
    for (key, multiplet, ch, paper) in cases
        sol = fixed_channel_solution(
            params, masses, multiplet; solver = solver_ho, nlevels = 4)
        L = GIModel.L_SYMBOLS[multiplet.L_label]
        ratios[key] = amp(ch, sm(radial_wave(sol, 1), L), sol.eigenvalues_GeV[1]) / paper
        @test 0.85 < ratios[key] < 1.15
    end

    # 3. the splitting patterns the central wave misses collapse at paper order:
    #    central waves give ratio-of-ratios eta_c/psi ≈ 0.78, chi_0c/chi_2c ≈ 0.68
    Sc0 = sm(radial_wave(fdsol, 1), 0)
    Mc0 = channel_solution(
        params, masses, 0;
        nlevels = 1,
        solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax),
    ).eigenvalues_GeV[1]
    central_eta_psi = (amp(:S0_2g, Sc0, Mc0) / 4.700) / (amp(:S1_3g, Sc0, Mc0) / 0.420)
    @test central_eta_psi < 0.85
    @test 0.90 < ratios[:eta_c] / ratios[:psi] < 1.10
    @test ratios[:chi_0c] / ratios[:chi_2c] > 0.80   # central-wave value ≈ 0.68
end

@testset "FD observable momentum cutoff is refinement-stable" begin
    # Refining an FD coordinate mesh must not silently coarsen the momentum
    # quadrature used by annihilation/leptonic observables. The production grid
    # already reaches the certified 60 GeV physical cutoff; finer grids retain
    # that range instead of extending pmax = π/h at fixed npoints.
    r_coarse = collect(range(0.05, 20.0; step = 0.05))
    r_fine = collect(range(0.025, 20.0; step = 0.025))
    trial(r) = MeshWave(r .* exp.(-r), r)
    coarse_momentum = GIModel._observable_momentum_wave(trial(r_coarse), 0, 900)
    fine_momentum = GIModel._observable_momentum_wave(trial(r_fine), 0, 900)
    @test last(coarse_momentum.p) ≈ min(π / 0.05, GIModel.FD_OBSERVABLE_PMAX_GEV)
    @test last(fine_momentum.p) == GIModel.FD_OBSERVABLE_PMAX_GEV
    @test length(coarse_momentum.p) == length(fine_momentum.p) == 900
    for mass in (0.22, 4.977)
        coarse = wavefunction_origin_smearing(trial(r_coarse), mass; L = 0)
        fine = wavefunction_origin_smearing(trial(r_fine), mass; L = 0)
        @test abs(coarse - fine) / abs(fine) < 1.0e-3
    end
end

@testset "Native HO full fixed-channel diagonalization" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    solver_ho = OscillatorSolver()
    ngrid, rmax = 900, 24.0

    # 1. The paper's treatment is FULL diagonalization in the finite HO basis,
    #    not first-order PT. The light ¹S₀ (pion) mass discriminates: full-diag
    #    keeps it resummed (≈0.15 GeV) like the fine-grid FD, while first-order
    #    PT over-raises it (≈0.28 GeV).
    nn = ConstituentMasses(mq["q"], mq["q"])
    pion = FineStructureMultiplet("S", 1, 0)
    certified_pion = fixed_channel_solution(
        params, nn, pion; solver = solver_ho, nlevels = 4)
    m_full = certified_pion.eigenvalues_GeV[1]
    m_fd = fixed_channel_solution(
        params, nn, pion;
        solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax),
        nlevels = 4,
    ).eigenvalues_GeV[1]
    @test 0.12 < m_full < 0.18        # resummed, near the paper's 0.15 GeV
    @test abs(m_full - m_fd) < 0.02   # matches the fine-grid FD resummation

    @test certified_pion.convergence.status == :converged
    @test certified_pion.convergence.energy_delta_GeV <=
          solver_ho.energy_tolerance_GeV

    # 2. full diagonalization still lands the charm gluonic singlet on the paper
    #    (the spin distortion the central wave misses), so both subtable regimes
    #    are served by ONE treatment.
    mc = mq["c"]; cc = ConstituentMasses(mc, mc)
    sol = fixed_channel_solution(
        params, cc, pion; solver = solver_ho, nlevels = 4)
    wave = radial_wave(sol, 1)
    S = wavefunction_origin_smearing(wave, mc; L = 0)
    eta_c = abs(gluonic_annihilation_amplitude(
        :S0_2g, S, GIModel.alpha_s_q(sol.eigenvalues_GeV[1]), mc)) * sqrt(1000) / 4.700
    @test 0.95 < eta_c < 1.20
end
