using Test

root = dirname(@__DIR__)
include(joinpath(root, "src", "GIModel", "GIModel.jl"))
using .GIModel

@testset "Table II digitization vs parameters TOML" begin
    script = joinpath(root, "scripts", "verify_table_ii_toml.py")
    p = run(`python3 $script`, wait = false)
    wait(p)
    @test success(p)
end

@testset "reference spectrum CSVs (required columns)" begin
    script = joinpath(root, "scripts", "validate_reference_spectra.py")
    p = run(`python3 $script`, wait = false)
    wait(p)
    @test success(p)
end

@testset "reference loading" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    @test params.masses["c"] ≈ 1.628
    @test params.masses["b"] ≈ 4.977
    @test params.b ≈ 0.18
    @test params.appendix_a_smearing == false
    @test params.epsilon_c ≈ -0.168
    @test params.fine_structure == true
    @test params.k_spin_orbit > 0.0

    ccbar = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    bbbar = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_bottomonium.csv"))
    @test length(ccbar) == 28
    @test length(bbbar) == 30
    @test ccbar[1].quark_content == "c cbar"
    @test ccbar[1].composition == "1^1S_0"
end

@testset "baseline solver shape" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    cc = solve_sector(params, "c"; maxn = 4, ngrid = 250, rmax = 20.0)
    bb = solve_sector(params, "b"; maxn = 4, ngrid = 250, rmax = 16.0)

    @test cc[(1, "S")] < cc[(2, "S")] < cc[(3, "S")]
    @test bb[(1, "S")] < bb[(2, "S")] < bb[(3, "S")]
    @test cc[(1, "S")] < cc[(1, "P")] < cc[(1, "D")]
    @test bb[(1, "S")] < bb[(1, "P")] < bb[(1, "D")]
end

@testset "reduced_mass" begin
    @test reduced_mass(1.5, 0.3) ≈ (1.5 * 0.3) / (1.5 + 0.3)
    @test reduced_mass(2.0, 2.0) ≈ 1.0
end

@testset "quark mass resolution" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m1, m2 = parse_quark_masses(params, "ccbar", "ignore")
    @test m1 ≈ m2 ≈ params.masses["c"]
    m1, m2 = parse_quark_masses(params, "charmonium", "c cbar")
    @test m1 ≈ m2 ≈ params.masses["c"]
    m1, m2 = parse_quark_masses(params, "charmed", "-c dbar; c ubar")
    @test m1 ≈ params.masses["c"] && m2 ≈ params.masses["d"]
    m1, m2 = parse_quark_masses(params, "charmed_strange", "c sbar")
    @test m1 ≈ params.masses["c"] && m2 ≈ params.masses["s"]
    m1, m2 = parse_quark_masses(params, "bottom_light", "b ubar; -b dbar")
    @test m1 ≈ params.masses["b"] && m2 ≈ params.masses["u"]
    m1, m2 = parse_quark_masses(params, "isoscalar", "ignore")
    @test m1 ≈ m2 ≈ 0.5 * (params.masses["u"] + params.masses["d"])
end

@testset "appendix_a_smearing code path (finite S-wave energy)" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        @test occursin("appendix_a_smearing = false", s)
        write(p, replace(s, "appendix_a_smearing = false" => "appendix_a_smearing = true"))
        params = load_parameters(p)
        @test params.appendix_a_smearing == true
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params, mc, mc, 0;
            nlevels = 2, ngrid = 120, rmax = 12.0, kinetic = :relativistic,
        )
        @test isfinite(vals[1]) && isfinite(vals[2])
        @test vals[1] < vals[2]
    end
end

@testset "Appendix A 3D smearing (constant preserves norm)" begin
    h = 0.02
    n = 2000
    r = collect(h:h:(h * n))
    σ = 1.0
    n_tail = max(0, Int(ceil(8 * σ / h)))
    r_ext = n_tail > 0 ? vcat(r, collect((r[end] + h):h:(r[end] + n_tail * h))) : r
    w = GIModel.smear_3d_radial(ones(length(r_ext)), r_ext, σ)
    w = w[1:n]
    @test maximum(abs.(w .- 1.0)) < 0.01
    @test w[1] ≈ 1.0 atol = 0.01
    @test w[div(n, 2)] ≈ 1.0 atol = 0.01
end

@testset "L·S and spin_dot algebra" begin
    @test GIModel.LdotS(1, 1, 0) ≈ -2.0
    @test GIModel.LdotS(1, 1, 1) ≈ -1.0
    @test GIModel.LdotS(1, 1, 2) ≈ 1.0
    @test GIModel.spin_dot(1) ≈ -0.75
    @test GIModel.spin_dot(3) ≈ 0.25
end

@testset "triplet fine-structure angular factors" begin
    for L in 1:4
        js = collect((L - 1):(L + 1))
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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _, umat, r = GIModel.channel_solution(params, m, m, 0; nlevels = 2, ngrid = 200, rmax = 20.0)
    h = r[2] - r[1]
    u_s = collect(umat[:, 1])
    @test GIModel.fine_structure_split(
        params, m, m, "S", 3, 1, u_s, r, h;
        k_spin_orbit = params.k_spin_orbit, k_tensor = params.k_tensor,
    ) == 0.0
    @test GIModel.fine_structure_split(
        params, m, m, "S", 1, 0, u_s, r, h;
        k_spin_orbit = 1.0, k_tensor = 1.0,
    ) == 0.0
    v_p, umat_p, r_p = GIModel.channel_solution(params, m, m, 1; nlevels = 2, ngrid = 200, rmax = 20.0)
    h_p = r_p[2] - r_p[1]
    u1p = collect(umat_p[:, 1])
    δ0 = GIModel.fine_structure_split(
        params, m, m, "P", 3, 0, u1p, r_p, h_p;
        k_spin_orbit = 1.0, k_tensor = 1.0,
    )
    δ1 = GIModel.fine_structure_split(
        params, m, m, "P", 3, 1, u1p, r_p, h_p;
        k_spin_orbit = 1.0, k_tensor = 1.0,
    )
    δ2 = GIModel.fine_structure_split(
        params, m, m, "P", 3, 2, u1p, r_p, h_p;
        k_spin_orbit = 1.0, k_tensor = 1.0,
    )
    @test isfinite(δ0) && isfinite(δ1) && isfinite(δ2)
    @test δ0 != δ1 || δ1 != δ2
    @test GIModel.fine_structure_split(
        params, m, m, "P", 1, 0, u1p, r_p, h_p;
        k_spin_orbit = 1.0, k_tensor = 1.0,
    ) == 0.0
end

@testset "erf_approx basic symmetries" begin
    @test GIModel.erf_approx(0.0) ≈ 0.0 atol = 1e-7
    for x in (0.05, 0.3, 0.8, 1.5)
        @test GIModel.erf_approx(x) + GIModel.erf_approx(-x) ≈ 0.0 atol = 1e-12
    end
end

@testset "Coulomb derivative consistency (erf_approx)" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    for r0 in (0.05, 0.2, 1.3)
        δ = 1e-6 * max(1.0, r0)
        num = (GIModel.alpha_s_r(r0 + δ) - GIModel.alpha_s_r(r0 - δ)) / (2δ)
        @test GIModel.alpha_s_prime_r(r0) ≈ num rtol = 1e-6 atol = 1e-10
    end

    V(r) = -(4 / 3) * GIModel.alpha_s_r(r) / r
    for r0 in (0.2, 1.0)
        δ = 1e-6 * max(1.0, r0)
        num = (V(r0 + δ) - V(r0 - δ)) / (2δ)
        @test GIModel.dV_coul_central_dr(r0, params) ≈ num rtol = 1e-6 atol = 1e-10
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

@testset "contact hyperfine: only S-waves" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    r = collect(0.05:0.05:1.0)
    u = exp.(-2.0 .* r)
    @test GIModel.contact_hyperfine_shift(params, m, m, "P", 3, u, r) == 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "D", 3, u, r) == 0.0
end

@testset "contact hyperfine shift uses u(r) normalization" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    r = collect(0.05:0.05:1.0)
    u = exp.(-2.0 .* r)
    base = GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u, r)
    @test base != 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, 3.0 .* u, r) ≈ base rtol = 1e-12 atol = 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, 0.2 .* u, r) ≈ base rtol = 1e-12 atol = 0.0
end

@testset "write_residual_report keyword alias (appendix_a_central)" begin
    rows = [
        (
            sector = "test",
            state = "1^1S_0",
            L = "S",
            n = 1,
            J = 0,
            multiplicity = 1,
            reference_GeV = 1.0,
            predicted_GeV = 1.0,
            residual_MeV = 0.0,
            confidence = "test-only",
        ),
    ]
    mktemp() do path, io
        close(io)
        GIModel.write_residual_report(path, "alias check", rows; appendix_a_central = true)
        txt = read(path, String)
        @test occursin("Appendix A 3D isotropic smearing", txt)
    end
end
