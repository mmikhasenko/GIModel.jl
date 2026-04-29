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
    @test params.coulomb_1d_smear == false
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

@testset "central: Coulomb+confinement = pointwise V" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    for r0 in (0.15, 0.4, 1.2, 3.0)
        a = GIModel.central_potential(r0, params)
        b = GIModel.static_coulomb_G(r0, params) + GIModel.static_confinement_S(r0, params)
        @test a ≈ b rtol = 1e-12
    end
end

@testset "central_potential_path (default = pointwise)" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    p = GIModel.central_potential_path(params)
    @test p.name == "pointwise_fd"
    @test params.appendix_a_smearing == false
    @test params.coulomb_1d_smear == false
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

@testset "coulomb_1d_smear code path (finite S-wave energy)" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        @test occursin("coulomb_1d_smear = false", s)
        write(p, replace(s, "coulomb_1d_smear = false" => "coulomb_1d_smear = true"))
        params = load_parameters(p)
        @test params.coulomb_1d_smear == true
        @test params.appendix_a_smearing == false
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params, mc, mc, 0;
            nlevels = 2, ngrid = 120, rmax = 12.0, kinetic = :relativistic,
        )
        @test isfinite(vals[1]) && isfinite(vals[2])
        @test vals[1] < vals[2]
    end
end

@testset "1D Coulomb smear: constant vector unchanged on uniform grid" begin
    h = 0.05
    n = 200
    r = collect(h:h:(h * n))
    c = GIModel.convolve_1d_gaussian_same_length(ones(n), r, h, 0.3)
    @test maximum(abs.(c .- 1.0)) < 1e-12
end

@testset "appendix_a_smearing wins over coulomb_1d in central_potential_path" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s2 = replace(s, "appendix_a_smearing = false" => "appendix_a_smearing = true")
        s2 = replace(s2, "coulomb_1d_smear = false" => "coulomb_1d_smear = true")
        write(p, s2)
        params = load_parameters(p)
        @test params.appendix_a_smearing == true
        @test params.coulomb_1d_smear == true
        path = GIModel.central_potential_path(params)
        @test path.name == "experimental_3d_convl_a7a8"
    end
end

@testset "Appendix A 3D smearing (constant preserves norm)" begin
    h = 0.02
    n = 2000
    r = collect(h:h:(h * n))
    σ = 1.0
    n_tail = σ > 0 ? max(0, Int(ceil(8 / (σ * h)))) : 0
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

@testset "fine_structure_components: decomposition sums correctly" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _v_p, umat_p, r_p = GIModel.channel_solution(params, m, m, 1; nlevels = 2, ngrid = 200, rmax = 20.0)
    h_p = r_p[2] - r_p[1]
    u1p = collect(umat_p[:, 1])
    for J in (0, 1, 2)
        comp = GIModel.fine_structure_components(
            params, m, m, "P", 3, J, u1p, r_p, h_p;
            k_spin_orbit = 1.0, k_tensor = 1.0,
        )
        @test comp.total ≈ comp.spin_orbit + comp.tensor atol = 1e-12
        @test comp.total ≈ GIModel.fine_structure_split(
            params, m, m, "P", 3, J, u1p, r_p, h_p;
            k_spin_orbit = 1.0, k_tensor = 1.0,
        ) atol = 1e-12
    end
end

@testset "fine structure uses u(r) normalization (scale invariant)" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _v_p, umat_p, r_p = GIModel.channel_solution(params, m, m, 1; nlevels = 2, ngrid = 200, rmax = 20.0)
    h_p = r_p[2] - r_p[1]
    u1p = collect(umat_p[:, 1])
    for J in (0, 1, 2)
        base = GIModel.fine_structure_components(
            params, m, m, "P", 3, J, u1p, r_p, h_p;
            k_spin_orbit = 1.0, k_tensor = 1.0,
        )
        scaled_hi = GIModel.fine_structure_components(
            params, m, m, "P", 3, J, 3.0 .* u1p, r_p, h_p;
            k_spin_orbit = 1.0, k_tensor = 1.0,
        )
        scaled_lo = GIModel.fine_structure_components(
            params, m, m, "P", 3, J, 0.2 .* u1p, r_p, h_p;
            k_spin_orbit = 1.0, k_tensor = 1.0,
        )
        @test scaled_hi.total ≈ base.total rtol = 1e-12 atol = 0.0
        @test scaled_lo.total ≈ base.total rtol = 1e-12 atol = 0.0
        @test scaled_hi.spin_orbit ≈ base.spin_orbit rtol = 1e-12 atol = 0.0
        @test scaled_lo.spin_orbit ≈ base.spin_orbit rtol = 1e-12 atol = 0.0
        @test scaled_hi.tensor ≈ base.tensor rtol = 1e-12 atol = 0.0
        @test scaled_lo.tensor ≈ base.tensor rtol = 1e-12 atol = 0.0
    end
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

@testset "contact hyperfine matches radial_expect_udr convention" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    r = collect(0.05:0.05:1.0)
    h = r[2] - r[1]
    u = exp.(-2.0 .* r)
    sigma = GIModel.contact_smearing_sigma(params, m, m)
    expectation = GIModel.radial_expect_udr(
        u,
        r,
        h,
        (ri, i) -> begin
            delta_sigma = sigma^3 / (pi^(3 / 2)) * exp(-(sigma * ri)^2)
            GIModel.alpha_s_r(ri) * delta_sigma
        end,
    )
    manual = (1.0 + params.epsilon_c) * (32 * pi / (9 * m * m)) * expectation * GIModel.spin_dot(3)
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u, r) ≈ manual rtol = 1e-12 atol = 0.0
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
        @test occursin("experimental (A7)", txt)
        @test occursin("3D isotropic smearing", txt)
    end
end

@testset "compare_sector returns shift breakdown fields" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    reference = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    rows = compare_sector(
        params, reference[1:1], "c";
        ngrid = 120, rmax = 12.0, kinetic = :relativistic, contact_hyperfine = true, use_fine_structure = true,
    )
    @test length(rows) == 1
    row = rows[1]
    @test hasproperty(row, :central_GeV)
    @test hasproperty(row, :contact_shift_GeV)
    @test hasproperty(row, :spin_orbit_shift_GeV)
    @test hasproperty(row, :tensor_shift_GeV)
    @test hasproperty(row, :fine_structure_shift_GeV)
    @test row.fine_structure_shift_GeV ≈ row.spin_orbit_shift_GeV + row.tensor_shift_GeV atol = 1e-12
    @test row.predicted_GeV ≈ row.central_GeV + row.contact_shift_GeV + row.fine_structure_shift_GeV atol = 1e-12
end
