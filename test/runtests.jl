using Pkg
Pkg.activate(dirname(@__DIR__))

using Test
using FiniteDifferences
using QuadGK
using GIModel

root = dirname(@__DIR__)

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
    @test params.appendix_a_derivative_g == false
    @test params.coulomb_1d_smear == false
    @test params.epsilon_c ≈ -0.168
    @test params.epsilon_t ≈ 0.025
    @test params.epsilon_so_vector ≈ -0.035
    @test params.epsilon_so_scalar ≈ 0.0
    @test params.fine_structure == true
    @test params.k_spin_orbit > 0.0

    ccbar =
        load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    bbbar = load_reference_spectrum(
        joinpath(root, "data", "reference_spectrum_bottomonium.csv"),
    )
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

@testset "Krylov eigensolver matches full eigensolver" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    for kinetic in (:relativistic, :nonrelativistic)
        full, _vec_full, _r_full = channel_solution(
            params,
            m,
            m,
            0;
            nlevels = 3,
            ngrid = 120,
            rmax = 16.0,
            kinetic = kinetic,
            eigensolver = :full,
        )
        krylov, _vec_krylov, _r_krylov = channel_solution(
            params,
            m,
            m,
            0;
            nlevels = 3,
            ngrid = 120,
            rmax = 16.0,
            kinetic = kinetic,
            eigensolver = :krylov,
        )
        @test krylov ≈ full rtol = 1e-10 atol = 1e-10
    end
end

@testset "central: Coulomb+confinement = pointwise V" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    for r0 in (0.15, 0.4, 1.2, 3.0)
        a = GIModel.central_potential(r0, params)
        b = GIModel.static_coulomb_G(r0, params) + GIModel.static_confinement_S(r0, params)
        @test a ≈ b rtol = 1e-12
    end
end

@testset "Coulomb central derivative matches finite difference" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))

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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    p = GIModel.central_potential_path(params)
    @test p.name == "appendix_a_momentum_sandwich"
    @test GIModel.central_potential_mode(params) == :appendix_a_momentum_sandwich
    @test params.appendix_a_smearing == false
    @test params.appendix_a_derivative_g == false
    @test params.appendix_a_closed_form == false
    @test params.appendix_a_momentum_sandwich == true
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
        s2 = replace(
            s,
            "appendix_a_momentum_sandwich = true" => "appendix_a_momentum_sandwich = false",
        )
        s2 = replace(s2, "appendix_a_smearing = false" => "appendix_a_smearing = true")
        write(p, s2)
        params = load_parameters(p)
        @test params.appendix_a_smearing == true
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params,
            mc,
            mc,
            0;
            nlevels = 2,
            ngrid = 120,
            rmax = 12.0,
            kinetic = :relativistic,
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
        s2 = replace(
            s,
            "appendix_a_momentum_sandwich = true" => "appendix_a_momentum_sandwich = false",
        )
        s2 = replace(s2, "coulomb_1d_smear = false" => "coulomb_1d_smear = true")
        write(p, s2)
        params = load_parameters(p)
        @test params.coulomb_1d_smear == true
        @test params.appendix_a_smearing == false
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params,
            mc,
            mc,
            0;
            nlevels = 2,
            ngrid = 120,
            rmax = 12.0,
            kinetic = :relativistic,
        )
        @test isfinite(vals[1]) && isfinite(vals[2])
        @test vals[1] < vals[2]
    end
end

@testset "appendix_a_derivative_g code path (finite S-wave energy)" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        @test occursin("appendix_a_derivative_g = false", s)
        s2 = replace(
            s,
            "appendix_a_momentum_sandwich = true" => "appendix_a_momentum_sandwich = false",
        )
        s2 = replace(
            s2,
            "appendix_a_derivative_g = false" => "appendix_a_derivative_g = true",
        )
        write(p, s2)
        params = load_parameters(p)
        @test params.appendix_a_derivative_g == true
        @test GIModel.central_potential_mode(params) == :appendix_a_derivative_g
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params,
            mc,
            mc,
            0;
            nlevels = 2,
            ngrid = 120,
            rmax = 12.0,
            kinetic = :relativistic,
        )
        @test isfinite(vals[1]) && isfinite(vals[2])
        @test vals[1] < vals[2]
    end
end

@testset "appendix_a_closed_form code path (finite S-wave energy)" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        @test occursin("appendix_a_closed_form = false", s)
        s2 = replace(
            s,
            "appendix_a_momentum_sandwich = true" => "appendix_a_momentum_sandwich = false",
        )
        s2 =
            replace(s2, "appendix_a_closed_form = false" => "appendix_a_closed_form = true")
        write(p, s2)
        params = load_parameters(p)
        @test params.appendix_a_closed_form == true
        @test GIModel.central_potential_mode(params) == :appendix_a_closed_form
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params,
            mc,
            mc,
            0;
            nlevels = 2,
            ngrid = 120,
            rmax = 12.0,
            kinetic = :relativistic,
        )
        @test isfinite(vals[1]) && isfinite(vals[2])
        @test vals[1] < vals[2]
    end
end

@testset "appendix_a_momentum_sandwich code path (finite S-wave energy)" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        @test occursin("appendix_a_momentum_sandwich = true", s)
        write(p, s)
        params = load_parameters(p)
        @test params.appendix_a_momentum_sandwich == true
        @test GIModel.central_potential_mode(params) == :appendix_a_momentum_sandwich
        mc = params.masses["c"]
        vals, _v, _r = GIModel.channel_solution(
            params,
            mc,
            mc,
            0;
            nlevels = 2,
            ngrid = 80,
            rmax = 10.0,
            kinetic = :relativistic,
        )
        @test isfinite(vals[1]) && isfinite(vals[2])
        @test vals[1] < vals[2]
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

@testset "central_potential_values named modes" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    mc = params.masses["c"]
    r, _h = GIModel.radial_grid(80, 8.0)
    pointwise = GIModel.central_potential_values(params, mc, mc, r; mode = :pointwise)
    @test pointwise ≈ [GIModel.central_potential(ri, params) for ri in r]
    @test GIModel.central_potential_values(params, mc, mc, r; mode = :coulomb_1d) ≈
          GIModel.coulomb_1d_smeared_central_values(params, mc, mc, r)
    derivative =
        GIModel.central_potential_values(params, mc, mc, r; mode = :appendix_a_derivative_g)
    @test all(isfinite, derivative)
    @test maximum(abs.(derivative .- pointwise)) > 0.0
    closed =
        GIModel.central_potential_values(params, mc, mc, r; mode = :appendix_a_closed_form)
    @test all(isfinite, closed)
    @test maximum(abs.(closed .- pointwise)) > 0.0
    sandwich_view = GIModel.central_potential_values(
        params,
        mc,
        mc,
        r;
        mode = :appendix_a_momentum_sandwich,
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
    @test_throws ErrorException GIModel.central_potential_values(
        params,
        mc,
        mc,
        r;
        mode = :unknown_mode,
    )
end

@testset "appendix_a_momentum_sandwich wins over diagonal central modes" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s2 = s
        s2 =
            replace(s2, "appendix_a_closed_form = false" => "appendix_a_closed_form = true")
        s2 = replace(
            s2,
            "appendix_a_derivative_g = false" => "appendix_a_derivative_g = true",
        )
        s2 = replace(s2, "appendix_a_smearing = false" => "appendix_a_smearing = true")
        s2 = replace(s2, "coulomb_1d_smear = false" => "coulomb_1d_smear = true")
        write(p, s2)
        params = load_parameters(p)
        path = GIModel.central_potential_path(params)
        @test path.name == "appendix_a_momentum_sandwich"
    end
end

@testset "appendix_a_smearing wins over coulomb_1d in central_potential_path" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s2 = replace(
            s,
            "appendix_a_momentum_sandwich = true" => "appendix_a_momentum_sandwich = false",
        )
        s2 = replace(s2, "appendix_a_smearing = false" => "appendix_a_smearing = true")
        s2 = replace(s2, "coulomb_1d_smear = false" => "coulomb_1d_smear = true")
        write(p, s2)
        params = load_parameters(p)
        @test params.appendix_a_smearing == true
        @test params.coulomb_1d_smear == true
        path = GIModel.central_potential_path(params)
        @test path.name == "experimental_3d_convl_a7a8"
    end
end

@testset "appendix_a_derivative_g wins over older central smearing toggles" begin
    mktempdir() do d
        p = joinpath(d, "p.toml")
        s = read(joinpath(root, "data", "parameters.provisional.toml"), String)
        s2 = replace(
            s,
            "appendix_a_momentum_sandwich = true" => "appendix_a_momentum_sandwich = false",
        )
        s2 = replace(
            s2,
            "appendix_a_derivative_g = false" => "appendix_a_derivative_g = true",
        )
        s2 = replace(s2, "appendix_a_smearing = false" => "appendix_a_smearing = true")
        s2 = replace(s2, "coulomb_1d_smear = false" => "coulomb_1d_smear = true")
        write(p, s2)
        params = load_parameters(p)
        @test params.appendix_a_derivative_g == true
        @test params.appendix_a_smearing == true
        @test params.coulomb_1d_smear == true
        path = GIModel.central_potential_path(params)
        @test path.name == "appendix_a_derivative_g"
    end
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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    σ = GIModel.contact_smearing_sigma(params, m, m)

    function quadgk_smeared_coulomb(R)
        total = 0.0
        for (α, γ) in zip(GIModel.ALPHA_COEFFS, GIModel.ALPHA_GAMMAS)
            integrand(rp) = begin
                g =
                    rp == 0 ? -8 * α * γ / (3 * sqrt(π)) :
                    -4 * α * GIModel.gi_erf(γ * rp) / (3 * rp)
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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _, umat, r =
        GIModel.channel_solution(params, m, m, 0; nlevels = 2, ngrid = 200, rmax = 20.0)
    h = r[2] - r[1]
    u_s = collect(umat[:, 1])
    @test GIModel.fine_structure_split(
        params,
        m,
        m,
        "S",
        3,
        1,
        u_s,
        r,
        h;
        k_spin_orbit = params.k_spin_orbit,
        k_tensor = params.k_tensor,
    ) == 0.0
    @test GIModel.fine_structure_split(
        params,
        m,
        m,
        "S",
        1,
        0,
        u_s,
        r,
        h;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    ) == 0.0
    v_p, umat_p, r_p =
        GIModel.channel_solution(params, m, m, 1; nlevels = 2, ngrid = 200, rmax = 20.0)
    h_p = r_p[2] - r_p[1]
    u1p = collect(umat_p[:, 1])
    δ0 = GIModel.fine_structure_split(
        params,
        m,
        m,
        "P",
        3,
        0,
        u1p,
        r_p,
        h_p;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    )
    δ1 = GIModel.fine_structure_split(
        params,
        m,
        m,
        "P",
        3,
        1,
        u1p,
        r_p,
        h_p;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    )
    δ2 = GIModel.fine_structure_split(
        params,
        m,
        m,
        "P",
        3,
        2,
        u1p,
        r_p,
        h_p;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    )
    @test isfinite(δ0) && isfinite(δ1) && isfinite(δ2)
    @test δ0 != δ1 || δ1 != δ2
    @test GIModel.fine_structure_split(
        params,
        m,
        m,
        "P",
        1,
        0,
        u1p,
        r_p,
        h_p;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    ) == 0.0
end

@testset "fine_structure_components: decomposition sums correctly" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _v_p, umat_p, r_p =
        GIModel.channel_solution(params, m, m, 1; nlevels = 2, ngrid = 200, rmax = 20.0)
    h_p = r_p[2] - r_p[1]
    u1p = collect(umat_p[:, 1])
    for J in (0, 1, 2)
        comp = GIModel.fine_structure_components(
            params,
            m,
            m,
            "P",
            3,
            J,
            u1p,
            r_p,
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
        )
        @test comp.spin_orbit ≈ comp.spin_orbit_vector + comp.spin_orbit_thomas atol = 1e-12
        @test comp.total ≈ comp.spin_orbit + comp.tensor atol = 1e-12
        @test comp.total ≈ GIModel.fine_structure_split(
            params,
            m,
            m,
            "P",
            3,
            J,
            u1p,
            r_p,
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
        ) atol = 1e-12
    end
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

        params0 = load_parameters(p0)
        mc = params0.masses["c"]
        _v_p, umat_p, r_p = GIModel.channel_solution(
            params0,
            mc,
            mc,
            1;
            nlevels = 2,
            ngrid = 200,
            rmax = 20.0,
        )
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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

        params0 = load_parameters(p0)
        mc = params0.masses["c"]
        _v_p, umat_p, r_p = GIModel.channel_solution(
            params0,
            mc,
            mc,
            1;
            nlevels = 2,
            ngrid = 200,
            rmax = 20.0,
        )
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
            h_p;
            k_spin_orbit = 0.0,
            k_tensor = 1.0,
        )
        @test comp.I_tk ≈ Itk rtol = 1e-12 atol = 0.0
        expected = (1.0 / (3.0 * mc * mc)) * Itk * GIModel.tensor_triplet_LJ(1, 2, 1)
        @test comp.tensor ≈ expected rtol = 1e-12 atol = 0.0
    end
end

@testset "contact_hyperfine_shift: normalization + spin algebra" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _vals, umat, r =
        GIModel.channel_solution(params, m, m, 0; nlevels = 2, ngrid = 200, rmax = 20.0)
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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    _v_p, umat_p, r_p =
        GIModel.channel_solution(params, m, m, 1; nlevels = 2, ngrid = 200, rmax = 20.0)
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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
            h_p;
            k_spin_orbit = 1.0,
            k_tensor = 1.0,
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
    @test GIModel.gi_erf(0.0) ≈ 0.0 atol = 1e-7
    for x in (0.05, 0.3, 0.8, 1.5)
        @test GIModel.gi_erf(x) + GIModel.gi_erf(-x) ≈ 0.0 atol = 1e-12
    end
end

@testset "GI erf profile second derivative consistency" begin
    for x in (0.1, 0.4, 1.1, 1.8)
        δ = 1e-6 * max(1.0, x)
        num = (GIModel.gi_erf_prime(x + δ) - GIModel.gi_erf_prime(x - δ)) / (2δ)
        ana = GIModel.gi_erf_second(x)
        @test ana ≈ num rtol = 2e-5 atol = 1e-9
        @test GIModel.gi_erf_second(-x) ≈ -ana rtol = 1e-12 atol = 1e-12
    end
end

@testset "Coulomb derivative consistency (GI erf profile)" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
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
        params = load_parameters(p)
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
        m = params.masses["c"]
        comp = GIModel.fine_structure_components(
            params,
            m,
            m,
            "P",
            3,
            2,
            collect(u),
            collect(r),
            h;
            k_spin_orbit = 1.0,
            k_tensor = 0.0,
        )
        @test comp.I_cm ≈ I_cm rtol = 1e-12 atol = 0.0
        inv2_cm = 0.5 * (1.0 / m^2 + 1.0 / m^2 + 2.0 / (m * m))
        ls = GIModel.LdotS(1, 1, 2)
        expected_vec = inv2_cm * ls * (1.0 + params.epsilon_so_vector) * I_cm
        @test comp.spin_orbit_vector ≈ expected_vec rtol = 1e-12 atol = 0.0

        I_tp = GIModel.radial_expect_udr(
            u,
            r,
            h,
            (ri, i) -> begin
                r0 = max(ri, 1.0e-8)
                (1.0 / (2.0 * r0)) * (params.b + GIModel.dV_coul_central_dr(r0, params))
            end,
        )
        @test comp.I_tp ≈ I_tp rtol = 1e-12 atol = 0.0
        inv2_tp = 0.5 * (1.0 / m^2 + 1.0 / m^2)
        expected_tp = (-inv2_tp) * ls * (1.0 + params.epsilon_so_scalar) * I_tp
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
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
    r = collect(0.05:0.05:1.0)
    u = exp.(-2.0 .* r)
    @test GIModel.contact_hyperfine_shift(params, m, m, "P", 3, u, r) == 0.0
    @test GIModel.contact_hyperfine_shift(params, m, m, "D", 3, u, r) == 0.0
end

@testset "contact smearing σ implements Appendix A (A9)" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m1, m2 = params.masses["c"], params.masses["b"]

    σ = GIModel.contact_smearing_sigma(params, m1, m2)
    σ_swapped = GIModel.contact_smearing_sigma(params, m2, m1)
    @test σ ≈ σ_swapped rtol = 0.0 atol = 0.0

    mass_factor = 4 * m1 * m2 / (m1 + m2)^2
    reduced_twice = 2 * m1 * m2 / (m1 + m2)
    σ_manual = sqrt(
        params.sigma0^2 * (0.5 + 0.5 * mass_factor^4) +
        params.smearing_s^2 * reduced_twice^2,
    )
    @test σ ≈ σ_manual rtol = 0.0 atol = 0.0

    m = params.masses["c"]
    σ_equal = GIModel.contact_smearing_sigma(params, m, m)
    σ_equal_manual = sqrt(params.sigma0^2 + params.smearing_s^2 * m^2)
    @test σ_equal ≈ σ_equal_manual rtol = 0.0 atol = 0.0
end

@testset "contact hyperfine shift uses u(r) normalization" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m = params.masses["c"]
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
    manual =
        (1.0 + params.epsilon_c) *
        (32 * pi / (9 * m * m)) *
        expectation *
        GIModel.spin_dot(3)
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u, r) ≈ manual rtol = 1e-12 atol =
        0.0
end

@testset "RadialChannelKey collapses nearly-equal masses" begin
    a = RadialChannelKey(1.628, 1.628, "S")
    b = RadialChannelKey(1.628 + 1e-20, 1.628 + 1e-20, "S")
    @test a == b
    @test hash(a) == hash(b)
    @test RadialChannelKey(1.6, 1.6, "S") != RadialChannelKey(1.6, 1.6, "P")
end

@testset "compute_sector + compare return shift breakdown fields" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    reference =
        load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    sub = reference[1:1]
    computed = compute_sector(
        params,
        sub,
        "c";
        ngrid = 120,
        rmax = 12.0,
        kinetic = :relativistic,
    )
    rows = compare(
        computed,
        sub;
        contact_hyperfine = true,
        use_fine_structure = true,
    )
    @test length(rows) == 1
    row = rows[1]
    @test hasproperty(row, :central_GeV)
    @test hasproperty(row, :contact_shift_GeV)
    @test hasproperty(row, :spin_orbit_shift_GeV)
    @test hasproperty(row, :tensor_shift_GeV)
    @test hasproperty(row, :fine_structure_shift_GeV)
    @test hasproperty(row, :m1_GeV)
    @test hasproperty(row, :m2_GeV)
    @test hasproperty(row, :fine_structure_mass_convention)
    @test isfinite(row.m1_GeV) && isfinite(row.m2_GeV)
    @test row.fine_structure_mass_convention in
          ("equal_mass", "unequal_mass_equal_share_LdotS", "disabled")
    @test row.fine_structure_shift_GeV ≈ row.spin_orbit_shift_GeV + row.tensor_shift_GeV atol =
        1e-12
    @test row.predicted_GeV ≈
          row.central_GeV + row.contact_shift_GeV + row.fine_structure_shift_GeV atol =
        1e-12
end

@testset "compute_sector: empty reference" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    empty_ref = ReferenceState[]
    computed = compute_sector(
        params,
        empty_ref,
        "c";
        ngrid = 80,
        rmax = 8.0,
        kinetic = :relativistic,
    )
    @test isempty(computed.channel_cache)
    @test isempty(compare(computed, empty_ref))
end

@testset "compare reuses SectorComputation (stable + superset cache)" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    reference =
        load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    prefix = reference[1:min(4, length(reference))]
    isempty(prefix) && error("charmonium reference unexpectedly empty")
    computed = compute_sector(
        params,
        prefix,
        "c";
        ngrid = 100,
        rmax = 10.0,
        kinetic = :relativistic,
    )
    r1 = compare(computed, prefix; contact_hyperfine = true, use_fine_structure = true)
    r2 = compare(computed, prefix; contact_hyperfine = true, use_fine_structure = true)
    @test length(r1) == length(prefix)
    @test length(r2) == length(prefix)
    for i = 1:length(r1)
        @test r1[i].predicted_GeV ≈ r2[i].predicted_GeV rtol = 0.0 atol = 1e-15
        @test r1[i].residual_MeV ≈ r2[i].residual_MeV rtol = 0.0 atol = 1e-12
    end

    short = prefix[1:min(2, length(prefix))]
    from_superset = compare(computed, short; contact_hyperfine = true, use_fine_structure = true)
    computed_short = compute_sector(
        params,
        short,
        "c";
        ngrid = 100,
        rmax = 10.0,
        kinetic = :relativistic,
    )
    direct = compare(computed_short, short; contact_hyperfine = true, use_fine_structure = true)
    @test length(from_superset) == length(direct)
    for i = 1:length(direct)
        @test from_superset[i].predicted_GeV ≈ direct[i].predicted_GeV rtol = 1e-12 atol = 0.0
        @test from_superset[i].central_GeV ≈ direct[i].central_GeV rtol = 1e-12 atol = 0.0
    end
end

@testset "compare(contact_hyperfine=false) drops contact shift for covered states" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    reference =
        load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    sub = reference[1:1]
    computed = compute_sector(
        params,
        sub,
        "c";
        ngrid = 120,
        rmax = 12.0,
        kinetic = :relativistic,
    )
    with_hf =
        compare(computed, sub; contact_hyperfine = true, use_fine_structure = false)
    no_hf =
        compare(computed, sub; contact_hyperfine = false, use_fine_structure = false)
    @test length(with_hf) == 1 && length(no_hf) == 1
    @test no_hf[1].contact_shift_GeV ≈ 0.0 atol = 1e-15
    @test with_hf[1].central_GeV ≈ no_hf[1].central_GeV rtol = 1e-12 atol = 0.0
    @test with_hf[1].predicted_GeV - no_hf[1].predicted_GeV ≈ with_hf[1].contact_shift_GeV atol =
        1e-12
end
