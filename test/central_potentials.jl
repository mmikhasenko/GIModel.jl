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
