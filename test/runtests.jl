# Invoked by `Pkg.test()` with the package environment already active.
# To run manually: `julia --project=. test/runtests.jl` from the repo root.

using Test
using FiniteDifferences
using LinearAlgebra
using QuadGK
using SpecialFunctions: erf
using GIModel

root = dirname(@__DIR__)

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
    @test params.fine_structure.k_spin_orbit > 0.0
end

@testset "baseline solver shape" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    cc = solve_sector(params, mq["c"]; maxn = 4, ngrid = 250, rmax = 20.0)
    bb = solve_sector(params, mq["b"]; maxn = 4, ngrid = 250, rmax = 16.0)

    @test cc[(1, "S")] < cc[(2, "S")] < cc[(3, "S")]
    @test bb[(1, "S")] < bb[(2, "S")] < bb[(3, "S")]
    @test cc[(1, "S")] < cc[(1, "P")] < cc[(1, "D")]
    @test bb[(1, "S")] < bb[(1, "P")] < bb[(1, "D")]
end

@testset "Krylov eigensolver matches full eigensolver" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    for kinetic in (:relativistic, :nonrelativistic)
        full, _vec_full, _r_full = channel_solution(
            params,
            ConstituentMasses(m, m),
            0;
            nlevels = 3,
            ngrid = 120,
            rmax = 16.0,
            kinetic = kinetic,
            eigensolver = :full,
        )
        krylov, _vec_krylov, _r_krylov = channel_solution(
            params,
            ConstituentMasses(m, m),
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
    @test params isa GIParameters{FiniteDifferenceBasis}
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

@testset "GIParameters basis dispatch keeps FD p² path explicit" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    r, h = GIModel.radial_grid(12, 3.0)
    p2_default = GIModel.p2_operator(mq["c"], 1, r, h)
    p2_param = GIModel.p2_operator(params, mq["c"], 1, r, h)
    p2_basis = GIModel.p2_operator(FiniteDifferenceBasis, mq["c"], 1, r, h)
    @test p2_param ≈ p2_default
    @test p2_basis ≈ p2_default
end

@testset "HarmonicOscillatorBasis channel solve returns mesh wavefunctions" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    params_ho = GIModel.with_basis(params, HarmonicOscillatorBasis)
    vals, vecs, r = channel_solution(
        params_ho,
        ConstituentMasses(mq["c"], mq["c"]),
        0;
        nlevels = 3,
        ngrid = 120,
        rmax = 14.0,
    )
    @test params_ho isa GIParameters{HarmonicOscillatorBasis}
    @test length(vals) == 3
    @test size(vecs) == (length(r), 3)
    @test vals[1] < vals[2] < vals[3]
    h = r[2] - r[1]
    for col in axes(vecs, 2)
        @test sum(abs2, vecs[:, col]) * h ≈ 1.0 rtol = 1e-10
    end
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
        vals, _v, _r = GIModel.channel_solution(
            params,
            ConstituentMasses(mc, mc),
            0;
            nlevels = 2,
            ngrid = ngrid,
            rmax = rmax,
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
    _, umat, r =
        GIModel.channel_solution(params, ConstituentMasses(m, m), 0; nlevels = 2, ngrid = 200, rmax = 20.0)
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
        k_spin_orbit = params.fine_structure.k_spin_orbit,
        k_tensor = params.fine_structure.k_tensor,
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
        GIModel.channel_solution(params, ConstituentMasses(m, m), 1; nlevels = 2, ngrid = 200, rmax = 20.0)
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
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    _v_p, umat_p, r_p =
        GIModel.channel_solution(params, ConstituentMasses(m, m), 1; nlevels = 2, ngrid = 200, rmax = 20.0)
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
    @test_throws ArgumentError GIModel.MixingBlock(
        "bad asymmetric block",
        [GIModel.BasisState(1, "S", 3, 1), GIModel.BasisState(1, "D", 3, 1)],
        [1.0 0.2; 0.1 2.0],
    )

    _vals_cc, umat_cc, r_cc = GIModel.channel_solution(
        params,
        ConstituentMasses(mc, mc),
        1;
        nlevels = 1,
        ngrid = 200,
        rmax = 20.0,
    )
    radial_cc = RadialWaveOnUniformMesh(umat_cc[:, 1], r_cc)
    off_cc = GIModel.spin_orbit_mixing_components(
        params,
        ConstituentMasses(mc, mc),
        "P",
        radial_cc;
        k_spin_orbit = params.fine_structure.k_spin_orbit,
    )
    @test off_cc.total == 0.0
    mix_cc = GIModel.same_j_mixing(3.5, 3.6, off_cc.total)
    @test mix_cc.masses ≈ [3.5, 3.6]
    @test mix_cc.theta_deg ≈ 0.0 atol = 1e-12
    @test mix_cc.block isa GIModel.MixingBlock
    @test mix_cc.block.mechanism == "antisymmetric_spin_orbit"

    vals_bc, umat_bc, r_bc = GIModel.channel_solution(
        params,
        ConstituentMasses(mb, mc),
        1;
        nlevels = 1,
        ngrid = 200,
        rmax = 24.0,
    )
    radial_bc = RadialWaveOnUniformMesh(umat_bc[:, 1], r_bc)
    off_bc = GIModel.spin_orbit_mixing_components(
        params,
        ConstituentMasses(mb, mc),
        "P",
        radial_bc;
        k_spin_orbit = params.fine_structure.k_spin_orbit,
    )
    @test isfinite(off_bc.total)
    @test off_bc.total != 0.0

    triplet_shift = GIModel.fine_structure_split(
        params,
        ConstituentMasses(mb, mc),
        FineStructureMultiplet("P", 3, 1),
        radial_bc;
        k_spin_orbit = params.fine_structure.k_spin_orbit,
        k_tensor = params.fine_structure.k_tensor,
    )
    mix_bc = GIModel.same_j_mixing(vals_bc[1], vals_bc[1] + triplet_shift, off_bc.total)
    @test isfinite(mix_bc.theta_deg)
    @test minimum(mix_bc.masses) < vals_bc[1] < maximum(mix_bc.masses)
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
        _v_p, umat_p, r_p = GIModel.channel_solution(
            params0,
            ConstituentMasses(mc, mc),
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

        params0, mq0 = load_parameters_and_quark_masses(p0)
        mc = mq0["c"]
        _v_p, umat_p, r_p = GIModel.channel_solution(
            params0,
            ConstituentMasses(mc, mc),
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

@testset "tensor off-diagonal same-J angular factor" begin
    @test GIModel.tensor_triplet_offdiag_sameJ(0, 1) ≈ 0.0
    @test GIModel.tensor_triplet_offdiag_sameJ(1, 1) ≈ sqrt(8.0)
    @test GIModel.tensor_triplet_offdiag_sameJ(2, 1) ≈ 6sqrt(6.0) / 5
    @test GIModel.tensor_triplet_offdiag_sameJ(1, 0) ≈ 0.0
end

@testset "contact_hyperfine_shift: normalization + spin algebra" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    _vals, umat, r =
        GIModel.channel_solution(params, ConstituentMasses(m, m), 0; nlevels = 2, ngrid = 200, rmax = 20.0)
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
    _v_p, umat_p, r_p =
        GIModel.channel_solution(params, ConstituentMasses(m, m), 1; nlevels = 2, ngrid = 200, rmax = 20.0)
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
            h;
            k_spin_orbit = 1.0,
            k_tensor = 0.0,
        )
        @test comp.I_cm ≈ I_cm rtol = 1e-12 atol = 0.0
        inv2_cm = 0.5 * (1.0 / m^2 + 1.0 / m^2 + 2.0 / (m * m))
        ls = GIModel.LdotS(1, 1, 2)
        expected_vec = inv2_cm * ls * (1.0 + params.factors.epsilon_so_vector) * I_cm
        @test comp.spin_orbit_vector ≈ expected_vec rtol = 1e-12 atol = 0.0

        I_tp = GIModel.radial_expect_udr(
            u,
            r,
            h,
            (ri, i) -> begin
                r0 = max(ri, 1.0e-8)
                (1.0 / (2.0 * r0)) * (params.potential.b + GIModel.dV_coul_central_dr(r0, params))
            end,
        )
        @test comp.I_tp ≈ I_tp rtol = 1e-12 atol = 0.0
        inv2_tp = 0.5 * (1.0 / m^2 + 1.0 / m^2)
        expected_tp = (-inv2_tp) * ls * (1.0 + params.factors.epsilon_so_scalar) * I_tp
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
        (1.0 + params.factors.epsilon_c) *
        (32 * pi / (9 * m * m)) *
        expectation *
        GIModel.spin_dot(3)
    @test GIModel.contact_hyperfine_shift(params, m, m, "S", 3, u, r) ≈ manual rtol = 1e-12 atol =
        0.0
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
    b = RadialChannelKey(1.628 + 1e-20, 1.628 + 1e-20, "S")
    @test a == b
    @test hash(a) == hash(b)
    @test RadialChannelKey(1.6, 1.6, "S") != RadialChannelKey(1.6, 1.6, "P")
    masses = GIModel.ConstituentMasses(1.6, 1.6)
    @test RadialChannelKey(masses, "S") == RadialChannelKey(1.6, 1.6, "S")
end

@testset "RadialWaveOnUniformMesh agrees with ChannelRadialSolution column" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    m = mq["c"]
    ev, vecs, r =
        GIModel.channel_solution(params, ConstituentMasses(m, m), 1; nlevels = 2, ngrid = 120, rmax = 16.0)
    sol = GIModel.ChannelRadialSolution(ev, vecs, r)
    wave = GIModel.RadialWaveOnUniformMesh(sol, 1)
    h = r[2] - r[1]
    @test wave.u ≈ vecs[:, 1]
    @test wave.r ≈ r
    @test wave.h ≈ h
    mult = GIModel.FineStructureMultiplet("P", 3, 2)
    masses = GIModel.ConstituentMasses(m, m)
    c1 = GIModel.fine_structure_components(
        params,
        masses,
        mult,
        wave;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    )
    c2 = GIModel.fine_structure_components(
        params,
        m,
        m,
        "P",
        3,
        2,
        collect(vecs[:, 1]),
        r,
        h;
        k_spin_orbit = 1.0,
        k_tensor = 1.0,
    )
    @test c1.total ≈ c2.total rtol = 1e-12 atol = 0.0
end

@testset "mixing mechanisms are comparison-layer markers" begin
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

@testset "Every solve returns u with the same normalization" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    ho = with_basis(params, HarmonicOscillatorBasis)

    # The central solve, in both bases and every orbital: integral u^2 dr = 1.
    # FD used to return Euclidean eigenvectors here (sum u^2 = 1) while HO
    # returned physical ones, differing by sqrt(h).
    for key in ("q", "s", "c", "b"), L in 0:2
        masses = ConstituentMasses(mq[key], mq[key])
        for prm in (params, ho)
            _e, u, r = channel_solution(prm, masses, L; nlevels = 2)
            h = r[2] - r[1]
            for col in axes(u, 2)
                @test isapprox(sum(abs2, view(u, :, col)) * h, 1.0; rtol = 1e-10)
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

@testset "Both bases solve H+V to the same quantity, same normalization" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    ho = with_basis(params, HarmonicOscillatorBasis)
    r, h = GIModel.radial_grid(450, 24.0)

    for key in ("q", "c", "b")
        masses = ConstituentMasses(mq[key], mq[key])
        V = Matrix(GIModel.contact_hyperfine_operator(params, masses, "S", 1, r))
        e_fd, u_fd, r_fd = resummed_channel_solution(params, masses, 0, V; nlevels = 2)
        e_ho, u_ho, r_ho = resummed_channel_solution(ho, masses, 0, V; nlevels = 2)

        # Same mesh out.
        @test r_fd == r_ho == collect(r)
        # Same normalization: physical, integral u^2 dr = 1, for BOTH bases.
        # FD used to return Euclidean eigenvectors (sum u^2 = 1), differing from
        # HO by exactly sqrt(h); anything quadratic in u was then off by h.
        for u in (u_fd, u_ho), col in axes(u, 2)
            @test isapprox(sum(abs2, u[:, col]) * h, 1.0; rtol = 1e-10)
        end
        # Same quantity: the two algorithms agree on the ground-state energy to
        # well under an MeV across light, charm and bottom.
        @test abs(e_ho[1] - e_fd[1]) < 1e-3
    end

    # The original exported name is the oscillator method of the unified solve.
    masses = ConstituentMasses(mq["c"], mq["c"])
    V = Matrix(GIModel.contact_hyperfine_operator(params, masses, "S", 1, r))
    @test ho_full_distorted_states(ho, masses, 0, V; nlevels = 2)[1] ==
          resummed_channel_solution(ho, masses, 0, V; nlevels = 2)[1]

    # V must live on the solver's mesh, in either basis.
    bad = zeros(10, 10)
    @test_throws Exception resummed_channel_solution(params, masses, 0, bad; nlevels = 1)
    @test_throws Exception resummed_channel_solution(ho, masses, 0, bad; nlevels = 1)
end

@testset "An unimplemented solve fails instead of degrading" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["q"], mq["q"])
    r, _ = GIModel.radial_grid(450, 24.0)

    # FD has the non-perturbative contact solve, and returns it.
    levels, vectors, _ = contact_hyperfine_nonperturbative_states(params, masses, "S", 1, r, 2)
    @test length(levels) == 2 && size(vectors, 2) == 2

    # Empty is still the right answer where the path is genuinely inactive:
    # not an S wave, or a multiplicity the contact term does not touch.
    @test contact_hyperfine_nonperturbative_states(params, masses, "P", 1, r, 2)[1] == Float64[]
    @test GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 2, r, 2) == Float64[]

    # But a basis with no implementation must NOT answer "empty". It used to,
    # and `add_spin_corrections` read that as "fall back to first-order PT",
    # which for the light 1S0 gave 0.2842 GeV against the resummed 0.0950 GeV
    # — a 3x wrong pion with no warning.
    ho = with_basis(params, HarmonicOscillatorBasis)
    @test_throws ArgumentError contact_hyperfine_nonperturbative_states(ho, masses, "S", 1, r, 2)
    @test_throws ArgumentError GIModel.contact_hyperfine_nonperturbative_levels(ho, masses, "S", 1, r, 2)
    @test_throws ArgumentError compute_spectrum(ho, Meson(mq, :q, :q); levels = spectrum_levels(1))

    # The FD front door is untouched and still resums.
    fd_pi = spectrum_state(
        compute_spectrum(params, Meson(mq, :q, :q); levels = spectrum_levels(1)), "1^1S_0")
    @test fd_pi.mass_GeV < 0.15        # resummed; first-order PT lands near 0.28
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
        ngrid = 450, rmax = 24.0, kinetic = :relativistic, eigensolver = :full,
        nlevels_per_channel = 6, contact_hyperfine = true, use_fine_structure = true,
        same_j_spin_orbit_mixing = true, tensor_mixing = true)) == masses(base)
    @test masses(compute_spectrum(params, meson; levels = levels,
        solver = RadialSolver(), terms = SpinTerms())) == masses(base)

    # RadialSolver is "how well": refining the mesh must not move a mass more
    # than the discretization error it removes (sub-MeV here).
    fine = compute_spectrum(params, meson; levels = levels,
        solver = RadialSolver(ngrid = 900, rmax = 32.0))
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
    ev, _, _ = channel_solution(params, meson.constituent_masses, 0;
        solver = RadialSolver(nlevels_per_channel = 3))
    @test length(ev) == 3
    tuned = RadialSolver(RadialSolver(ngrid = 900); rmax = 32.0)
    @test tuned.ngrid == 900 && tuned.rmax == 32.0 && tuned.kinetic === :relativistic

    # The deprecated loose keywords still work, still win, and now say so.
    # Clean paths must stay silent — the shims fold legacy keywords in once at
    # the public entry point, so internal forwarding cannot self-trigger.
    @test_logs compute_spectrum(params, meson; levels = levels)
    @test_logs compute_spectrum(params, meson; levels = levels,
        solver = RadialSolver(), terms = SpinTerms())
    @test_logs central_spectrum(params, meson; levels = levels)
    @test_logs channel_solution(params, meson.constituent_masses, 0;
        solver = RadialSolver())
    @test (@test_logs (:warn,) match_mode = :any masses(
        compute_spectrum(params, meson; levels = levels, ngrid = 450))) == masses(base)
    @test (@test_logs (:warn,) match_mode = :any masses(
        compute_spectrum(params, meson; levels = levels, tensor_mixing = true))) == masses(base)
    # The silent-override edge the warning exists to expose: the loose keyword
    # beats the object, so this is the 450 answer, not the 900 one.
    @test (@test_logs (:warn,) match_mode = :any masses(compute_spectrum(params, meson;
        levels = levels, solver = RadialSolver(ngrid = 900), ngrid = 450))) == masses(base)

    # Invalid settings are construction errors, not silent fallbacks.
    @test_throws ArgumentError RadialSolver(kinetic = :newtonian)
    @test_throws ArgumentError RadialSolver(eigensolver = :lanczos)
    @test_throws ArgumentError RadialSolver(ngrid = 1)
    @test_throws ArgumentError RadialSolver(rmax = 0)
    @test_throws ArgumentError RadialSolver(nlevels_per_channel = 0)
end

@testset "Resolution walls are detected, not silent" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))

    # Grid points spanned by the ground state's RMS radius. The FD grid is fixed,
    # so a compact enough state falls between points and comes back silently
    # under-resolved (at 30 GeV the hyperfine splitting collapses to 0.0025).
    function points_across(m; ngrid = 450, rmax = 24.0)
        vals, vecs, r = channel_solution(params, ConstituentMasses(m, m), 0;
            nlevels = 2, ngrid = ngrid, rmax = rmax)
        h = r[2] - r[1]
        u = vecs[:, 1]
        rms = sqrt(sum(abs2.(u) .* r .^ 2) * h / (sum(abs2, u) * h))
        return rms / h
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

    # Same story in the oscillator basis: the fixed HO_BETA_GRID rails, i.e. the
    # variational optimum lands on the last candidate, only well above bottomonium.
    function best_beta(m)
        ph = with_basis(params, HarmonicOscillatorBasis)
        mm = ConstituentMasses(m, m)
        r, h = GIModel.radial_grid(450, 24.0)
        best = nothing
        for b in GIModel.HO_BETA_GRID
            H, _ = GIModel.oscillator_hamiltonian_for_beta(ph, mm, 0, r, h, b; nbasis = 24)
            v, _ = GIModel.lowest_eigenpairs(Matrix(H), 2; eigensolver = :full)
            (isnothing(best) || v[end] < best[2]) && (best = (b, v[end]))
        end
        return best[1]
    end
    @test best_beta(mq["b"]) < last(GIModel.HO_BETA_GRID)     # bottomonium is fine
    @test best_beta(30.0) == last(GIModel.HO_BETA_GRID)       # railed

    # And the detector is wired to warn (maxlog=1, so this is the first trigger).
    @test_logs (:warn,) match_mode = :any channel_solution(
        params, ConstituentMasses(30.0, 30.0), 0; nlevels = 2)
end

@testset "Isoscalar coherence factor is stated, not string-matched" begin
    r = collect(range(0.05, 6.0; length = 64))
    wave = RadialWaveOnUniformMesh(exp.(-r), r)
    mk(label, coh) = pseudoscalar_annihilation_basis_input(label, 0.22, 0.9, wave;
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

    # The deprecated label inference still reproduces the historical rule for
    # unmigrated callers (archived forensics), and warns.
    @test GIModel._label_implies_coherent("1 ns")
    @test GIModel._label_implies_coherent("2 n nbar")
    @test !GIModel._label_implies_coherent("1 ss")
    @test (@test_logs (:warn,) match_mode = :any pseudoscalar_annihilation_basis_input(
        "1 ns", 0.22, 0.9, wave)).isoscalar_coherent
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
    spec = compute_spectrum(params, meson; levels = levels, ngrid = 120, rmax = 12.0)
    @test length(spec.states) == 4
    for s in spec.states
        if isempty(s.mixings)
            @test s.mass_GeV ≈ s.central_GeV + s.contact_shift_GeV + s.fine_structure_shift_GeV atol = 1e-12
        end
        @test s.fine_structure_shift_GeV ≈ s.spin_orbit_shift_GeV + s.tensor_shift_GeV atol = 1e-12
    end
    # 2^3S_1 and 1^3D_1 form a tensor block; trace is conserved
    mixed = [s for s in spec.states if !isempty(s.mixings)]
    @test length(mixed) == 2
    @test all(m.mechanism == "tensor_mixing" for s in mixed for m in s.mixings)
    @test sum(s.mass_GeV for s in mixed) ≈
          sum(s.mixings[end].unmixed_GeV for s in mixed) atol = 1e-10
    @test mixed[1].mixings[end].partner_masses_GeV == mixed[2].mixings[end].partner_masses_GeV
    # lookup by quantum numbers
    s = spectrum_state(spec, 1, "P", 3, 2)
    @test s.label == "1^3P_2"
    @test_throws ArgumentError spectrum_state(spec, 3, "S", 1, 0)
    # level beyond the per-channel budget fails loudly
    @test_throws ArgumentError compute_spectrum(
        params, meson;
        levels = [BasisState(7, "S", 1, 0)], ngrid = 80, rmax = 8.0,
    )
end

@testset "compute_spectrum same-J mixing gated by flavor content" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    pair = [BasisState(1, "P", 1, 1), BasisState(1, "P", 3, 1)]
    us = compute_spectrum(params, Meson(mq, :u, :s); levels = pair, ngrid = 120, rmax = 12.0)
    so_mixed = [s for s in us.states if any(m.mechanism == "antisymmetric_spin_orbit" for m in s.mixings)]
    @test length(so_mixed) == 2
    @test all(s.fine_structure_mass_convention == "unequal_mass_same_j_mixed" for s in so_mixed)
    @test sum(s.mass_GeV for s in so_mixed) ≈
          sum(s.mixings[end].unmixed_GeV for s in so_mixed) atol = 1e-10
    @test any(abs(s.mixings[end].offdiag_GeV) > 0 for s in so_mixed)
    # equal flavor: the antisymmetric matrix element vanishes, no block forms
    cc = compute_spectrum(params, Meson(mq, :c, :c); levels = pair, ngrid = 120, rmax = 12.0)
    @test all(isempty(s.mixings) for s in cc.states)
end

@testset "annihilation blocks built from two spectra (no reference data)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    ps_levels = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
    nn = compute_spectrum(params, Meson(mq, :q, :q);
        levels = ps_levels, ngrid = 120, rmax = 12.0, use_fine_structure = false)
    ss = compute_spectrum(params, Meson(mq, :s, :s);
        levels = ps_levels, ngrid = 120, rmax = 12.0, use_fine_structure = false)

    # the calibrated rank-one block reproduces its targets by construction
    targets = (0.520, 0.960, 1.440, 1.630)
    cal = pseudoscalar_annihilation_block(CalibratedP1Annihilation(), params, nn, ss;
        targets = targets)
    @test sort(cal.masses) ≈ sort(collect(targets)) atol = 1e-10
    @test all(cal.weights_GeV .>= 0.0)
    # explicit targets are mandatory - the digitized values live in GIPaper
    @test_throws ArgumentError pseudoscalar_annihilation_block(
        CalibratedP1Annihilation(), params, nn, ss)

    p1 = pseudoscalar_annihilation_block(PaperP1Annihilation(), params, nn, ss)
    p2 = pseudoscalar_annihilation_block(PaperP2Annihilation(), params, nn, ss)
    @test all(isfinite, p1.masses)
    @test all(isfinite, p2.masses)
    @test p1.masses != p2.masses
    # light nn̄ inputs carry the sqrt(2) flavor-coherence label
    input = annihilation_basis_input(nn, ps_levels[1])
    @test input.label == "1 ns"
    @test annihilation_basis_input(ss, ps_levels[2]).label == "2 ss"

    # general Eq. (16) block for one channel across the two flavors
    s1 = BasisState(1, "S", 3, 1)
    nn3 = compute_spectrum(params, Meson(mq, :q, :q);
        levels = [s1], ngrid = 120, rmax = 12.0, use_fine_structure = false)
    ss3 = compute_spectrum(params, Meson(mq, :s, :s);
        levels = [s1], ngrid = 120, rmax = 12.0, use_fine_structure = false)
    block = isoscalar_annihilation_block(params, nn3, ss3, s1;
        amplitude_A = params.annihilation.s1_A)
    @test length(block.masses) == 2
    @test all(isfinite, block.masses)
    # trace conservation: eigenvalue sum equals diagonal sum plus block trace
    diag_sum = spectrum_state(nn3, s1).mass_GeV + spectrum_state(ss3, s1).mass_GeV
    @test sum(block.masses) ≈ diag_sum + sum(diag(block.annihilation_matrix_GeV)) atol = 1e-10
end

@testset "compute_spectrum builds phase-fixed HO wave caches" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    spec = compute_spectrum(params, Meson(mq, :q, :q);
        levels = [BasisState(1, "S", 1, 0)], ngrid = 80, rmax = 8.0)
    key = RadialChannelKey(spec.meson.constituent_masses, "S")
    @test haskey(spec.computation.ho_wave_cache, key)
    ho_sol = spec.computation.ho_wave_cache[key]
    # GI annihilation phase convention: Φ(0) ∝ ∫ r u(r) dr > 0 for every level
    phase(sol, n) = sum(sol.r .* view(sol.eigenvectors, :, n))
    @test ho_sol.eigenvectors[1, 1] > 0
    @test phase(ho_sol, 1) > 0
    @test phase(ho_sol, 2) > 0
    # opt out of the HO basis entirely
    fd_only = compute_spectrum(params, Meson(mq, :q, :q);
        levels = [BasisState(1, "S", 1, 0)], ngrid = 80, rmax = 8.0,
        annihilation_wave_basis = :fd)
    @test isempty(fd_only.computation.ho_wave_cache)
end

@testset "staged spectrum: central -> corrected -> mixed" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    us = Meson(mq, :u, :s)
    levels = spectrum_levels(2; L_labels = ("S", "P"))
    kwargs = (ngrid = 250, rmax = 16.0)

    central = central_spectrum(params, us; levels = levels, kwargs...)
    @test central isa CentralSpectrum
    @test eltype(central.states) === CentralState
    @test parameters(central) === central.computation.params
    @test all(isfinite(s.central_GeV) for s in central.states)

    corrected = add_spin_corrections(central)
    @test corrected isa CorrectedSpectrum
    @test corrected.computation === central.computation
    for s in corrected.states
        # property forwarding into the wrapped central state
        @test s.central_GeV == s.central.central_GeV
        @test s.mass_GeV ≈ s.central_GeV + s.contact_shift_GeV + s.fine_structure_shift_GeV
    end

    mixed = add_intra_meson_mixing(corrected)
    @test mixed isa MixedSpectrum
    @test mixed.computation === central.computation

    # composition is compute_spectrum
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
        add_spin_corrections(central; contact_hyperfine = false, use_fine_structure = false),
    )
    @test all(s.mass_GeV == s.central_GeV for s in bare.states)
    @test all(isempty(s.mixings) for s in bare.states)
    @test all(s.fine_structure_mass_convention == "disabled" for s in bare.states)

    # stage order is enforced by dispatch
    @test !hasmethod(add_spin_corrections, Tuple{MixedSpectrum})
    @test !hasmethod(add_intra_meson_mixing, Tuple{CentralSpectrum})
end

@testset "nonperturbative contact states expose hyperfine-distinct S waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = Meson(mq, :q, :q).constituent_masses
    central = central_spectrum(params, Meson(mq, :q, :q); levels = spectrum_levels(1))
    r = central.computation.channel_cache[RadialChannelKey(masses, "S")].r
    h = r[2] - r[1]

    lvl1, vec1, r1 = contact_hyperfine_nonperturbative_states(params, masses, "S", 1, r, 2)
    lvl3, vec3, r3 = contact_hyperfine_nonperturbative_states(params, masses, "S", 3, r, 2)
    @test r1 == r3 == r                     # shared mesh with the central solve
    @test size(vec1, 2) >= 1 && size(vec3, 2) >= 1
    # levels agree with the energy-only accessor
    @test lvl1 ≈ GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 1, r, 2)
    @test lvl3 ≈ GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 3, r, 2)

    # ^1S_0 (pi) is more compact than ^3S_1 (rho): smaller <r^2>, lower energy
    r2(u) = radial_cross_expect_udr(u, u, r, h, (x, _i) -> x^2)
    @test r2(vec1[:, 1]) < r2(vec3[:, 1])
    @test lvl1[1] < lvl3[1]

    # A basis with no implementation now throws. This reverses an earlier
    # deliberate choice ("inactive path (non-FD basis) returns empties, not an
    # error"), because "inactive" conflated two different things: a
    # configuration the contact term genuinely does not touch (non-S wave,
    # sandwich off -- still empty, asserted above) and a basis nobody wrote the
    # solve for. Callers read empty as "not available" and substitute
    # first-order PT, which for the light 1S0 gives 0.2842 GeV against the
    # resummed 0.0950 GeV. No caller relied on the empty-for-HO return.
    @test_throws ArgumentError contact_hyperfine_nonperturbative_states(
        with_basis(params, HarmonicOscillatorBasis), masses, "S", 1, r, 2)
end

@testset "Table VII gluonic annihilation (Eq. 17 S_L)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mb = mq["b"]
    masses = Meson(mq, :b, :b).constituent_masses
    tomev(a) = abs(a) * sqrt(1000)   # GeV^1/2 -> MeV^1/2

    # bottomonium (most paper-faithful) S-wave: eta_b / Upsilon central wave
    vals, vecs, r = channel_solution(params, masses, 0; nlevels = 2, ngrid = 900, rmax = 24.0)
    S0 = wavefunction_origin_smearing(RadialWaveOnUniformMesh(vecs[:, 1], r), mb; L = 0)
    a0 = GIModel.alpha_s_q(vals[1])
    # zero-parameter amplitudes vs paper (eta_b -> 2g = 2.5, Upsilon -> 3g = 0.21)
    @test 0.85 < tomev(gluonic_annihilation_amplitude(:S0_2g, S0, a0, mb)) / 2.5  < 1.15
    @test 0.85 < tomev(gluonic_annihilation_amplitude(:S1_3g, S0, a0, mb)) / 0.21 < 1.25

    # width is amplitude squared, for every channel
    for ch in GLUONIC_CHANNELS
        @test gluonic_annihilation_width(ch, S0, a0, mb) ≈
              gluonic_annihilation_amplitude(ch, S0, a0, mb)^2
    end

    # P-wave chi_2b via S1 (paper chi_2b -> 2g = 0.35)
    valsP, vecsP, rP = channel_solution(params, masses, 1; nlevels = 1, ngrid = 900, rmax = 24.0)
    S1 = wavefunction_origin_smearing(RadialWaveOnUniformMesh(vecsP[:, 1], rP), mb; L = 1)
    @test 0.8 < tomev(gluonic_annihilation_amplitude(:P2_2g, S1, GIModel.alpha_s_q(valsP[1]), mb)) / 0.35 < 1.2

    # S_L is normalization-invariant (the wave is renormalized internally)
    S0_scaled = wavefunction_origin_smearing(RadialWaveOnUniformMesh(3.0 .* vecs[:, 1], r), mb; L = 0)
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
        _, _, r = channel_solution(params, masses, 0; nlevels = 2, ngrid = 1000, rmax = 24.0)
        lv, vec, r2 = contact_hyperfine_nonperturbative_states(params, masses, "S", 3, r, 2)
        return lv[n], RadialWaveOnUniformMesh(outer!(copy(vec[:, n])), r2)
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
    w_pi = RadialWaveOnUniformMesh(u_pi, r); w_rho = RadialWaveOnUniformMesh(u_rho, r)

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
        w_pi, RadialWaveOnUniformMesh(u_pi, r .+ 1.0), 1.0, 1.0, m)
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
    _, _, r = channel_solution(params, ConstituentMasses(mc, mc), 0; nlevels = 2, ngrid = 1000, rmax = 24.0)
    lv, vec, r2 = contact_hyperfine_nonperturbative_states(params, ConstituentMasses(mc, mc), "S", 1, r, 2)
    wηc = RadialWaveOnUniformMesh(outer!(copy(vec[:, 1])), r2)
    Aηc = tokeV(two_photon_amplitude(:P, wηc, mc, lv[1], 4 / 9))
    @test 0.85 < Aηc / 2.6 < 1.25

    # A2 -> gamma gamma: ³P₂ light isovector, q_eff = (e_u²−e_d²)/√2, paper −1.2 keV^½
    mqk = mq["q"]
    valsP, vecsP, rP = channel_solution(params, ConstituentMasses(mqk, mqk), 1; nlevels = 1, ngrid = 1000, rmax = 24.0)
    wA2 = RadialWaveOnUniformMesh(outer!(copy(vecsP[:, 1])), rP)
    AA2 = tokeV(two_photon_amplitude(:P2, wA2, mqk, valsP[1], (4 / 9 - 1 / 9) / sqrt(2)))
    @test AA2 < 0                       # −√(4/5) prefactor => negative amplitude
    @test 0.8 < abs(AA2) / 1.2 < 1.2

    @test_throws ArgumentError two_photon_amplitude(:bogus, wηc, mc, lv[1], 4 / 9)
end

@testset "Table VII charge radii (part d)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mu, ms = mq["q"], mq["s"]

    function ps_wave(m1, m2)
        _, _, r = channel_solution(params, ConstituentMasses(m1, m2), 0; nlevels = 2, ngrid = 1200, rmax = 26.0)
        lv, vec, r2 = contact_hyperfine_nonperturbative_states(params, ConstituentMasses(m1, m2), "S", 1, r, 2)
        return RadialWaveOnUniformMesh(vec[:, 1], r2)
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
    @test charge_radius_squared(RadialWaveOnUniformMesh(3.0 .* w.u, w.r), mu, 2 / 3, mu, 1 / 3) ≈
          charge_radius_squared(w, mu, 2 / 3, mu, 1 / 3)
end

@testset "Table VII mixed eta/eta' two-photon (P1 coherent sum)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    mu, ms = mq["q"], mq["s"]
    outer!(u) = (pk = maximum(abs, u); i = findlast(x -> abs(x) > 0.2pk, u);
                 (i !== nothing && u[i] < 0) && (u .*= -1); u)
    Qnn = (4 / 9 + 1 / 9) / sqrt(2); Qss = 1 / 9

    function psfam(m1, m2)
        _, _, r = channel_solution(params, ConstituentMasses(m1, m2), 0; nlevels = 3, ngrid = 1000, rmax = 24.0)
        lv, vec, r2 = contact_hyperfine_nonperturbative_states(params, ConstituentMasses(m1, m2), "S", 1, r, 3)
        return [RadialWaveOnUniformMesh(outer!(copy(vec[:, n])), r2) for n in 1:2]
    end
    NN, SS = psfam(mu, mu), psfam(ms, ms)
    comp = [(NN[1], mu, Qnn), (SS[1], ms, Qss), (NN[2], mu, Qnn), (SS[2], ms, Qss)]

    psl = [GIModel.BasisState(1, "S", 1, 0), GIModel.BasisState(2, "S", 1, 0)]
    vl = [GIModel.BasisState(1, "S", 3, 1)]
    nn_spec = GIModel.compute_spectrum(params, Meson(:q, :q, ConstituentMasses(mu, mu)); levels = vcat(psl, vl))
    ss_spec = GIModel.compute_spectrum(params, Meson(:s, :s, ConstituentMasses(ms, ms)); levels = vcat(psl, vl))
    psb = GIModel.pseudoscalar_annihilation_block(GIModel.PaperP1Annihilation(), params, nn_spec, ss_spec)

    ggamp(col, Mphys) = sum(psb.vectors[k, col] *
        two_photon_amplitude(:P, comp[k][1], comp[k][2], Mphys, comp[k][3]) for k in 1:4) * sqrt(1e6)
    Aη = ggamp(1, 0.548)
    Aη′ = ggamp(2, 0.958)
    # both positive, and the η<η' ordering that IDEAL mixing gets backwards
    @test Aη > 0 && Aη′ > 0
    @test abs(Aη) < abs(Aη′)
    # magnitudes in the paper's ballpark (η≈0.5, η'≈1.3 keV^½)
    @test 0.6 < abs(Aη) / 0.5 < 1.6
    @test 0.6 < abs(Aη′) / 1.3 < 1.4
end

@testset "W6 paper-order spin-distorted waves (HO first order)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    params_ho = with_basis(params, HarmonicOscillatorBasis)
    mc = mq["c"]
    masses = ConstituentMasses(mc, mc)
    outer!(u) = (pk = maximum(abs, u); i = findlast(x -> abs(x) > 0.2pk, u);
                 (i !== nothing && u[i] < 0) && (u .*= -1); u)
    ngrid, rmax = 900, 24.0
    r, h = GIModel.radial_grid(ngrid, rmax)
    tomev(a) = abs(a) * sqrt(1000)
    sm(u, rr, L) = wavefunction_origin_smearing(
        RadialWaveOnUniformMesh(outer!(copy(u)), rr), mc; L = L)

    # 1. basis-fidelity control: HO central S_L matches FD to <2% for charm —
    #    the 15-20% gluonic row residuals were never a basis artifact
    _, fdv, fdr = channel_solution(params, masses, 0; nlevels = 2, ngrid = ngrid, rmax = rmax)
    _, hov, hor = channel_solution(params_ho, masses, 0; nlevels = 2, ngrid = ngrid, rmax = rmax)
    @test 0.98 < sm(hov[:, 1], hor, 0) / sm(fdv[:, 1], fdr, 0) < 1.02

    # 2. paper-order treatment: first-order PT in the HO central eigenbasis with
    #    the calibrated spin blocks lands the charm gluonic rows on the paper
    V1 = GIModel.contact_hyperfine_operator(params, masses, "S", 1, r)
    V3 = GIModel.contact_hyperfine_operator(params, masses, "S", 3, r)
    VP0 = fine_structure_grid_operator(params, masses, 0, r, h)
    VP2 = fine_structure_grid_operator(params, masses, 2, r, h)
    amp(ch, S, M) = tomev(gluonic_annihilation_amplitude(ch, S, GIModel.alpha_s_q(M), mc))
    ratios = Dict{Symbol,Float64}()
    for (key, L, V, ch, paper) in ((:eta_c, 0, V1, :S0_2g, 4.700), (:psi, 0, V3, :S1_3g, 0.420),
                                   (:chi_0c, 1, VP0, :P0_2g, 2.500), (:chi_2c, 1, VP2, :P2_2g, 0.880))
        vals, waves, rr = ho_first_order_distorted_states(params_ho, masses, L, V;
            nlevels = 4, ngrid = ngrid, rmax = rmax)
        ratios[key] = amp(ch, sm(waves[:, 1], rr, L), vals[1]) / paper
        @test 0.85 < ratios[key] < 1.15
    end

    # 3. the splitting patterns the central wave misses collapse at paper order:
    #    central waves give ratio-of-ratios eta_c/psi ≈ 0.78, chi_0c/chi_2c ≈ 0.68
    Sc0 = sm(fdv[:, 1], fdr, 0)
    Mc0 = channel_solution(params, masses, 0; nlevels = 1, ngrid = ngrid, rmax = rmax)[1][1]
    central_eta_psi = (amp(:S0_2g, Sc0, Mc0) / 4.700) / (amp(:S1_3g, Sc0, Mc0) / 0.420)
    @test central_eta_psi < 0.85
    @test 0.90 < ratios[:eta_c] / ratios[:psi] < 1.10
    @test ratios[:chi_0c] / ratios[:chi_2c] > 0.80   # central-wave value ≈ 0.68
end

@testset "W6 paper-order full diagonalization (light-mass discriminator)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    params_ho = with_basis(params, HarmonicOscillatorBasis)
    ngrid, rmax = 900, 24.0
    r, h = GIModel.radial_grid(ngrid, rmax)
    outer!(u) = (pk = maximum(abs, u); i = findlast(x -> abs(x) > 0.2pk, u);
                 (i !== nothing && u[i] < 0) && (u .*= -1); u)

    # 1. The paper's treatment is FULL diagonalization in the finite HO basis,
    #    not first-order PT. The light ¹S₀ (pion) mass discriminates: full-diag
    #    keeps it resummed (≈0.10 GeV) like the fine-grid FD, while first-order
    #    PT over-raises it (≈0.28 GeV).
    nn = ConstituentMasses(mq["q"], mq["q"])
    Vpi = GIModel.contact_hyperfine_operator(params, nn, "S", 1, r)
    m_full = ho_full_distorted_states(params_ho, nn, 0, Vpi; nlevels = 4, ngrid = ngrid, rmax = rmax)[1][1]
    m_pt = ho_first_order_distorted_states(params_ho, nn, 0, Vpi; nlevels = 4, ngrid = ngrid, rmax = rmax)[1][1]
    m_fd = GIModel.lowest_eigenpairs(
        Symmetric(Matrix(GIModel.relativistic_hamiltonian(params, nn, 0; ngrid = ngrid, rmax = rmax)[1]) + Matrix(Vpi)), 1)[1][1]
    @test m_full < 0.15               # resummed, light pion
    @test abs(m_full - m_fd) < 0.02   # matches the fine-grid FD resummation
    @test m_pt > 0.20                 # first-order PT over-raises it

    # 2. full diagonalization still lands the charm gluonic singlet on the paper
    #    (the spin distortion the central wave misses), so both subtable regimes
    #    are served by ONE treatment.
    mc = mq["c"]; cc = ConstituentMasses(mc, mc)
    V1 = GIModel.contact_hyperfine_operator(params, cc, "S", 1, r)
    v, w, rr = ho_full_distorted_states(params_ho, cc, 0, V1; nlevels = 4, ngrid = ngrid, rmax = rmax)
    S = wavefunction_origin_smearing(RadialWaveOnUniformMesh(outer!(copy(w[:, 1])), rr), mc; L = 0)
    eta_c = abs(gluonic_annihilation_amplitude(:S0_2g, S, GIModel.alpha_s_q(v[1]), mc)) * sqrt(1000) / 4.700
    @test 0.95 < eta_c < 1.20
end
