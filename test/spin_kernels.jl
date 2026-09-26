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

@testset "public Appendix-A fine-structure radial kernels" begin
    params, mq = load_parameters_and_quark_masses(
        joinpath(root, "data", "parameters.provisional.toml"),
    )
    masses = ConstituentMasses(mq["q"], mq["c"])
    swapped = ConstituentMasses(mq["c"], mq["q"])
    for r in (0.05, 0.4, 1.5, 4.0)
        k = fine_structure_radial_kernels(params, masses, r)
        ks = fine_structure_radial_kernels(params, swapped, r)
        @test all(isfinite, values(k))
        @test k.vector_11 == ks.vector_22
        @test k.vector_22 == ks.vector_11
        @test k.vector_12 == ks.vector_12
        @test k.scalar_11 == ks.scalar_22
        @test k.scalar_22 == ks.scalar_11
        @test k.tensor_12 == ks.tensor_12
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

@testset "smeared contact acts in L>0 states (GI Eqs. 23-26)" begin
    # A15 has no L restriction and GI Eqs. (23)-(26) keep the contact S in every
    # P-wave level: ^1P_1 gets -3S/4, ^3P_J get +S/4, with S "small but nonzero
    # due to relativistic smearing". Omitting it put the ^1P_1/^3P_1 diagonal gap
    # of every unequal-mass P wave ~20-40 MeV too high.
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = ConstituentMasses(mq["c"], mq["u"])
    only_contact(on) = SpinTerms(contact_hyperfine = on, fine_structure = false,
        same_j_spin_orbit = false, tensor = false)
    for solver in (FiniteDifferenceSolver(), OscillatorSolver())
        shift(mult) = begin
            multiplet = FineStructureMultiplet("P", mult, 1)
            on = fixed_channel_solution(params, masses, multiplet;
                solver, terms = only_contact(true), nlevels = 1)
            off = fixed_channel_solution(params, masses, multiplet;
                solver, terms = only_contact(false), nlevels = 1)
            first_order = GIModel.contact_hyperfine_shift_active(params, masses, multiplet,
                radial_wave(off, 1))
            (on.eigenvalues_GeV[1] - off.eigenvalues_GeV[1], first_order)
        end
        (singlet, singlet_first), (triplet, triplet_first) = shift(1), shift(3)
        S = -4singlet_first / 3
        @test 0.020 < S < 0.030                     # c ubar 1P: S ≈ 25 MeV
        @test triplet_first ≈ S / 4 rtol = 0.03     # same S, spin factor +1/4
        @test singlet ≈ singlet_first rtol = 0.05   # resummed ≈ first order
        @test triplet ≈ triplet_first rtol = 0.05
    end
    # Independent Python evaluation (numerical A7-A8 smearing, momentum factors
    # applied by a j1 Hankel transform) on the production FD 1^1P_1 wave
    # (contact included in its fixed-sector Hamiltonian): S = 25.064 MeV.
    spec = fixed_spectrum(params, Meson(mq, :c, :u); levels = [BasisState(1, "P", 1, 1)])
    wave = radial_wave(spec, spectrum_state(spec, "1^1P_1"))
    S = -4GIModel.contact_hyperfine_shift_active(params, masses, FineStructureMultiplet("P", 1, 1), wave) / 3
    @test 1000S ≈ 25.064 atol = 0.02
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
