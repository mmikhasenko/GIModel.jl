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
