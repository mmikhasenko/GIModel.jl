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
    f_V = sqrt(5.55e-6 / ((4π / 3) * QuarkModelTransitions.ALPHA_EM^2 * 3.0969))
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
        L = orbital_angular_momentum(multiplet.L_label)
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
    coarse_momentum = observable_momentum_wave(trial(r_coarse), 0; npoints = 900)
    fine_momentum = observable_momentum_wave(trial(r_fine), 0; npoints = 900)
    @test last(coarse_momentum.p) ≈ min(π / 0.05, 60.0)
    @test last(fine_momentum.p) == 60.0
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

@testset "Observable handoff from GIModel waves" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    lv = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
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
    sol = channel_solution(params, ConstituentMasses(mq["q"], mq["s"]), 0; nlevels=1)
    wave = radial_wave(sol, 1)
    scaled = MeshWave(13.7 .* wave.u, wave.r)
    radius(w) = charge_radius_squared(w, mq["q"], 2//3, mq["s"], 1//3)
    @test isapprox(radius(wave), radius(scaled); rtol=1e-10)
end
