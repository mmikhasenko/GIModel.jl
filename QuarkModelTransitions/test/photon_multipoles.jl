const QMTm = GIModel.QuarkModelTransitions

# Pure component state with an explicit wave and mass.
mstate(label, flavors, L, S, J, wave, M) = PhysicalState(label, M, [(
    basis = BasisState(1, L, 2S + 1, J; flavors = flavors),
    coefficient = 1.0, wave = wave,
)])
fraction(amplitude, name) = multipole_fractions(amplitude)[name]
amp_of(amplitude, name) = only(m.amplitude for m in amplitude.multipoles if m.multipole == name)

@testset "Photon multipoles: Siegert's theorem" begin
    # Exact NR oscillator, H = p^2/2mu + mu Omega^2 r^2/2, beta^2 = mu Omega: the
    # convection current's electric multipoles equal -sqrt((k+1)/2k) (omega/q)
    # times the charge-density multipoles, omega = E_i - E_f = k Omega.
    m1, m2, Omega = 0.4, 1.7, 0.3
    beta = sqrt(m1 * m2 / (m1 + m2) * Omega)
    masses = QuarkMassTable("c" => m1, "b" => m2)
    op = MultipolePhotonEmission(masses)
    final = mstate("1S", (:c, :b), "S", 0, 0, OscillatorWave(0, beta, [1.0]), 1.0)
    for (k, L) in ((1, "P"), (2, "D"), (3, "F"))
        initial = mstate("1$L", (:c, :b), L, 0, k, OscillatorWave(k, beta, [1.0]), 2.0)
        q = 1e-3
        conv, _, rho = QMTm._photon_component_sums(final, op, initial, q)
        current = QMTm._multipole_projection(conv, k, 0, 1)[k]
        density = QMTm._multipole_projection(rho, k, 0, 0)[k]
        @test current ≈ QMTm._siegert_factor(k) * (k * Omega / q) * density rtol = 1e-5
    end
end

@testset "Photon multipoles: M1 and E1 closed forms" begin
    masses = QuarkMassTable("c" => 1.628, "b" => 4.977)
    op = MultipolePhotonEmission(masses)
    us, ut = OscillatorWave(0, 0.70, [1.0, 0.2]), OscillatorWave(0, 0.66, [1.0, 0.15])
    up = OscillatorWave(1, 0.55, [1.0, -0.1])
    ud = OscillatorWave(2, 0.50, [1.0])
    ov(a, b, f) = radial_overlap(a, b, f)
    j0(x) = QMTm._spherical_bessel_j(0, x)
    for (f1, f2) in ((:c, :c), (:c, :b))
        m = (masses[string(f1)], masses[string(f2)]); M = sum(m)
        e = (f1 == :c ? 2/3 : -1/3, f2 == :c ? 2/3 : -1/3)
        # M1, exact recoil j0(k m_j r / M) for each constituent (Godfrey-Moats)
        V = mstate("V", (f1, f2), "S", 1, 1, ut, M + 0.2)
        P = mstate("P", (f1, f2), "S", 0, 0, us, M + 0.1)
        a = matrix_element(P, op, V)
        q = a.momentum_GeV
        mu = e[1] / m[1] * ov(us, ut, r -> j0(m[2] / M * q * r)) +
             e[2] / m[2] * ov(us, ut, r -> j0(m[1] / M * q * r))
        @test decay_width(a) ≈ 1000 * QMTm.ALPHA_EM / 3 * q^3 * mu^2 rtol = 1e-9
        @test [x.multipole for x in a.multipoles] == [:M1]
        # E1 with the centre-of-mass charge at small q (magnetization adds O(q/m))
        effective = (m[2] * e[1] + m[1] * e[2]) / M
        for (initial, final) in ((mstate("P2", (f1, f2), "P", 1, 2, up, M + 1e-4), V),
                                 (mstate("D1", (f1, f2), "D", 1, 1, ud, M + 1e-4),
                                  mstate("P0", (f1, f2), "P", 1, 0, up, M)))
            b = QMTm._multipole_matrix_element(final, op, initial, 1e-4)
            Li = GIModel.orbital_angular_momentum(initial.components[1].basis.L_label)
            Lf = GIModel.orbital_angular_momentum(final.components[1].basis.L_label)
            R = QMTm._photon_angular_rate(Lf, 1, final.J, Li, 1, initial.J, QMTm._E1_OPERATORS)
            wf, wi = final.components[1].wave, initial.components[1].wave
            standard = 1000 * 4 / 3 * QMTm.ALPHA_EM * 1e-12 * effective^2 * R * ov(wf, wi, r -> r)^2
            @test amp_of(b, :E1)^2 ≈ standard rtol = 1e-3
        end
    end
end

@testset "Photon multipoles: magnetization E1 and Karl-Meshkov-Rosner M2/E1" begin
    mc = 1.628
    masses = QuarkMassTable("c" => mc)
    op = MultipolePhotonEmission(masses)
    lead = PhotonEmission(masses; electric_exponent = 0.0)
    s1, s2 = OscillatorWave(0, 0.7, [1.0]), OscillatorWave(0, 0.7, [0.0, 1.0])
    pw = OscillatorWave(1, 0.6, [1.0])
    q = 1e-3
    psi = mstate("psi", (:c, :c), "S", 1, 1, s1, 3.0)
    psi2 = mstate("psi2", (:c, :c), "S", 1, 1, s2, 3.0)
    for J in 0:2
        chi = mstate("chi$J", (:c, :c), "P", 1, J, pw, 3.0)
        a = QMTm._multipole_matrix_element(psi, op, chi, q)
        leading = QMTm._photon_matrix_element(psi, lead, chi, q)
        # E1 rate times 1 + c_J q/(2m) with c_J = -<L.S> = 2 - J(J+1)/2
        @test amp_of(a, :E1)^2 / abs2(leading.value) - 1 ≈ (2 - J * (J + 1) / 2) * q / (2mc) rtol = 1e-2
        if J > 0
            kmr = J == 1 ? -1.0 : -3 / sqrt(5)
            @test fraction(a, :M2) / (q / (4mc)) ≈ kmr rtol = 1e-3
            # psi(2S) -> gamma chi_cJ: CLEO's b_2 = -(this convention's M2 fraction)
            b = QMTm._multipole_matrix_element(chi, op, psi2, q)
            @test -fraction(b, :M2) / (q / (4mc)) ≈ -kmr rtol = 1e-3
        end
        J == 2 && @test amp_of(a, :E3) == 0.0     # no E3 between pure P and S waves
    end
end

@testset "Photon multipoles: spin flip agrees with PhotonEmission" begin
    # Equal and unequal masses: validates the sqrt(120)/sqrt(72) denominators and
    # the centre-of-mass emitter weights of the leading kernel.
    masses = QuarkMassTable("q" => 0.22, "s" => 0.419, "c" => 1.628, "b" => 4.977)
    op = MultipolePhotonEmission(masses)
    lead = PhotonEmission(masses; electric_exponent = 0.0)
    sw, pw = OscillatorWave(0, 0.5, [1.0, 0.1]), OscillatorWave(1, 0.45, [1.0])
    q = 1e-3
    for flavors in ((:u, :d), (:u, :s), (:c, :s), (:u, :b)), J in (1, 2)
        P = mstate("P", flavors, "P", 1, J, pw, 3.0)
        S = mstate("S", flavors, "S", 0, 0, sw, 2.0)
        exact = decay_width(QMTm._multipole_matrix_element(S, op, P, q))
        leading = decay_width(QMTm._photon_matrix_element(S, lead, P, q))
        @test exact ≈ leading rtol = 1e-3
    end
end

@testset "Photon multipoles: representations and selection rules" begin
    masses = QuarkMassTable("c" => 1.628)
    op = MultipolePhotonEmission(masses)
    sw, pw = OscillatorWave(0, 0.7, [1.0, 0.1]), OscillatorWave(1, 0.6, [1.0])
    r = collect(range(0.02, 30.0; length = 1500))
    mesh(w) = sample_wave(w, r)
    chi = mstate("chi", (:c, :c), "P", 1, 2, pw, 3.556)
    psi = mstate("psi", (:c, :c), "S", 1, 1, sw, 3.097)
    a = matrix_element(psi, op, chi)
    b = matrix_element(mstate("psi", (:c, :c), "S", 1, 1, mesh(sw), 3.097), op,
                       mstate("chi", (:c, :c), "P", 1, 2, mesh(pw), 3.556))
    @test decay_width(b) ≈ decay_width(a) rtol = 1e-3
    @test fraction(b, :M2) ≈ fraction(a, :M2) rtol = 1e-3
    # C parity: c cbar 3S1 -> 1S0 has M1 and no E2
    eta = mstate("eta", (:c, :c), "S", 0, 0, sw, 2.98)
    @test [m.multipole for m in matrix_element(eta, op, psi).multipoles] == [:M1]
    @test_throws ArgumentError matrix_element(mstate("x", (:c, :c), "S", 0, 0, sw, 3.5), op, eta)
    @test_throws ArgumentError MultipolePhotonEmission(masses; electric_form = :velocity)
    @test decay_width(chi, op, psi) == 0.0
end
