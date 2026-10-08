@testset "Radiative charge normalization and recoil" begin
    # Normalize flavor states before applying the electromagnetic operator.
    # The Table VI perfect-mixing eta coefficients must follow from the
    # state vector, not be multiplied into both state and operator.
    @test neutral_m1_charge(:q) / sqrt(2) ≈ 1 / (3sqrt(2))
    @test neutral_m1_charge(:s) * (-1 / sqrt(2)) ≈ sqrt(2) / 3
    @test neutral_m1_charge(:q; isovector_left = true) / sqrt(2) ≈ 1 / sqrt(2)
    @test neutral_m1_charge(:q; isovector_left = true, isovector_right = true) ≈ 1 / 3
    @test neutral_m1_charge(:c) == 4 / 3
    @test neutral_m1_charge(:b) == -2 / 3

    w1 = OscillatorWave(0, 0.5, [1.0, 0.0])
    w2 = OscillatorWave(0, 0.5, [0.0, -1.0])
    m = 1.8
    mw1, mw2 = GIModel.momentum_wave(w1), GIModel.momentum_wave(w2)
    direct = m1_transition_moment(mw1, mw2, m, m, [(4 / 3, m)])
    @test m1_recoil_moment(w1, w2, m, 4 / 3, 0.0) ≈ direct
    q1, q2 = 0.1, 0.2
    r1 = m1_recoil_moment(w1, w2, m, 4 / 3, q1) - direct
    r2 = m1_recoil_moment(w1, w2, m, 4 / 3, q2) - direct
    @test abs(r1) > 1e-6
    @test r2 ≈ 4r1
    @test m1_recoil_moment(w2, w1, m, 4 / 3, q1) ≈
          m1_recoil_moment(w1, w2, m, 4 / 3, q1)

    @test photon_recoil_form_factor(0.0) == 1.0
    @test photon_recoil_form_factor(1.6; beta = 0.4) ≈ exp(-1)
    @test m1_radiative_width(1.0, 0.0) == 0.0
    @test m1_radiative_width(-1.0, 0.1) == m1_radiative_width(1.0, 0.1)
    @test m1_radiative_width(1.0, 0.2) ≈ 8m1_radiative_width(1.0, 0.1)
    @test m1_radiative_width(1.0, 0.1; parent_spin = 0) ≈
          3m1_radiative_width(1.0, 0.1; parent_spin = 1)
end

@testset "Photon multipole angular factors and powers" begin
    @test e1_angular_coefficient(2; parent_is_S = true)^2 /
          e1_angular_coefficient(2)^2 ≈ 5 / 3
    @test e1_angular_coefficient(0; parent_is_S = true)^2 /
          e1_angular_coefficient(0)^2 ≈ 1 / 3
    # E1 does not act on spin: 1P1 -> 1S0 and 3PJ -> 3S1 share 1/3.
    @test e1_angular_coefficient(1; singlet = true) ≈ e1_angular_coefficient(2)
    @test e1_angular_coefficient(1; singlet = true, parent_is_S = true)^2 /
          e1_angular_coefficient(1; singlet = true)^2 ≈ 3
    # The factors are derived from L⊗S coupling; the old hand-typed values
    # remain as oracles.
    for J in 0:2
        @test e1_angular_coefficient(J) ≈ 1 / 3
        @test e1_angular_coefficient(J; parent_is_S = true) ≈ sqrt((2J + 1) / 3) / 3
    end
    @test e1_angular_coefficient(1; singlet = true) ≈ 1 / 3
    @test e1_angular_coefficient(1; singlet = true, parent_is_S = true) ≈ sqrt(3) / 3
    @test GIModel.QuarkModelTransitions._spin_flip_denominator(2) ≈ sqrt(120)
    @test GIModel.QuarkModelTransitions._spin_flip_denominator(1) ≈ sqrt(72)
    @test_throws ArgumentError GIModel.QuarkModelTransitions._spin_flip_denominator(0)
    sw = OscillatorWave(0, 0.5, [1.0])
    pw = OscillatorWave(1, 0.5, [1.0])
    a = spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 2, 0.2)
    b = spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 2, 0.4)
    @test b ≈ 2^2.5 * a # q^2 times sqrt(q), at fixed waves
    @test spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 1, 0.2) ≈
          sqrt(120 / 72) * a
    @test spin_flip_photon_amplitude(sw, pw, [(2/3, 0.22), (1/3, 0.22)], 2, 0.2) ≈ a

end

@testset "Photon emission amplitudes and widths" begin
    mq = QuarkMassTable("q" => 0.22, "c" => 1.628, "b" => 4.977)
    sw = OscillatorWave(0, 0.5, [1.0, 0.0])
    sw2 = OscillatorWave(0, 0.5, [0.8, 0.6])
    pw = OscillatorWave(1, 0.5, [1.0])

    vector = PhysicalState("psi", 3.10, [(
        basis = BasisState(1, "S", 3, 1; flavors = (:c, :c)),
        coefficient = 1.0, wave = sw,
    )])
    pseudoscalar = PhysicalState("eta_c", 2.98, [(
        basis = BasisState(1, "S", 1, 0; flavors = (:c, :c)),
        coefficient = 1.0, wave = sw2,
    )])
    m1 = PhotonEmission(mq)
    amplitude = matrix_element(pseudoscalar, m1, vector)
    q = photon_momentum(vector.mass_GeV, pseudoscalar.mass_GeV)
    expected = m1_transition_moment(
        momentum_wave(sw2, 0), momentum_wave(sw, 0), mq["c"], mq["c"],
        [(4 / 3, mq["c"])],
    )
    @test amplitude.value ≈ expected
    @test decay_width(amplitude) ≈
          1000m1_radiative_width(expected, q; parent_spin = 1)
    @test decay_width(pseudoscalar, m1, vector) ≈ decay_width(amplitude)
    @test decay_width(vector, m1, pseudoscalar) == 0.0

    tensor = PhysicalState("chi_c2", 3.56, [(
        basis = BasisState(1, "P", 3, 2; flavors = (:c, :c)),
        coefficient = 1.0, wave = pw,
    )])
    e1 = PhotonEmission(mq)
    e1_amplitude = matrix_element(vector, e1, tensor)
    qe = photon_momentum(tensor.mass_GeV, vector.mass_GeV)
    expected_e1 = e1_transition_amplitude(
        sw, momentum_wave(sw, 0), pw, momentum_wave(pw, 1), mq["c"],
        qvalue -> (4 / 3) * e1_angular_coefficient(2) * qvalue,
        1.0, 0.0; q = qe,
    )
    @test e1_amplitude.value ≈ expected_e1
    @test decay_width(e1_amplitude) ≈ abs2(expected_e1)

    # Same radial waves and photon momentum: h_c -> eta_c gamma equals
    # chi_c2 -> J/psi gamma.
    h_c = PhysicalState("h_c", tensor.mass_GeV, [(
        basis = BasisState(1, "P", 1, 1; flavors = (:c, :c)),
        coefficient = 1.0, wave = pw,
    )])
    eta_c = PhysicalState("eta_c", vector.mass_GeV, [(
        basis = BasisState(1, "S", 1, 0; flavors = (:c, :c)),
        coefficient = 1.0, wave = sw,
    )])
    @test matrix_element(eta_c, e1, h_c).value ≈ expected_e1

    # Unequal masses: the dipole is taken about the centre of mass, giving the
    # Eichten-Quigg effective charge (e_c m_b + e_b m_c)/(m_c+m_b) for c bbar.
    # With the m/E exponent off, E1 no longer depends on the emitter mass.
    bc(label, L, J, wave, mass) = PhysicalState(label, mass, [(
        basis = BasisState(1, L, 3, J; flavors = (:c, :b)),
        coefficient = 1.0, wave = wave,
    )])
    bc_tensor, bc_vector = bc("Bc2", "P", 2, pw, 6.75), bc("Bc*", "S", 1, sw, 6.33)
    flat = PhotonEmission(mq; electric_exponent = 0.0)
    mc, mb = mq["c"], mq["b"]
    effective_charge = (2/3 * mb - 1/3 * mc) / (mc + mb)
    expected_bc = e1_transition_amplitude(
        sw, momentum_wave(sw, 0), pw, momentum_wave(pw, 1), mc,
        qvalue -> 2effective_charge * e1_angular_coefficient(2) * qvalue,
        bc_tensor.mass_GeV, bc_vector.mass_GeV; exponent = 0.0,
    )
    @test matrix_element(bc_vector, flat, bc_tensor).value ≈ expected_bc rtol = 1e-12

    singlet = PhysicalState("pi", 2.98, [(
        basis = BasisState(1, "S", 1, 0; flavors = (:u, :d)),
        coefficient = 1.0, wave = sw,
    )])
    light_tensor = PhysicalState("a2", 3.56, [(
        basis = BasisState(1, "P", 3, 2; flavors = (:u, :d)),
        coefficient = 1.0, wave = pw,
    )])
    m2 = PhotonEmission(mq)
    m2_amplitude = matrix_element(singlet, m2, light_tensor)
    @test m2_amplitude.value ≈ spin_flip_photon_amplitude(
        sw, pw, [(2 / 3, mq["q"]), (1 / 3, mq["q"])], 2,
        photon_momentum(light_tensor.mass_GeV, singlet.mass_GeV),
    )

    light_axial = PhysicalState("a1", 3.40, [(
        basis = BasisState(1, "P", 3, 1; flavors = (:u, :d)),
        coefficient = 1.0, wave = pw,
    )])
    spin_flip_e1 = matrix_element(singlet, PhotonEmission(mq), light_axial)
    @test spin_flip_e1.value ≈ spin_flip_photon_amplitude(
        sw, pw, [(2/3, mq["q"]), (1/3, mq["q"])], 1,
        photon_momentum(light_axial.mass_GeV, singlet.mass_GeV)) rtol=1e-10

    excited_vector = PhysicalState("psi(2S)", 3.68, [(
        basis = BasisState(2, "S", 3, 1; flavors = (:c, :c)),
        coefficient = 1.0, wave = sw2,
    )])
    recoil = PhotonEmission(mq; recoil_order = 2)
    recoil_amplitude = matrix_element(pseudoscalar, recoil, excited_vector)
    qr = photon_momentum(excited_vector.mass_GeV, pseudoscalar.mass_GeV)
    @test recoil_amplitude.value ≈ m1_recoil_moment(sw2, sw2, mq["c"], 4/3, qr) rtol=1e-10

end
