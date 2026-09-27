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
    @test_throws ArgumentError neutral_m1_charge(:s; isovector_left = true)

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
    @test_throws ArgumentError m1_recoil_moment(w1, w2, m, 1.0, -0.1)

    @test photon_recoil_form_factor(0.0) == 1.0
    @test photon_recoil_form_factor(1.6; beta = 0.4) ≈ exp(-1)
    @test_throws ArgumentError photon_recoil_form_factor(1.0; beta = 0.0)
    @test m1_radiative_width(1.0, 0.0) == 0.0
    @test m1_radiative_width(-1.0, 0.1) == m1_radiative_width(1.0, 0.1)
    @test m1_radiative_width(1.0, 0.2) ≈ 8m1_radiative_width(1.0, 0.1)
    @test m1_radiative_width(1.0, 0.1; parent_spin = 0) ≈
          3m1_radiative_width(1.0, 0.1; parent_spin = 1)
    @test_throws ArgumentError m1_radiative_width(1.0, -0.1)
end

@testset "Four-flavor physical-state annihilation bridge" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    solver = FiniteDifferenceSolver(ngrid = 90, rmax = 12.0, nlevels_per_channel = 2)
    flavors = [(:q, :q), (:s, :s), (:c, :c), (:b, :b)]
    levels = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0),
              BasisState(1, "S", 3, 1)]
    spectra = [compute_spectrum(params, Meson(mq, f...);
        levels = levels, solver = solver) for f in flavors]
    basis = [BasisState(n, "S", 1, 0; flavors = f)
             for n in 1:2 for f in flavors if !(n == 2 && f == (:b, :b))]
    amplitudes = Dict(("S", 3, 1) => params.annihilation.s1_A)
    final = add_isoscalar_annihilation(params, spectra;
        pseudoscalar = PaperP1Annihilation(), pseudoscalar_basis = basis,
        amplitudes = amplitudes, annihilation_radial_levels = (1,))
    eta = spectrum_state(final, BasisState(1, "S", 1, 0; flavors = (:q, :q)))
    components = physical_components(final, eta)
    @test length(components) == 7
    @test Set(c.basis.flavors for c in components) == Set(flavors)
    @test sum(abs2(c.coefficient) for c in components) ≈ 1.0
    @test all(c.wave isa MeshWave for c in components)
    @test all(abs(c.coefficient) > 0 for c in components)
    @test all(isempty(s.mixings) for spec in spectra for s in spec.states)
    # Exactly the same two-channel model through the old and new interfaces.
    old = add_isoscalar_annihilation(params, spectra[1], spectra[2]; amplitudes = amplitudes)
    new = add_isoscalar_annihilation(params, spectra[1:2]; amplitudes = amplitudes)
    @test [s.mass_GeV for s in old.states] ≈ [s.mass_GeV for s in new.states]
    # The final state's coefficient composition agrees with the independently
    # assembled seven-channel kernel, including charm and bottom terms.
    inputs = [annihilation_basis_input(spectra[findfirst(==(b.flavors), flavors)], b)
              for b in basis]
    expected = isoscalar_pseudoscalar_annihilation_solution(
        PaperP1Annihilation(), params, inputs)
    actual = last(eta.mixings).result
    @test actual.block.matrix ≈ expected.block.matrix
    @test eta.mass_GeV ≈ first(expected.masses)
    for state in final.states
        isempty(state.mixings) && continue
        components = physical_components(final, state)
        precursor = only(c for c in components if c.basis == state.basis)
        @test precursor.coefficient > 0
    end
    @test_throws ArgumentError add_isoscalar_annihilation(params, [spectra[1], spectra[1]])
    @test_throws ArgumentError add_isoscalar_annihilation(params, spectra;
        pseudoscalar = PaperP1Annihilation())
    @test_throws ArgumentError add_isoscalar_annihilation(params, spectra;
        pseudoscalar_basis = basis)
end

@testset "Photon multipole angular factors and powers" begin
    @test e1_angular_coefficient(2; parent_is_S = true)^2 /
          e1_angular_coefficient(2)^2 ≈ 5 / 3
    @test e1_angular_coefficient(0; parent_is_S = true)^2 /
          e1_angular_coefficient(0)^2 ≈ 1 / 3
    @test e1_angular_coefficient(1; singlet = true) == sqrt(2.0)
    @test_throws ArgumentError e1_angular_coefficient(0; singlet = true)
    sw = OscillatorWave(0, 0.5, [1.0])
    pw = OscillatorWave(1, 0.5, [1.0])
    a = spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 2, 0.2)
    b = spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 2, 0.4)
    @test b ≈ 2^2.5 * a # q^2 times sqrt(q), at fixed waves
    @test spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 1, 0.2) ≈
          sqrt(60) / 6 * a
    @test spin_flip_photon_amplitude(sw, pw, [(2/3, 0.22), (1/3, 0.22)], 2, 0.2) ≈ a
    @test_throws ArgumentError spin_flip_photon_amplitude(sw, pw, [(1.0, 0.22)], 0, 0.2)

    # Attractive general annihilation gives opposite-sign flavor admixtures;
    # each physical tensor state nevertheless keeps its dominant flavor positive.
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    inputs = [GIModel.pseudoscalar_annihilation_basis_input(
        BasisState(1, "P", 3, 2; flavors = (f, f)), m, M, pw;
        isoscalar_coherent = f == :q)
        for (f, m, M) in [(:q, mq["q"], 1.3), (:s, mq["s"], 1.5)]]
    result = isoscalar_general_annihilation_solution(params, inputs;
        amplitude_A = -0.8, L = 1, multiplicity = 3, J = 2)
    @test result.vectors[2, 2] > 0
    @test result.vectors[1, 2] < 0
    @test result.vectors' * result.vectors ≈ Matrix{Float64}(I, 2, 2)
end

@testset "PhotonEmission matrix-element API" begin
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
    @test amplitude isa RadiativeAmplitude
    @test amplitude.transition_class isa DirectM1
    @test photon_transition_class(pseudoscalar, vector) isa DirectM1
    @test amplitude.multipole == :M1
    @test amplitude.value ≈ expected
    @test length(amplitude.terms) == 2
    @test amplitude.provenance.recoil_order == 0
    @test amplitude.provenance.prescription == :hybrid_mock_meson
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
    @test e1_amplitude.transition_class isa AllowedE1
    @test e1_amplitude.value ≈ expected_e1
    @test decay_width(e1_amplitude) ≈ abs2(expected_e1)

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
    @test m2_amplitude.transition_class isa SpinFlipM2
    @test m2_amplitude.value ≈ spin_flip_photon_amplitude(
        sw, pw, [(2 / 3, mq["q"]), (1 / 3, mq["q"])], 2,
        photon_momentum(light_tensor.mass_GeV, singlet.mass_GeV),
    )

    light_axial = PhysicalState("a1", 3.40, [(
        basis = BasisState(1, "P", 3, 1; flavors = (:u, :d)),
        coefficient = 1.0, wave = pw,
    )])
    spin_flip_e1 = matrix_element(singlet, PhotonEmission(mq), light_axial)
    @test spin_flip_e1.transition_class isa SpinFlipE1
    @test spin_flip_e1.multipole == :E1

    excited_vector = PhysicalState("psi(2S)", 3.68, [(
        basis = BasisState(2, "S", 3, 1; flavors = (:c, :c)),
        coefficient = 1.0, wave = sw2,
    )])
    @test photon_transition_class(pseudoscalar, excited_vector) isa HinderedM1
    recoil = PhotonEmission(mq; recoil_order = 2)
    recoil_amplitude = matrix_element(pseudoscalar, recoil, excited_vector)
    @test recoil_amplitude.transition_class isa HinderedM1
    @test recoil_amplitude.provenance.recoil_order == 2

    unresolved = PhysicalState("qbarq", 2.0, [(
        basis = BasisState(1, "S", 3, 1; flavors = (:q, :q)),
        coefficient = 1.0, wave = sw,
    )])
    unresolved_final = PhysicalState("qbarq0", 1.0, [(
        basis = BasisState(1, "S", 1, 0; flavors = (:q, :q)),
        coefficient = 1.0, wave = sw,
    )])
    @test_throws ArgumentError matrix_element(unresolved_final, m1, unresolved)
    @test_throws ArgumentError matrix_element(vector, PhotonEmission(mq; recoil_order = 2), tensor)
    @test_throws ArgumentError PhotonEmission(mq; recoil_order = 1)
    @test_throws ArgumentError matrix_element(
        ReferenceState("eta_c", 2.98; J = 0, parity = -1), m1, vector,
    )
end
