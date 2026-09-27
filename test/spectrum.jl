@testset "Isoscalar coherence factor is stated, not string-matched" begin
    r = collect(range(0.05, 6.0; length = 64))
    wave = MeshWave(exp.(-r), r)
    mk(label, coh) = pseudoscalar_annihilation_basis_input(
        BasisState(1, "S", 1, 0; label = label), 0.22, 0.9, wave;
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

    @test_throws UndefKeywordError pseudoscalar_annihilation_basis_input(
        BasisState(1, "S", 1, 0; label = "1 ns"), 0.22, 0.9, wave)
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
    spec = compute_spectrum(
        params, meson;
        levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
    )
    @test length(spec.states) == 4
    for s in spec.states
        if isempty(s.mixings)
            @test s.mass_GeV ≈ s.central_GeV + s.contact_shift_GeV + s.fine_structure_shift_GeV atol = 1e-12
        end
        @test s.fine_structure_shift_GeV ≈ s.spin_orbit_shift_GeV + s.tensor_shift_GeV atol = 1e-12
    end
    # Every requested radial state in the compatible ^3S_1/^3D_1 sectors enters
    # one tensor block; trace and the shared eigensystem are conserved.
    mixed = [s for s in spec.states if !isempty(s.mixings)]
    @test length(mixed) == 3
    @test all(m.mechanism == "tensor_mixing" for s in mixed for m in s.mixings)
    @test sum(s.mass_GeV for s in mixed) ≈
          sum(s.mixings[end].unmixed_GeV for s in mixed) atol = 1e-10
    @test all(
        s.mixings[end].partner_masses_GeV == mixed[1].mixings[end].partner_masses_GeV for
        s in mixed
    )
    @test all(s.mixings[end].result === mixed[1].mixings[end].result for s in mixed)
    @test fieldnames(StateMixing) == (:result, :eigenstate, :unmixed_GeV)
    @test all(length(physical_components(spec, s)) == 3 for s in mixed)
    components = physical_components(spec, mixed[1])
    coherent_r2 = sum(
        left.coefficient * right.coefficient *
        radial_overlap(left.wave, right.wave, r -> r^2) for
        left in components for right in components if
        left.basis.L_label == right.basis.L_label &&
        left.basis.multiplicity == right.basis.multiplicity &&
        left.basis.J == right.basis.J
    )
    @test radial_expect(spec, mixed[1], r -> r^2) ≈ coherent_r2 atol = 1e-12
    @test_throws ArgumentError radial_wave(spec, mixed[1])
    # A consumer that is specific to one spectroscopic channel receives the
    # coherent radial-n projection of the physical state, not one arbitrarily
    # selected precursor wave.  The remaining D component does not contribute
    # to an S-wave annihilation kernel.
    annihilation_input = GIModel.annihilation_basis_input(
        spec, BasisState(1, "S", 3, 1),
    )
    @test length(annihilation_input.radial_components) == 2
    @test all(isfinite(first(component)) for component in annihilation_input.radial_components)
    @test isfinite(GIModel._sL_smearing_factor(
        MomentumIntegralSmearing(120), annihilation_input, 0,
    ))
    # lookup by quantum numbers
    s = spectrum_state(spec, 1, "P", 3, 2)
    @test s.label == "1^3P_2"
    @test_throws ArgumentError spectrum_state(spec, 3, "S", 1, 0)
    # level beyond the per-channel budget fails loudly
    @test_throws ArgumentError compute_spectrum(
        params, meson;
        levels = [BasisState(7, "S", 1, 0)],
        solver = FiniteDifferenceSolver(ngrid = 80, rmax = 8.0),
    )
end

@testset "compute_spectrum same-J mixing gated by flavor content" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    pair = [BasisState(1, "P", 1, 1), BasisState(1, "P", 3, 1)]
    us = compute_spectrum(
        params, Meson(mq, :u, :s);
        levels = pair,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
    )
    so_mixed = [s for s in us.states if any(m.mechanism == "antisymmetric_spin_orbit" for m in s.mixings)]
    @test length(so_mixed) == 2
    @test all(s.fine_structure_mass_convention == "unequal_mass_same_j_mixed" for s in so_mixed)
    @test sum(s.mass_GeV for s in so_mixed) ≈
          sum(s.mixings[end].unmixed_GeV for s in so_mixed) atol = 1e-10
    @test any(abs(s.mixings[end].offdiag_GeV) > 0 for s in so_mixed)
    # equal flavor: the antisymmetric matrix element vanishes, no block forms
    cc = compute_spectrum(
        params, Meson(mq, :c, :c);
        levels = pair,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
    )
    @test all(isempty(s.mixings) for s in cc.states)
end

@testset "annihilation blocks built from two spectra (no reference data)" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    ps_levels = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
    nn = compute_spectrum(
        params, Meson(mq, :q, :q);
        levels = ps_levels,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )
    ss = compute_spectrum(
        params, Meson(mq, :s, :s);
        levels = ps_levels,
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )

    # the calibrated rank-one block reproduces its targets by construction
    targets = (0.520, 0.960, 1.440, 1.630)
    cal = pseudoscalar_annihilation_block(CalibratedP1Annihilation(), params, nn, ss;
        targets = targets)
    @test sort(cal.masses) ≈ sort(collect(targets)) atol = 1e-10
    calibrated_diagonal = [
        spectrum_state(nn, ps_levels[1]).mass_GeV,
        spectrum_state(ss, ps_levels[1]).mass_GeV,
        spectrum_state(nn, ps_levels[2]).mass_GeV,
        spectrum_state(ss, ps_levels[2]).mass_GeV,
    ]
    calibrated_update = cal.block.matrix - Diagonal(calibrated_diagonal)
    @test all(diag(calibrated_update) .>= 0.0)
    # explicit targets are mandatory - the digitized values live in GIPaper
    @test_throws ArgumentError pseudoscalar_annihilation_block(
        CalibratedP1Annihilation(), params, nn, ss)

    p1 = pseudoscalar_annihilation_block(PaperP1Annihilation(), params, nn, ss)
    p2 = pseudoscalar_annihilation_block(PaperP2Annihilation(), params, nn, ss)
    @test cal isa MixingResult
    @test p1 isa MixingResult
    @test p2 isa MixingResult
    @test isempty(p1.pole_matrices)
    @test length(p2.pole_matrices) == length(p2.masses)
    @test all(isfinite, p1.masses)
    @test all(isfinite, p2.masses)
    @test p1.masses != p2.masses
    pseudoscalar_basis = [
        annihilation_basis_input(nn, ps_levels[1]),
        annihilation_basis_input(ss, ps_levels[1]),
        annihilation_basis_input(nn, ps_levels[2]),
        annihilation_basis_input(ss, ps_levels[2]),
    ]
    @test_throws ArgumentError GIModel.isoscalar_pseudoscalar_annihilation_solution(
        PaperP2Annihilation(), params, pseudoscalar_basis; maxiter = 0,
    )
    @test_throws ArgumentError GIModel.isoscalar_pseudoscalar_annihilation_solution(
        PaperP2Annihilation(), params, pseudoscalar_basis; tol = 0.0,
    )
    # light nn̄ inputs carry the sqrt(2) flavor-coherence label
    input = annihilation_basis_input(nn, ps_levels[1])
    @test input.label == "1 ns"
    @test annihilation_basis_input(ss, ps_levels[2]).label == "2 ss"

    # general Eq. (16) block for one channel across the two flavors
    s1 = BasisState(1, "S", 3, 1)
    nn3 = compute_spectrum(
        params, Meson(mq, :q, :q);
        levels = [s1],
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )
    ss3 = compute_spectrum(
        params, Meson(mq, :s, :s);
        levels = [s1],
        solver = FiniteDifferenceSolver(ngrid = 120, rmax = 12.0),
        terms = SpinTerms(fine_structure = false),
    )
    block = isoscalar_annihilation_block(params, nn3, ss3, s1;
        amplitude_A = params.annihilation.s1_A)
    @test length(block.masses) == 2
    @test all(isfinite, block.masses)
    # trace conservation: eigenvalue sum equals diagonal sum plus block trace
    diag_sum = spectrum_state(nn3, s1).mass_GeV + spectrum_state(ss3, s1).mass_GeV
    annihilation_matrix = block.block.matrix - Diagonal([
        spectrum_state(nn3, s1).mass_GeV,
        spectrum_state(ss3, s1).mass_GeV,
    ])
    @test sum(block.masses) ≈ diag_sum + sum(diag(annihilation_matrix)) atol = 1e-10
end

@testset "PA-16 final isoscalar spectrum owns flavor annihilation" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    solver = FiniteDifferenceSolver(ngrid = 80, rmax = 8.0, nlevels_per_channel = 2)
    levels = [
        BasisState(1, "S", 1, 0),
        BasisState(2, "S", 1, 0),
        BasisState(1, "S", 3, 1),
        BasisState(2, "S", 3, 1),
        BasisState(1, "D", 3, 1),
    ]
    nn = compute_spectrum(params, Meson(mq, :q, :q); levels = levels, solver = solver)
    ss = compute_spectrum(params, Meson(mq, :s, :s); levels = levels, solver = solver)
    final = add_isoscalar_annihilation(
        params, nn, ss;
        amplitudes = Dict(("S", 3, 1) => params.annihilation.s1_A),
    )

    @test final isa MixedSpectrum
    @test [(m.flavor1, m.flavor2) for m in final.channels] == [(:q, :q), (:s, :s)]
    @test length(final.states) == length(nn.states) + length(ss.states)
    @test all(!isnothing(state.basis.flavors) for state in final.states)
    @test_throws ArgumentError spectrum_state(final, 1, "S", 3, 1)

    qS = spectrum_state(final, BasisState(1, "S", 3, 1; flavors = (:q, :q)))
    sS = spectrum_state(final, BasisState(1, "S", 3, 1; flavors = (:s, :s)))
    q_ann = last(qS.mixings)
    s_ann = last(sS.mixings)
    @test q_ann.result === s_ann.result
    @test q_ann.mechanism == "general_eq16_annihilation"
    @test qS.mass_GeV == q_ann.result.masses[q_ann.eigenstate]
    @test sS.mass_GeV == s_ann.result.masses[s_ann.eigenstate]
    @test_throws ArgumentError radial_wave(final, qS)

    # Tensor S/D mixing precedes flavor mixing. The final state must flatten
    # both transformations, preserve flavor identity, and remain normalized.
    components = physical_components(final, qS)
    @test Set(component.basis.flavors for component in components) ==
          Set([(:q, :q), (:s, :s)])
    @test Set(component.basis.L_label for component in components) == Set(["S", "D"])
    @test sum(abs2(component.coefficient) for component in components) ≈ 1.0 atol = 1e-12
    @test radial_expect(final, qS, _ -> 1.0) ≈ 1.0 atol = 1e-12

    # The four-state radial/flavor block uses the same final-spectrum and
    # composition contract: unique native identities and normalized weights.
    p1_final = add_isoscalar_annihilation(
        params, nn, ss; pseudoscalar = PaperP1Annihilation(),
    )
    p1_state = first(sort(
        [
            state for state in p1_final.states if
            any(m -> m.mechanism == "paper_p1_pseudoscalar_annihilation", state.mixings)
        ];
        by = state -> state.mass_GeV,
    ))
    p1_components = physical_components(p1_final, p1_state)
    p1_identities = [
        (component.basis.n, component.basis.L_label, component.basis.multiplicity,
         component.basis.J, component.basis.flavors) for component in p1_components
    ]
    @test length(p1_identities) == length(unique(p1_identities))
    @test sum(abs2(component.coefficient) for component in p1_components) ≈ 1.0 atol = 1e-12
    weights = Dict(
        component.basis => float(i) for (i, component) in enumerate(p1_components)
    )
    @test physical_state_amplitude(p1_final, p1_state) do component
        weights[component.basis]
    end ≈ sum(
        component.coefficient * weights[component.basis] for
        component in p1_components
    )
    @test physical_transition_amplitude(p1_final, p1_state, p1_state) do left, right
        left.basis.flavors == right.basis.flavors ?
        radial_overlap(left.wave, right.wave, _ -> 1.0) : 0.0
    end ≈ 1.0 atol = 1e-10

    # No reference ordering or implicit amplitude exists in the model stage.
    combined_only = add_isoscalar_annihilation(params, nn, ss)
    @test all(
        mixing -> !occursin("annihilation", mixing.mechanism),
        Iterators.flatten(state.mixings for state in combined_only.states),
    )
    mismatched_ss = compute_spectrum(
        params, Meson(mq, :s, :s);
        levels = levels,
        solver = FiniteDifferenceSolver(ngrid = 90, rmax = 8.0, nlevels_per_channel = 2),
    )
    @test_throws ArgumentError add_isoscalar_annihilation(params, nn, mismatched_ss)
    @test_throws ArgumentError add_isoscalar_annihilation(
        params, nn, ss; pseudoscalar_targets = [0.5, 1.0, 1.5, 2.0],
    )
    @test_throws ArgumentError compute_spectrum(
        params, Meson(mq, :q, :q);
        levels = [
            BasisState(1, "S", 3, 1),
            BasisState(1, "S", 3, 1; label = "duplicate"),
        ],
        solver = solver,
    )

    # Architecture boundary: core model code must never learn paper-row types,
    # target masses, or comparison assignment.
    core_bridge = read(joinpath(root, "src", "flavor_mixing.jl"), String)
    @test !occursin("ReferenceState", core_bridge)
    @test !occursin("reference_GeV", core_bridge)
    @test !occursin("compare_reference", core_bridge)

    # The general builder retains the actual requested orbital channel in its
    # shared basis; it must never relabel a P-wave block as an S wave.
    p_level = BasisState(1, "P", 3, 2)
    nnP = compute_spectrum(params, Meson(mq, :q, :q); levels = [p_level], solver = solver)
    ssP = compute_spectrum(params, Meson(mq, :s, :s); levels = [p_level], solver = solver)
    p_block = isoscalar_annihilation_block(
        params, nnP, ssP, p_level; amplitude_A = params.annihilation.a_3p2,
    )
    @test all(basis.L_label == "P" for basis in p_block.block.basis)
    @test all(basis.flavors in ((:q, :q), (:s, :s)) for basis in p_block.block.basis)
end

@testset "Annihilation phase convention, on the spectrum's own waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    lv = [BasisState(1, "S", 1, 0), BasisState(2, "S", 1, 0)]
    solver = FiniteDifferenceSolver(ngrid = 80, rmax = 8.0)
    spec = compute_spectrum(params, Meson(mq, :q, :q); levels = lv, solver = solver)

    # A spectrum holds exactly one wave per level, from its own solver. The
    # second, oscillator-basis cache that used to sit beside it existed to work
    # around a normalization bug (FD Euclidean vs HO physical, a factor
    # 1/sqrt(h)); with that fixed the two solvers agree and the cache is gone.
    @test fieldnames(GIModel.SectorComputation) == (:params, :solver, :channel_cache)

    # The phase convention is NOT part of that fossil and must survive: the
    # eigensolver returns arbitrary column signs, and Table III amplitude signs
    # depend on Phi(0) proportional to the integral of r*u(r) being positive.
    for n in 1:2
        input = GIModel.annihilation_basis_input(spec, BasisState(n, "S", 1, 0))
        coefficient, w = only(input.radial_components)
        @test coefficient == 1.0
        @test sum(w.r .* w.u) * w.h > 0
    end

    # Applied without mutating the cached wave. Channel solutions already have
    # a deterministic outer-lobe phase; the annihilation convention may choose
    # the same or opposite sign for a particular radial excitation.
    key = RadialChannelKey(only(spec.channels).constituent_masses, "S", 1, 0)
    cached = spec.computation.channel_cache[key]
    before = [copy(w.u) for w in cached.waves]
    raw = radial_wave(spec, "1^1S_0")
    fixed = only(GIModel.annihilation_basis_input(
        spec, BasisState(1, "S", 1, 0),
    ).radial_components)[2]
    @test abs.(raw.u) ≈ abs.(fixed.u)
    @test sum(fixed.r .* fixed.u) > 0
    @test [w.u for w in cached.waves] == before


end

@testset "central diagnostic is independent of fixed -> mixed production" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    us = Meson(mq, :u, :s)
    levels = spectrum_levels(2; L_labels = ("S", "P"))
    kwargs = (solver = FiniteDifferenceSolver(ngrid = 250, rmax = 16.0),)

    central = central_spectrum(params, us; levels = levels, kwargs...)
    @test central isa CentralSpectrum
    @test eltype(central.states) === CentralState
    @test parameters(central) === central.computation.params
    @test all(isfinite(s.central_GeV) for s in central.states)

    corrected = fixed_spectrum(params, us; levels = levels, kwargs...)
    @test corrected isa CorrectedSpectrum
    @test corrected.computation !== central.computation
    @test all(k.multiplicity == 0 for k in keys(central.computation.channel_cache))
    @test all(k.multiplicity != 0 for k in keys(corrected.computation.channel_cache))
    for s in corrected.states
        @test s.basis isa BasisState
        @test !hasfield(CorrectedState, :central)
        @test s.mass_GeV ≈ s.central_GeV + s.contact_shift_GeV + s.fine_structure_shift_GeV
    end

    mixed = add_intra_meson_mixing(corrected)
    @test mixed isa MixedSpectrum
    @test mixed.computation === corrected.computation

    # compute_spectrum composes only the production stages
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
        fixed_spectrum(
            params,
            us;
            levels = levels,
            kwargs...,
            terms = SpinTerms(contact_hyperfine = false, fine_structure = false),
        ),
    )
    @test all(s.mass_GeV == s.central_GeV for s in bare.states)
    @test all(isempty(s.mixings) for s in bare.states)
    @test all(s.fine_structure_mass_convention == "disabled" for s in bare.states)

    # stage order is enforced by dispatch
    @test !isdefined(GIModel, :add_spin_corrections)
    @test !hasmethod(add_intra_meson_mixing, Tuple{CentralSpectrum})
end

@testset "nonperturbative contact states expose hyperfine-distinct S waves" begin
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    masses = Meson(mq, :q, :q).constituent_masses
    central = central_spectrum(params, Meson(mq, :q, :q); levels = spectrum_levels(1))
    central_wave = radial_wave(central.computation.channel_cache[RadialChannelKey(masses, "S")], 1)
    r = central_wave.r

    sol1 = contact_hyperfine_nonperturbative_states(params, masses, "S", 1, r, 2)
    sol3 = contact_hyperfine_nonperturbative_states(params, masses, "S", 3, r, 2)
    @test radial_wave(sol1, 1).r == radial_wave(sol3, 1).r == r
    @test !isempty(sol1.waves) && !isempty(sol3.waves)
    # levels agree with the energy-only accessor
    @test sol1.eigenvalues_GeV ≈ GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 1, r, 2)
    @test sol3.eigenvalues_GeV ≈ GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 3, r, 2)

    # ^1S_0 (pi) is more compact than ^3S_1 (rho): smaller <r^2>, lower energy
    @test radial_expect(radial_wave(sol1, 1), x -> x^2) <
          radial_expect(radial_wave(sol3, 1), x -> x^2)
    @test sol1.eigenvalues_GeV[1] < sol3.eigenvalues_GeV[1]

    # The oscillator basis resums the same operator in its own space and lands
    # on the same answer, so this wrapper is basis-generic. "Empty" is reserved
    # for the genuinely inactive cases (non-S wave, sandwich off); it is never
    # "this basis has no implementation", which is what used to make callers
    # substitute first-order PT and report 0.2842 GeV for the light 1S0.
    sol_ho = contact_hyperfine_nonperturbative_states(
        params, masses, "S", 1, 2; solver = OscillatorSolver())
    @test abs(sol_ho.eigenvalues_GeV[1] - sol1.eigenvalues_GeV[1]) < 1e-3
end

@testset "state accessors take labels and BasisStates alike" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    spec = compute_spectrum(params, Meson(mq, :c, :c);
        levels = spectrum_levels(1; L_labels = ("S", "D")))
    level = BasisState(1, "S", 1, 0)
    @test radial_wave(spec, level).u == radial_wave(spec, "1^1S_0").u
    @test [c.coefficient for c in physical_components(spec, BasisState(1, "S", 3, 1))] ==
          [c.coefficient for c in physical_components(spec, "1^3S_1")]
    @test physical_state_amplitude(c -> 1.0, spec, level) ==
          physical_state_amplitude(c -> 1.0, spec, "1^1S_0")
end

@testset "helpful errors and compact displays" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    err = try
        compute_isoscalar_spectrum(params, Meson(mq, :q, :q), Meson(mq, :s, :s);
            levels = [BasisState(1, "S", 1, 0)], pseudoscalar = PaperP1Annihilation())
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("pseudoscalar_basis", err.msg)

    spec = compute_spectrum(params, Meson(mq, :c, :s);
        levels = spectrum_levels(1; L_labels = ("P",)), solver = OscillatorSolver())
    state = spectrum_state(spec, "1^1P_1")
    shown = sprint(show, MIME"text/plain"(), state)
    @test occursin("antisymmetric_spin_orbit", shown) && length(shown) < 1000
    @test length(sprint(show, MIME"text/plain"(), fixed_spectrum(params, Meson(mq, :c, :s);
        levels = spectrum_levels(1; L_labels = ("P",))))) < 1000
    @test length(sprint(show, first(values(spec.computation.channel_cache)))) < 200
    @test length(sprint(show, MIME"text/plain"(), params)) < 1000
end

@testset "Single-channel and isoscalar defaults use the same solver" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    levels = [BasisState(1, "S", 1, 0)]
    nn, ss = Meson(mq, :q, :q), Meson(mq, :s, :s)
    iso = compute_isoscalar_spectrum(params, nn, ss; levels)
    ordinary = compute_spectrum(params, nn; levels)
    @test iso.computation.solver == ordinary.computation.solver == FiniteDifferenceSolver()
    @test spectrum_state(iso, BasisState(1, "S", 1, 0; flavors = (:q, :q))).mass_GeV ==
          only(ordinary.states).mass_GeV
end

@testset "Public convergence records retain source-channel scope" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    levels = [BasisState(1, "S", 1, 0)]
    for solver in (FiniteDifferenceSolver(ngrid=40), OscillatorSolver(nbasis=12, beta_grid=[0.6], converge=false))
        spec = compute_spectrum(params, Meson(mq,:c,:c); levels, solver)
        records = convergence(spec, "1^1S_0")
        @test length(records) == 1
        @test records == convergence(spec, only(spec.states))
        @test records[1].basis == only(spec.states).basis
        @test solver isa OscillatorSolver ? records[1].certificate.status == :unchecked : isnothing(records[1].certificate)
    end
end
