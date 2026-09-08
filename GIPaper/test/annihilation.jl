@testset "isoscalar pseudoscalar annihilation block is opt-in" begin
    reference = load_ref(reference_spectrum_path("isoscalar"))
    sub = reference[1:4]
    common = (contact_hyperfine = true, use_fine_structure = false, ngrid = 120, rmax = 12.0)
    plain = compare_reference(params, mq, sub; common...)
    mixed = compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :calibrated_p1, strange_mass_GeV = mq["s"],
        common...,
    )
    @test maximum(abs(row.residual_MeV) for row in plain) > 0.3e3
    @test maximum(abs(row.residual_MeV) for row in mixed) < 1e-8
    @test all(row.annihilation_scheme == "calibrated_p1" for row in mixed)
    @test any(abs(row.annihilation_shift_GeV) > 0.1 for row in mixed)
    # the scheme requires the strange partner mass
    @test_throws ArgumentError compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :calibrated_p1, common...,
    )

    p1 = compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :paper_p1, strange_mass_GeV = mq["s"], common...,
    )
    p2 = compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :paper_p2, strange_mass_GeV = mq["s"], common...,
    )
    p1_short = compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :p1, strange_mass_GeV = mq["s"], common...,
    )
    p2_short = compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :p2, strange_mass_GeV = mq["s"], common...,
    )
    @test all(row.annihilation_scheme == "paper_p1" for row in p1)
    @test all(row.annihilation_scheme == "paper_p2" for row in p2)
    @test all(isfinite(row.predicted_GeV) for row in p1)
    @test all(isfinite(row.predicted_GeV) for row in p2)
    @test [row.predicted_GeV for row in p1] != [row.predicted_GeV for row in p2]
    @test [row.predicted_GeV for row in p1_short] ≈ [row.predicted_GeV for row in p1]
    @test [row.predicted_GeV for row in p2_short] ≈ [row.predicted_GeV for row in p2]
end

@testset "isoscalar general_s1 annihilation splits omega/phi" begin
    reference = load_ref(reference_spectrum_path("isoscalar"))
    sub = [r for r in reference if r.L == "S" && r.n == 1 && r.multiplicity == 3 && r.J == 1]
    @test length(sub) == 2
    common = (contact_hyperfine = true, use_fine_structure = false, ngrid = 120, rmax = 12.0)
    plain = compare_reference(params, mq, sub; common...)
    mixed = compare_reference(
        params, mq, sub;
        isoscalar_pseudoscalar_annihilation = :general_s1, strange_mass_GeV = mq["s"],
        common...,
    )
    @test all(row.annihilation_scheme == "none" for row in plain)
    @test all(row.annihilation_scheme == "general_s1" for row in mixed)
    # mixing must produce a real splitting
    @test abs(mixed[1].predicted_GeV - mixed[2].predicted_GeV) >
          abs(plain[1].predicted_GeV - plain[2].predicted_GeV)
    # Each row now retains its actual q-qbar or s-sbar precursor rather than
    # pretending both physical states started from the q-qbar diagonal.
    @test sum(row.predicted_GeV for row in mixed) ≈
          sum(row.isoscalar_annihilation_unmixed_GeV for row in mixed) +
          sum(row.annihilation_shift_GeV for row in mixed) atol = 1e-10
    @test Set(row.m1_GeV for row in mixed) == Set([mq["q"], mq["s"]])
end

@testset "table_iii isoscalar scheme: ideal mixing + general Eq.(16)" begin
    reference = load_ref(reference_spectrum_path("isoscalar"))
    sub = [
        r for r in reference if r.n == 1 && (
            (r.L == "S" && r.multiplicity == 3 && r.J == 1) ||  # omega/phi
            (r.L == "P" && r.multiplicity == 3 && r.J == 2) ||  # f2/f2'
            (r.L == "P" && r.multiplicity == 1 && r.J == 1)     # h1/h1' (ideal)
        )
    ]
    @test length(sub) == 6
    rows = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = true,
        isoscalar_pseudoscalar_annihilation = :table_iii, strange_mass_GeV = mq["s"],
        ngrid = 120, rmax = 12.0,
    )
    @test length(rows) == 6
    by_ref(pred) = sort([r for r in rows if pred(r)]; by = r -> r.reference_GeV)

    # omega/phi: general Eq. (16) with A(^3S_1)=+2.5 and the three-gluon bracket.
    s1 = by_ref(r -> r.L == "S")
    @test all(r.isoscalar_annihilation_scheme == "general_eq16" for r in s1)
    omega, phi = s1
    # near-ideal mixing with a small ss admixture in the omega (Table III: -0.02)
    @test abs(omega.isoscalar_component_ns) > 0.99
    @test -0.15 < omega.isoscalar_component_ss < 0.0
    # the annihilation shift on the omega is small and positive (paper: +10 MeV)
    @test 0.0 < omega.annihilation_shift_GeV < 0.05
    # phi sits near the strange-channel diagonal, well above the omega
    @test phi.predicted_GeV - omega.predicted_GeV > 0.15

    # f2/f2': general Eq. (16) with A(^3P_2)=-0.8 and the two-gluon bracket.
    p2 = by_ref(r -> r.L == "P" && r.multiplicity == 3)
    @test all(r.isoscalar_annihilation_scheme == "general_eq16" for r in p2)
    f2, f2p = p2
    # negative amplitude pushes the f2 down and gives it a positive ss component
    @test f2.annihilation_shift_GeV < 0.0
    @test 0.0 < f2.isoscalar_component_ss < 0.2
    @test abs(f2p.isoscalar_component_ss) > 0.97

    # h1/h1': ideally mixed, heavier row gets the strange-channel prediction.
    h1 = by_ref(r -> r.L == "P" && r.multiplicity == 1)
    @test all(r.isoscalar_annihilation_scheme == "ideal" for r in h1)
    @test h1[1].annihilation_shift_GeV == 0.0
    @test h1[1].isoscalar_component_ns == 1.0
    @test h1[2].isoscalar_component_ss == 1.0
    @test 0.15 < h1[2].predicted_GeV - h1[1].predicted_GeV < 0.35
    @test h1[2].m1_GeV ≈ mq["s"]
    # assigned isoscalar rows become scoreable in the non-mixing summary
    @test !GIPaper.mixing_prone_state(h1[2])
    @test !GIPaper.mixing_prone_state(omega)
end
