@testset "compare_reference returns shift breakdown fields" begin
    reference = load_ref(reference_spectrum_path("charmonium"))
    rows = compare_reference(
        params, mq, reference[1:1];
        contact_hyperfine = true, use_fine_structure = true,
        ngrid = 120, rmax = 12.0,
    )
    @test length(rows) == 1
    row = rows[1]
    for field in (
        :central_GeV, :contact_shift_GeV, :spin_orbit_shift_GeV, :tensor_shift_GeV,
        :fine_structure_shift_GeV, :annihilation_shift_GeV, :annihilation_scheme,
        :m1_GeV, :m2_GeV, :fine_structure_mass_convention,
    )
        @test hasproperty(row, field)
    end
    @test isfinite(row.m1_GeV) && isfinite(row.m2_GeV)
    @test row.fine_structure_mass_convention in
          ("equal_mass", "unequal_mass_equal_share_LdotS", "disabled")
    @test row.fine_structure_shift_GeV ≈ row.spin_orbit_shift_GeV + row.tensor_shift_GeV atol = 1e-12
    @test row.predicted_GeV ≈
          row.central_GeV + row.contact_shift_GeV + row.fine_structure_shift_GeV +
          row.annihilation_shift_GeV atol = 1e-12
end

@testset "compare_reference applies same-J antisymmetric spin-orbit mixing" begin
    reference = load_ref(reference_spectrum_path("strange"))
    sub = [
        row for row in reference if
        row.n == 1 && row.L == "P" && row.J == 1 && row.multiplicity in (1, 3)
    ]
    @test length(sub) == 2
    plain = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = true,
        antisymmetric_spin_orbit_mixing = false, ngrid = 120, rmax = 12.0,
    )
    mixed = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = true,
        antisymmetric_spin_orbit_mixing = true, ngrid = 120, rmax = 12.0,
    )
    @test all(row.same_j_mixing_scheme == "none" for row in plain)
    @test all(row.same_j_mixing_scheme == "antisymmetric_spin_orbit" for row in mixed)
    @test all(row.fine_structure_mass_convention == "unequal_mass_same_j_mixed" for row in mixed)
    @test any(abs(row.same_j_offdiag_GeV) > 0 for row in mixed)
    @test sum(row.predicted_GeV for row in mixed) ≈ sum(row.predicted_GeV for row in plain) rtol = 1e-12
    @test all(!GIPaper.mixing_prone_state(row) for row in mixed)
    @test any(abs(mixed[i].predicted_GeV - plain[i].predicted_GeV) > 1e-6 for i in eachindex(mixed))
    # :model_order assigns the same eigenvalue set, possibly to different rows
    model_order = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = true,
        mixed_assignment = :model_order, ngrid = 120, rmax = 12.0,
    )
    @test sort([r.predicted_GeV for r in model_order]) ≈ sort([r.predicted_GeV for r in mixed]) atol = 1e-12
    @test_throws ArgumentError compare_reference(
        params, mq, sub; mixed_assignment = :bogus, ngrid = 120, rmax = 12.0,
    )

    # With more than one radial excitation the production eigensystem is one
    # shared 4x4 block. Comparison rows must project their exact basis rows,
    # never reinterpret the first two rows/columns as an independent 2x2 block.
    multi_sub = [
        row for row in reference if
        row.n in (1, 2) && row.L == "P" && row.J == 1 && row.multiplicity in (1, 3)
    ]
    @test length(multi_sub) == 4
    multi_plain = compare_reference(
        params, mq, multi_sub;
        contact_hyperfine = true, use_fine_structure = true,
        antisymmetric_spin_orbit_mixing = false, ngrid = 120, rmax = 12.0,
    )
    multi_mixed = compare_reference(
        params, mq, multi_sub;
        contact_hyperfine = true, use_fine_structure = true,
        antisymmetric_spin_orbit_mixing = true, ngrid = 120, rmax = 12.0,
    )
    @test all(isfinite(row.same_j_mixing_angle_deg) for row in multi_mixed)
    @test all(isfinite(row.same_j_component_singlet) for row in multi_mixed)
    @test all(isfinite(row.same_j_component_triplet) for row in multi_mixed)
    @test sum(row.predicted_GeV for row in multi_mixed) ≈
          sum(row.predicted_GeV for row in multi_plain) rtol = 1e-12
end

@testset "compare_reference applies tensor same-J triplet L/L' mixing" begin
    reference = load_ref(reference_spectrum_path("isovector"))
    sub = [
        row for row in reference if
        (row.n == 2 && row.L == "S" && row.J == 1 && row.multiplicity == 3) ||
        (row.n == 1 && row.L == "D" && row.J == 1 && row.multiplicity == 3)
    ]
    @test length(sub) == 2
    plain = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = true,
        tensor_mixing = false, ngrid = 120, rmax = 12.0,
    )
    mixed = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = true,
        tensor_mixing = true, ngrid = 120, rmax = 12.0,
    )
    @test all(row.tensor_mixing_scheme == "none" for row in plain)
    @test all(row.tensor_mixing_scheme == "tensor_mixing" for row in mixed)
    @test any(abs(row.tensor_offdiag_GeV) > 0 for row in mixed)
    @test sum(row.predicted_GeV for row in mixed) ≈ sum(row.predicted_GeV for row in plain) rtol = 1e-12
    @test all(!GIPaper.mixing_prone_state(row) for row in mixed)
end

@testset "nonmixing deviation summary excludes mixing-prone rows" begin
    rows = [
        (sector = "charmed", L = "S", n = 1, J = 0, multiplicity = 1, residual_MeV = -30.0),
        (sector = "charmed", L = "P", n = 1, J = 1, multiplicity = 1, residual_MeV = 500.0),
        (sector = "charmed", L = "P", n = 1, J = 1, multiplicity = 3, residual_MeV = -400.0),
        (sector = "charmonium", L = "D", n = 1, J = 1, multiplicity = 3, residual_MeV = 300.0),
        (sector = "charmonium", L = "P", n = 1, J = 2, multiplicity = 3, residual_MeV = 10.0),
        (sector = "isoscalar", L = "S", n = 1, J = 0, multiplicity = 1, residual_MeV = 900.0),
    ]
    @test GIPaper.mixing_prone_state(rows[2])
    @test GIPaper.mixing_prone_state(rows[3])
    @test GIPaper.mixing_prone_state(rows[4])
    @test GIPaper.mixing_prone_state(rows[6])
    summary = GIPaper.nonmixing_deviation_summary(rows)
    by_sector = Dict(row.sector => row for row in summary)
    @test by_sector["charmed"].n_included == 1
    @test by_sector["charmed"].n_excluded_mixing_prone == 2
    @test by_sector["charmed"].mean_abs_deviation_MeV ≈ 30.0
    @test by_sector["charmed"].max_abs_deviation_MeV ≈ 30.0
    @test by_sector["charmonium"].n_included == 1
    @test by_sector["charmonium"].n_excluded_mixing_prone == 1
    @test by_sector["isoscalar"].n_included == 0
    @test isnan(by_sector["isoscalar"].mean_abs_deviation_MeV)
end

@testset "compare_reference is deterministic and subset-stable" begin
    reference = load_ref(reference_spectrum_path("charmonium"))
    sub = reference[1:min(4, length(reference))]
    r1 = compare_reference(params, mq, sub; ngrid = 100, rmax = 10.0)
    r2 = compare_reference(params, mq, sub; ngrid = 100, rmax = 10.0)
    @test length(r1) == length(r2) == length(sub)
    for i in eachindex(r1)
        @test r1[i].predicted_GeV ≈ r2[i].predicted_GeV rtol = 0.0 atol = 1e-14
        @test r1[i].central_GeV ≈ r2[i].central_GeV rtol = 0.0 atol = 1e-14
    end
    # a subset of rows reproduces the same unmixed predictions
    short = sub[1:min(2, length(sub))]
    direct = compare_reference(params, mq, short; ngrid = 100, rmax = 10.0)
    for (i, row) in enumerate(direct)
        j = findfirst(
            r -> r.n == row.n && r.L == row.L && r.multiplicity == row.multiplicity && r.J == row.J,
            r1,
        )
        @test row.central_GeV ≈ r1[j].central_GeV rtol = 1e-12
    end
    # empty input, empty output
    @test isempty(compare_reference(params, mq, GIPaper.ReferenceState[]))
end

@testset "compare_reference(contact_hyperfine=false) drops contact shift" begin
    reference = load_ref(reference_spectrum_path("charmonium"))
    sub = reference[1:1]
    with_hf = compare_reference(
        params, mq, sub;
        contact_hyperfine = true, use_fine_structure = false, ngrid = 120, rmax = 12.0,
    )
    no_hf = compare_reference(
        params, mq, sub;
        contact_hyperfine = false, use_fine_structure = false, ngrid = 120, rmax = 12.0,
    )
    @test length(with_hf) == 1 && length(no_hf) == 1
    @test no_hf[1].contact_shift_GeV ≈ 0.0 atol = 1e-15
    # Full diagonalization lets the radial wave relax when contact is enabled,
    # so the central expectation changes too. The contributions in one solved
    # state, rather than a difference of two eigenproblems, sum to its mass.
    @test with_hf[1].central_GeV != no_hf[1].central_GeV
    @test with_hf[1].predicted_GeV ≈
          with_hf[1].central_GeV + with_hf[1].contact_shift_GeV atol = 1e-12
end
