# Invoked by `Pkg.test()` with the GIPaper environment active.
# To run manually: `julia --project=GIPaper GIPaper/test/runtests.jl` from the repo root.

using Test
using GIModel
using GIPaper

# Until GIModel drops its legacy comparison exports, some names are exported by
# both packages; qualify the GIPaper ones explicitly.
const load_ref = GIPaper.load_reference_spectrum

params, mq = load_parameters_and_quark_masses(GIPaper.model_parameters_path())

@testset "data CSV checks" begin
    script = joinpath(dirname(@__DIR__), "scripts", "data_checks.py")
    p = run(`python3 $script validate`, wait = false)
    wait(p)
    @test success(p)
end

@testset "reference loading" begin
    ccbar = load_ref(reference_spectrum_path("charmonium"))
    bbbar = load_ref(reference_spectrum_path("bottomonium"))
    @test length(ccbar) == 28
    @test length(bbbar) == 30
    @test ccbar[1].quark_content == "c cbar"
    @test ccbar[1].composition == "1^1S_0"
    @test_throws ArgumentError reference_spectrum_path("nonexistent_sector")
end

@testset "reference_meson resolution (loud, no fallback)" begin
    row(sector, content) = GIPaper.ReferenceState(sector, content, "1^1S_0", 1, 1, "S", 0, 1.0, "high")
    m = reference_meson(mq, row("ccbar", "ignore"))
    @test m.constituent_masses.m1_GeV ≈ m.constituent_masses.m2_GeV ≈ mq["c"]
    m = reference_meson(mq, row("charmonium", "c cbar"))
    @test (m.flavor1, m.flavor2) == (:c, :c)
    # u/d/n/q all resolve to the one LightQuark the model can represent, so the
    # light leg reads `:q` regardless of how the CSV spells it. Masses are
    # unchanged (m_u = m_d = m_q in the parameter set).
    m = reference_meson(mq, row("charmed", "-c dbar; c ubar"))
    @test (m.flavor1, m.flavor2) == (:c, :q)
    @test m.constituent_masses.m1_GeV ≈ mq["c"]
    @test m.constituent_masses.m2_GeV ≈ mq["d"]
    m = reference_meson(mq, row("charmed_strange", "c sbar"))
    @test (m.flavor1, m.flavor2) == (:c, :s)
    m = reference_meson(mq, row("bottom_light", "b ubar; -b dbar"))
    @test (m.flavor1, m.flavor2) == (:b, :q)
    m = reference_meson(mq, row("strange", "-u sbar; -d sbar"))
    @test (m.flavor1, m.flavor2) == (:q, :s)
    m = reference_meson(mq, row("isoscalar", "n nbar / s sbar mixed"))
    @test (m.flavor1, m.flavor2) == (:q, :q)
    @test m.constituent_masses.m1_GeV ≈ 0.5 * (mq["u"] + mq["d"])
    # failures are loud and name the offending row
    @test_throws ArgumentError reference_meson(mq, row("unknown_sector", "c cbar"))
    @test_throws ArgumentError reference_meson(mq, row("charmed", "gibberish"))
    @test_throws ArgumentError reference_meson(mq, row("charmed", "x ybar"))
    @test_throws ArgumentError reference_meson(mq, row("charmed", "c d"))  # antiquark must end in bar
end

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
    # sum of masses is conserved up to the recorded annihilation shifts
    @test sum(row.predicted_GeV for row in mixed) ≈
          sum(row.predicted_GeV for row in plain) + sum(row.annihilation_shift_GeV for row in mixed) atol = 1e-10
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

@testset "write_residual_report renders rows" begin
    reference = load_ref(reference_spectrum_path("charmonium"))
    rows = compare_reference(params, mq, reference[1:4]; ngrid = 100, rmax = 10.0)
    mktempdir() do d
        path = joinpath(d, "report.md")
        GIPaper.write_residual_report(path, "Residuals: test", rows)
        text = read(path, String)
        @test occursin("# Residuals: test", text)
        @test occursin("Contribution Breakdown", text)
        @test occursin("Spin-Averaged Diagnostics", text)
    end
end
