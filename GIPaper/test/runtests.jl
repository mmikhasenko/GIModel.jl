using Test, CSV, TOML, GIModel, GIPaper
using GIModel.QuarkModelTransitions
const QMT = GIModel.QuarkModelTransitions

@testset "GIPaper reference comparisons" begin
    params, mq = load_parameters_and_quark_masses(model_parameters_path())

    @testset "Paper data" begin
        sectors = ("isovector", "strange", "isoscalar", "charmonium",
                   "charmed", "bottomonium", "b_flavored")
        @test all(!isempty(load_reference_spectrum(reference_spectrum_path(s))) for s in sectors)
        photons = load_table_vi()
        @test length(photons) == 79
        @test Set(keys(load_table_vi_states())) == Set(id for r in photons for id in r.id[2:3])
        @test Set(keys(load_table_policy())) ⊇ Set(["table_v", "table_vi", "table_vii"])
    end

    @testset "Paper solver policy" begin
        select(solver=nothing; ngrid=nothing) = GIPaper._comparison_solver(
            solver; ngrid, rmax=nothing, kinetic=nothing, eigensolver=nothing)
        @test select() isa OscillatorSolver
        fd = FiniteDifferenceSolver(ngrid=80)
        @test select(fd) === fd
        @test select(; ngrid=80) isa FiniteDifferenceSolver
        @test_throws ArgumentError select(OscillatorSolver(); ngrid=80)
    end

    @testset "Charmonium spectrum versus Fig. 6" begin
        reference = filter(r -> r.n == 1 && r.L == "S",
            load_reference_spectrum(reference_spectrum_path("charmonium")))
        rows = compare_reference(params, mq, reference;
            contact_hyperfine = true, use_fine_structure = true)
        @test length(rows) == 2
        # A coarse smoke-test envelope, not a claim of paper-level precision.
        # Full per-state residuals are produced by checks/run_all_spectrum_checks.jl.
        @test maximum(abs(r.residual_MeV) for r in rows) < 100
    end

    @testset "Table V transition rows" begin
        rows = CSV.File(joinpath(paper_data_dir(), "table_v_strong_decays.csv"))
        fraction = mq["c"] / (mq["c"] + mq["d"])
        channels = load_table_v(rows; heavy_fraction_for = _ -> fraction)
        @test !isempty(channels)
        momentum(parent, a, b) = QMT.decay_momentum(
            experimental_mass("V", parent), experimental_mass("V", a), experimental_mass("V", b))
        policy = load_table_policy()["table_v"]
        model = QMT.calibrate_strong_decay_model(momentum("rho", "pi", "pi"),
            momentum("B", "omega", "pi"); convention = QMT.LeadingS0(),
            rho_amplitude = policy["rho_amplitude"], B_amplitude = policy["B_amplitude"])
        # Independent, non-calibration row; 5% allows the pinned modern kinematics.
        reference = only(r for r in rows if r.decay == "phi -> K Kbar" && r.section == "1^3S_1")
        channel = only(load_table_v((reference,)))
        amplitude = QMT.decay_amplitude(model, channel, momentum("phi", "K", "K");
            convention = QMT.LeadingS0()).total
        @test amplitude ≈ parse(Float64, reference.amp_MeV) rtol = 0.05
    end
end
