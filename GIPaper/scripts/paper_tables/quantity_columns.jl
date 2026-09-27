include(joinpath(@__DIR__, "common.jl"))
using TOML

# The canonical Table III/V/VI/VII calculations are GIPaper's residual reports;
# this registry records what each of their numerical columns means.
const REPORT_SOURCES = Dict(
    "table_iii" => ("docs/residual_reports/table_iii_mixing_audit.md", "scripts/audit_table_iii_mixings.jl"),
    "table_v" => ("docs/residual_reports/table_v_reproduction.md", "scripts/reproduce_table_v.jl"),
    "table_vi" => ("docs/residual_reports/table_vi_photon_decays.md", "scripts/audit_table_vi_photon_decays.jl"),
    "table_vii" => ("docs/residual_reports/table_vii_annihilation_em.md", "scripts/audit_table_vii.jl"),
)

"""Column provenance for already calculated report values (no invented GI reference).

Reports retain their printed precision. This registry distinguishes predictions,
derived kinematics, and algebraic inputs, even without a paper comparison column.
Sources and drivers are relative to the GIPaper package root.
"""
function compute_quantity_columns()
    definitions = [
        ("table_v", "q MeV", "derived kinematics", "MeV", "Two-body momentum from the audit's parent and daughter masses; src records the mass convention."),
        ("table_vi", "q MeV", "derived kinematics", "MeV", "Photon momentum from shared PDG 2026 experimental masses."),
        ("table_vii", "M (GeV)", "experimental mass input", "GeV", "Shared PDG 2026 experimental mass."),
        ("table_vii", "α_s(M)", "derived coupling", "", "GIModel.alpha_s_q(M), the GI Gaussian coupling evaluated at Q=M (paper notation alpha_s(M^2))."),
        ("table_vii", "S_L", "wavefunction integral", "GeV^(3/2)", "Eq. (17) wavefunction_origin_smearing; S_L includes the paper's momentum smearing."),
        ("table_vii", "M̃ (GeV)", "wavefunction integral", "GeV", "mock_meson_mass: expectation of E1+E2 over the same solved wave."),
        ("table_vii", "value", "wavefunction integral", "", "leptonic_decay_factor: the indicated P_P, V_V, V'_V or P'_A1 factor before its charge/flavor coefficient."),
        ("table_vii", "q_eff", "algebraic input", "", "Effective squared charge sum a_i e_i^2 in the stated pure/ideal flavor convention; not a fitted prediction."),
        ("table_vii", "M_P (GeV)", "physical mass input", "GeV", "Shared PDG 2026 experimental mass for the two-photon kinematics."),
        ("table_vii", "M model (GeV)", "model prediction", "GeV", "Computed eigenvalue of the final four-component P1-mixed state; the amplitude uses the separate physical M_P input."),
        [("table_vii", b, "mixing amplitude", "", "Signed native P1 component; its square is a probability, not the amplitude itself.") for b in ("1nn", "1ss", "2nn", "2ss")]...,
        ("table_vii", "f_model", "model prediction", "", "Absolute model vector decay constant used in Eq. (D8)."),
        ("table_vii", "width (D8)", "derived width", "GeV", "dilepton_vector_width(f,M) = (4pi/3) alpha_EM^2 M f^2; uses the physical mass, without a QCD radiative factor."),
        ("table_iii", "mass GeV", "model prediction", "GeV", "Eigenvalue of the stated annihilation block."),
    ]
    entries = [Dict("object"=>o, "column"=>c, "role"=>r, "unit"=>u, "method"=>m,
        "source"=>REPORT_SOURCES[o][1], "driver"=>REPORT_SOURCES[o][2])
        for (o, c, r, u, m) in definitions]
    path = joinpath(PAPER_TABLES_DIR, "quantity_columns.toml")
    mkpath(PAPER_TABLES_DIR)
    open(path, "w") do io
        TOML.print(io, Dict("column"=>entries); sorted=true)
    end
    return path
end
abspath(PROGRAM_FILE) == (@__FILE__) && compute_quantity_columns()
