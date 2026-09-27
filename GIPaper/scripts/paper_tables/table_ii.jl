using TOML
include(joinpath(@__DIR__, "common.jl"))

const TABLE_II_MODEL_KEYS = Dict(
    "m_ud_avg" => ("masses", "m_ud_avg_MeV"), "m_s" => ("masses", "m_s_MeV"),
    "m_c" => ("masses", "m_c_MeV"), "m_b" => ("masses", "m_b_MeV"),
    "b" => ("potential", "b_GeV2"),
    "c" => ("potential", "c_MeV"), "sigma0" => ("relativistic_smearing", "sigma0_GeV"),
    "s" => ("relativistic_smearing", "s"),
    "epsilon_c" => ("relativistic_factors", "epsilon_c"),
    "epsilon_t" => ("relativistic_factors", "epsilon_t"),
    "epsilon_so_vector" => ("relativistic_factors", "epsilon_so_vector"),
    "epsilon_so_scalar" => ("relativistic_factors", "epsilon_so_scalar"),
)

# Lambda only motivates the fixed Gaussian running-coupling profile; it is not a
# runtime input (docs/model_inputs.md), so it has no active value.
function active_table_ii_value(active, parameter)
    parameter == "alpha_s_critical" && return GIModel.alpha_s_q(0.0)
    parameter == "Lambda" && return "—"
    section, key = TABLE_II_MODEL_KEYS[parameter]
    return active[section][key]
end

function compute_table_ii()
    source = joinpath(GIPAPER_DIR, "data", "raw", "digitized_tables",
        "table_ii_parameters", "table_ii_parameters.csv")
    active = TOML.parsefile(GIModel.default_parameters_path())
    lines = ["parameter\tactive_value\tunit\tgi_value"]
    for row in filter(!isempty, readlines(source)[2:end])
        fields = split(row, ','; limit=10)
        parameter = fields[5]
        gi_value = parse(Float64, fields[7])
        unit = fields[8]
        value = active_table_ii_value(active, parameter)
        push!(lines, join((parameter, value, unit, gi_value), '\t'))
    end
    return write_paper_table("table_ii.tsv", join(lines, '\n'))
end

abspath(PROGRAM_FILE) == (@__FILE__) && compute_table_ii()
