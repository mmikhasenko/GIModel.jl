# Table-II-style quark masses from the parameters TOML (`[masses]` block).
#
# Public API (exported from GIModel.jl):
#   QuarkMassTable, load_quark_masses, load_parameters_and_quark_masses

"""Flavor-keyed constituent masses from TOML (GeV): `u`, `d`, `q`, `s`, `c`, `b`."""
const QuarkMassTable = Dict{String,Float64}

function quark_masses_from_raw(raw)::QuarkMassTable
    m = raw["masses"]
    return Dict{String,Float64}(
        "u" => m["m_ud_avg_MeV"] / 1000,
        "d" => m["m_ud_avg_MeV"] / 1000,
        "q" => m["m_ud_avg_MeV"] / 1000,
        "s" => m["m_s_MeV"] / 1000,
        "c" => m["m_c_MeV"] / 1000,
        "b" => m["m_b_MeV"] / 1000,
    )
end

function load_quark_masses(path::AbstractString)::QuarkMassTable
    return quark_masses_from_raw(TOML.parsefile(path))
end

function load_parameters_and_quark_masses(path::AbstractString)
    raw = TOML.parsefile(path)
    return gi_parameters_from_raw(raw), quark_masses_from_raw(raw)
end
