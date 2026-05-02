# Map reference CSV (sector, composition_raw) to (m_quark, m_antiquark) in GeV.
#
# Public API (exported from GIModel.jl):
#   parse_quark_masses, resolve_constituent_masses, attach_constituent_masses

function first_content_segment(composition_raw::AbstractString)
    sc = String(composition_raw)
    p = split(sc, ";", keepempty = false)
    length(p) >= 1 || return strip(sc)
    return String(strip(p[1]))
end

function _str_strip_quark(s::String)
    t = String(lowercase(strip(s)))
    t = lstrip(t, Char['-', ' ', '\t'])
    return t
end

function mass_one_quark_token(masses::QuarkMassTable, tok::String)
    t = _str_strip_quark(String(tok))
    t = replace(t, "bar" => "", "quark" => "")
    t = strip(t)
    t == "u" && return masses["u"]
    t == "d" && return masses["d"]
    t == "s" && return masses["s"]
    t == "c" && return masses["c"]
    t == "b" && return masses["b"]
    t in ("n",) && return masses["q"]
    t in ("q",) && return masses["q"]
    error("unknown quark token `$tok` (parsed=`$t`)")
end

function mass_anti_token(masses::QuarkMassTable, tok::String, sector::String)
    t = _str_strip_quark(String(tok))
    if !endswith(t, "bar")
        error("antiquark field must end in `bar` in sector $sector, got `$tok`")
    end
    core = String(t[1:prevind(t, end, 3)])
    return mass_one_quark_token(masses, core)
end

function parse_quark_masses(
    masses::QuarkMassTable,
    sector::String,
    composition_raw::AbstractString,
)
    if sector == "ccbar" || sector == "charmonium"
        m = masses["c"]
        return m, m
    end
    if sector == "bbbar" || sector == "bottomonium"
        m = masses["b"]
        return m, m
    end
    if sector == "isovector" || sector == "isoscalar"
        m = 0.5 * (masses["u"] + masses["d"])
        return m, m
    end
    if sector in (
        "charmed",
        "charmed_strange",
        "strange",
        "bottom_light",
        "bottom_strange",
        "bottom_charm",
    )
        seg = first_content_segment(composition_raw)
        parts = [String(s) for s in eachsplit(seg, isspace) if !isempty(s) && s != "—"]
        length(parts) < 2 &&
            error("need 2+ tokens in first segment of `$composition_raw` (sector $sector)")
        m1 = mass_one_quark_token(masses, String(parts[1]))
        m2 = mass_anti_token(masses, String(parts[2]), sector)
        return m1, m2
    end
    error("unknown sector `$sector` for mass resolution")
end

"""
    resolve_constituent_masses(quark_masses, state::ReferenceState, fallback_mass)

Return [`ConstituentMasses`](@ref) using [`parse_quark_masses`](@ref) on the row's sector and
quark content. On failure (unknown sector, malformed content), return equal masses
`(fallback_mass, fallback_mass)`.

See also [`attach_constituent_masses`](@ref).
"""
function resolve_constituent_masses(
    quark_masses::QuarkMassTable,
    s::ReferenceState,
    fallback_mass::Real,
)
    try
        m1, m2 = parse_quark_masses(quark_masses, String(s.sector), String(s.quark_content))
        return ConstituentMasses(m1, m2)
    catch
        return ConstituentMasses(fallback_mass, fallback_mass)
    end
end

"""
    attach_constituent_masses(quark_masses, states::Vector{ReferenceState}, fallback_mass)

Resolve [`ConstituentMasses`](@ref) **once** per row for [`compute_sector`](@ref) /
[`compare`](@ref). Use `fallback_mass` when parsing fails for a row.

`quark_masses` comes from [`load_quark_masses`](@ref) / [`load_parameters_and_quark_masses`](@ref).
"""
function attach_constituent_masses(
    quark_masses::QuarkMassTable,
    states::Vector{ReferenceState},
    fallback_mass::Real,
)
    mf = float(fallback_mass)
    return ReferenceStateWithMasses[
        ReferenceStateWithMasses(s, resolve_constituent_masses(quark_masses, s, mf)) for s in states
    ]
end
