# Map reference CSV (sector, composition_raw) to (m_quark, m_antiquark) in GeV.

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

function mass_one_quark_token(params::GIParameters, tok::String)
    t = _str_strip_quark(String(tok))
    t = replace(t, "bar" => "", "quark" => "")
    t = strip(t)
    t == "u" && return params.masses["u"]
    t == "d" && return params.masses["d"]
    t == "s" && return params.masses["s"]
    t == "c" && return params.masses["c"]
    t == "b" && return params.masses["b"]
    t in ("n",) && return params.masses["q"]
    t in ("q",) && return params.masses["q"]
    error("unknown quark token `$tok` (parsed=`$t`)")
end

function mass_anti_token(params::GIParameters, tok::String, sector::String)
    t = _str_strip_quark(String(tok))
    if !endswith(t, "bar")
        error("antiquark field must end in `bar` in sector $sector, got `$tok`")
    end
    core = String(t[1:prevind(t, end, 3)])
    return mass_one_quark_token(params, core)
end

function parse_quark_masses(
    params::GIParameters,
    sector::String,
    composition_raw::AbstractString,
)
    if sector == "ccbar" || sector == "charmonium"
        m = params.masses["c"]
        return m, m
    end
    if sector == "bbbar" || sector == "bottomonium"
        m = params.masses["b"]
        return m, m
    end
    if sector == "isovector" || sector == "isoscalar"
        m = 0.5 * (params.masses["u"] + params.masses["d"])
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
        # First segment: e.g. "b ubar", "c sbar", "-c dbar", "-u sbar", "b cbar"
        seg = first_content_segment(composition_raw)
        parts = [String(s) for s in eachsplit(seg, isspace) if !isempty(s) && s != "—"]
        length(parts) < 2 && error("need 2+ tokens in first segment of `$composition_raw` (sector $sector)")
        m1 = mass_one_quark_token(params, String(parts[1]))
        m2 = mass_anti_token(params, String(parts[2]), sector)
        return m1, m2
    end
    error("unknown sector `$sector` for mass resolution")
end
