# Map a reference CSV row (sector string + quark-content string) to a
# GIModel.Meson. Replaces the old masses_from_content.jl: parsing failures
# throw ArgumentError naming the row — there is no silent fallback mass.
#
# Public API (exported from GIPaper.jl): reference_meson, quark_for

"""
    quark_for(quark_masses, flavor::Symbol) -> AbstractQuark

Map a paper flavor label to a GIModel quark object. This is the comparison
layer's job: the paper names flavors, the engine only needs mass and charge.

`:u`, `:d`, `:n` and `:q` all resolve to the same [`LightQuark`](@ref) — the GI
parameter set gives u and d one mass and nothing in the dynamics branches on
which, so that distinction is CSV spelling rather than physics. (The digitized
tables are inconsistent about it anyway: `_SECTOR_FLAVORS` labels isovector and
isoscalar `:q` while the open-flavor rows spell the same light quark `d` or `u`.)
"""
function quark_for(quark_masses::QuarkMassTable, flavor::Symbol)
    flavor in (:u, :d, :n, :q) && return LightQuark(quark_masses["q"])
    flavor === :s && return StrangeQuark(quark_masses["s"])
    flavor === :c && return HeavyQuark{:up}(quark_masses["c"], :c)
    flavor === :b && return HeavyQuark{:down}(quark_masses["b"], :b)
    throw(ArgumentError(
        "no quark for flavor `$flavor` (expected :u/:d/:n/:q, :s, :c or :b)",
    ))
end

# Sectors whose flavor content is fixed by the sector name alone.
const _SECTOR_FLAVORS = Dict{String,Tuple{Symbol,Symbol}}(
    "ccbar" => (:c, :c),
    "charmonium" => (:c, :c),
    "bbbar" => (:b, :b),
    "bottomonium" => (:b, :b),
    "isovector" => (:q, :q),
    "isoscalar" => (:q, :q),
)

# Open-flavor sectors parse the first `;`-segment of the quark-content string.
const _OPEN_FLAVOR_SECTORS = (
    "charmed",
    "charmed_strange",
    "strange",
    "bottom_light",
    "bottom_strange",
    "bottom_charm",
    "b_flavored",
)

function _quark_token_flavor(token::AbstractString, state::ReferenceState)
    t = lowercase(strip(String(token)))
    t = lstrip(t, Char['-', ' ', '\t'])
    t = strip(replace(t, "bar" => "", "quark" => ""))
    flavor = get(
        Dict("u" => :u, "d" => :d, "s" => :s, "c" => :c, "b" => :b, "n" => :q, "q" => :q),
        t,
        nothing,
    )
    isnothing(flavor) && throw(ArgumentError(
        "unknown quark token `$token` in row (sector=$(state.sector), quark_content=`$(state.quark_content)`)",
    ))
    return flavor
end

function _open_flavor_pair(state::ReferenceState)
    segment = strip(first(split(state.quark_content, ";", keepempty = false)))
    parts = [String(s) for s in eachsplit(segment, isspace) if !isempty(s) && s != "—"]
    length(parts) >= 2 || throw(ArgumentError(
        "need 2+ tokens in first segment of `$(state.quark_content)` (sector $(state.sector))",
    ))
    endswith(lowercase(strip(lstrip(parts[2], Char['-', ' ']))), "bar") || throw(ArgumentError(
        "antiquark field must end in `bar` in sector $(state.sector), got `$(parts[2])`",
    ))
    return _quark_token_flavor(parts[1], state), _quark_token_flavor(parts[2], state)
end

"""
    reference_meson(quark_masses, state::ReferenceState) -> Meson

Resolve the flavor content of a reference row: closed-flavor sectors
(`charmonium`, `bottomonium`, `isovector`, `isoscalar`, …) map by sector name;
open-flavor sectors parse the first `;`-segment of the quark-content string
(e.g. `"c dbar; c ubar"` → `Meson(mq, :c, :d)`).

Unknown sectors, malformed content, or unknown quark tokens throw
`ArgumentError` naming the row. If a CSV row trips this, fix the CSV — there is
deliberately no fallback mass.
"""
function reference_meson(quark_masses::QuarkMassTable, state::ReferenceState)
    flavors = get(_SECTOR_FLAVORS, state.sector, nothing)
    if !isnothing(flavors)
        return Meson((quark_for(quark_masses, f) for f in flavors)...)
    end
    state.sector in _OPEN_FLAVOR_SECTORS || throw(ArgumentError(
        "unknown sector `$(state.sector)` for flavor resolution (quark_content=`$(state.quark_content)`)",
    ))
    f1, f2 = _open_flavor_pair(state)
    return Meson(quark_for(quark_masses, f1), quark_for(quark_masses, f2))
end
