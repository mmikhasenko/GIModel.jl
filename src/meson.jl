# Meson specification by quark flavors, resolved against the TOML quark-mass table.
#
# Public API (exported from GIModel.jl): Meson, is_equal_flavor

"""Canonical flavor symbols accepted by [`Meson`](@ref), mapped to [`QuarkMassTable`](@ref) keys."""
const MESON_FLAVOR_KEYS = Dict{Symbol,String}(
    :u => "u",
    :d => "d",
    :s => "s",
    :c => "c",
    :b => "b",
    :q => "q",
)

# `:n` is a common alias for the light average `q` in the paper's `n nbar` notation.
_canonical_flavor(flavor::Symbol) = flavor == :n ? :q : flavor

function _flavor_mass_GeV(quark_masses::QuarkMassTable, flavor::Symbol)
    key = get(MESON_FLAVOR_KEYS, flavor, nothing)
    isnothing(key) && throw(ArgumentError(
        "unknown quark flavor `$flavor`; expected one of $(sort(collect(keys(MESON_FLAVOR_KEYS)))) (or the alias :n for :q)",
    ))
    haskey(quark_masses, key) || throw(ArgumentError(
        "quark mass table has no entry for flavor `$flavor` (key `$key`)",
    ))
    return quark_masses[key]
end

"""
    Meson(quark_masses, flavor1, flavor2)
    Meson(flavor1, flavor2, constituent_masses)

A ``q_1 \\bar q_2`` meson channel specified by quark flavors. Flavors are
`:u`, `:d`, `:s`, `:c`, `:b`, or `:q` (the light `u`/`d` average; `:n` is an
accepted alias). The primary constructor resolves constituent masses from a
[`QuarkMassTable`](@ref) and throws `ArgumentError` for unknown flavors — there
is deliberately no silent fallback mass.

The second form takes explicit [`ConstituentMasses`](@ref) for parameter scans.
The antiquark flavor is stored unadorned (`Meson(mq, :c, :u)` is `c ubar`).
"""
struct Meson
    flavor1::Symbol
    flavor2::Symbol
    constituent_masses::ConstituentMasses
    function Meson(flavor1::Symbol, flavor2::Symbol, masses::ConstituentMasses)
        new(_canonical_flavor(flavor1), _canonical_flavor(flavor2), masses)
    end
end

function Meson(quark_masses::QuarkMassTable, flavor1::Symbol, flavor2::Symbol)
    f1 = _canonical_flavor(flavor1)
    f2 = _canonical_flavor(flavor2)
    masses = ConstituentMasses(
        _flavor_mass_GeV(quark_masses, f1),
        _flavor_mass_GeV(quark_masses, f2),
    )
    return Meson(f1, f2, masses)
end

"""
    is_equal_flavor(meson) -> Bool

`true` for self-conjugate flavor content (`c cbar`, `q qbar`, …). Gates the
same-`J` antisymmetric spin-orbit mixing, which vanishes for equal masses.
"""
is_equal_flavor(m::Meson) = m.flavor1 == m.flavor2

reduced_mass(m::Meson) = reduced_mass(m.constituent_masses)

flavor_label(m::Meson) = string(m.flavor1, m.flavor2)
