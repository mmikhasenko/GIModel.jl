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

A ``q_1 \\bar q_2`` channel with constituent masses. Flavors are `:u`, `:d`,
`:s`, `:c`, `:b`, or `:q` (light average; `:n` is an alias).
The antiquark flavor is unadorned: `:c, :u` means charm–antiup.

Supply a [`QuarkMassTable`](@ref) or explicit [`ConstituentMasses`](@ref).
Unknown or missing table flavors throw `ArgumentError`.

## Example

```julia
using GIModel
mq = QuarkMassTable("c" => 1.628, "b" => 4.977)
meson = Meson(mq, :c, :b)
```

Read flavors through `meson.flavor1` and `meson.flavor2`; masses through
`meson.constituent_masses.m1_GeV` and `.m2_GeV`.

```julia
reduced_mass(meson)  # reduced mass in GeV
flavor_label(meson)  # channel label, e.g. "cb"
is_equal_flavor(meson)
propertynames(meson)
```

For bound-state masses, supply `params` from [`load_parameters`](@ref):

```julia
spec = compute_spectrum(params, meson; levels=spectrum_levels(1))
state = spectrum_state(spec, "1^1P_1")
state.mass_GeV
```

## Related

- [`reduced_mass`](@ref): reduced mass in GeV.
- [`flavor_label`](@ref): channel label.
- [`is_equal_flavor`](@ref): flavor symmetry.
- [`load_quark_masses`](@ref): masses from TOML.
- [`load_parameters`](@ref): interaction parameters.
- [`spectrum_levels`](@ref): level selection.
- [`compute_spectrum`](@ref): meson spectrum.
- [`spectrum_state`](@ref): one computed level.
- [`physical_components`](@ref): wavefunction components.
"""
struct Meson{M<:ConstituentMasses}
    flavor1::Symbol
    flavor2::Symbol
    constituent_masses::M
    function Meson(flavor1::Symbol, flavor2::Symbol, masses::M) where {M<:ConstituentMasses}
        new{M}(_canonical_flavor(flavor1), _canonical_flavor(flavor2), masses)
    end
end

"""
    Meson(q1::AbstractQuark, q2::AbstractQuark)

Build a [`Meson`](@ref) from two quarks, using their masses and
[`flavor_symbol`](@ref). The second quark specifies the antiquark flavor.

## Example

```julia
using GIModel
Meson(HeavyQuark{:up}(1.628, :c), HeavyQuark{:up}(1.628, :c))
```

Inspect `.constituent_masses`, then use [`compute_spectrum`](@ref) with
interaction parameters to calculate bound-state masses.

## Related

- [`reduced_mass`](@ref): reduced mass in GeV.
- [`flavor_label`](@ref): channel label.
- [`is_equal_flavor`](@ref): flavor symmetry.
- [`load_parameters`](@ref): interaction parameters.
"""
function Meson(q1::AbstractQuark, q2::AbstractQuark)
    return Meson(
        flavor_symbol(q1),
        flavor_symbol(q2),
        ConstituentMasses(mass_GeV(q1), mass_GeV(q2)),
    )
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
Different flavors with manually equal masses still return `false`.

## Related

- [`flavor_label`](@ref): channel label for a [`Meson`](@ref).
- [`compute_spectrum`](@ref): spectrum with allowed mixing.
"""
is_equal_flavor(m::Meson) = m.flavor1 == m.flavor2

"""
    reduced_mass(meson::Meson)

Return `m1 * m2 / (m1 + m2)` in GeV from the [`Meson`](@ref)'s constituent masses.
This is a two-body input quantity, not the predicted bound-state mass.

## Related

- [`compute_spectrum`](@ref): bound-state masses.
- [`spectrum_state`](@ref): select a predicted level.
"""
reduced_mass(m::Meson) = reduced_mass(m.constituent_masses)

"""
    flavor_label(meson::Meson) -> String

Return the compact flavor-pair label, e.g. `"cb"` for a charm–antibottom
[`Meson`](@ref). Read `meson.flavor1` and `meson.flavor2` for individual symbols.

## Related

- [`is_equal_flavor`](@ref): flavor symmetry.
- [`compute_spectrum`](@ref): levels for this channel.
"""
flavor_label(m::Meson) = string(m.flavor1, m.flavor2)

function Base.show(io::IO, ::MIME"text/plain", m::Meson)
    print(
        io, "Meson: ", m.flavor1, " ", m.flavor2, "bar   (m1 = ",
        m.constituent_masses.m1_GeV, ", m2 = ", m.constituent_masses.m2_GeV,
        " GeV, reduced = ", round(reduced_mass(m), digits = 4), " GeV)",
    )
    is_equal_flavor(m) && print(io, "\n  self-conjugate: same-J antisymmetric spin-orbit mixing vanishes")
    return nothing
end

Base.show(io::IO, m::Meson) = print(io, "Meson(", m.flavor1, " ", m.flavor2, "bar)")
