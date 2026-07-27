# Quark types. The hierarchy encodes which flavor distinctions the GI model can
# actually resolve — see the `AbstractQuark` docstring.
#
# Public API (exported from GIModel.jl):
#   AbstractQuark, LightQuark, StrangeQuark, HeavyQuark, charge, flavor_symbol, mass_GeV

"""
    AbstractQuark

A constituent quark: a **mass** plus whatever flavor identity the model can
resolve. In the GI model the medium (`b`, `c`, `σ₀`, `s`, `α_s(Q²)`) is
flavor-blind and the mass is the only dynamical input — quarks of equal mass
give bit-identical spectra regardless of name. Flavor identity therefore does
exactly two jobs: it carries **electric charge**, and it supplies the name whose
equality gates self-conjugate behavior ([`is_equal_flavor`](@ref)).

The three concrete types record *which distinctions are physical here*:

| type | mass | [`charge`](@ref) |
|:--|:--|:--|
| [`LightQuark`](@ref)   | `m_u = m_d` | **undefined** — the model has no isospin breaking |
| [`StrangeQuark`](@ref) | `m_s`       | `-1//3` |
| [`HeavyQuark`](@ref)   | free dial   | `+2//3` / `-1//3` from the type tag |

Every subtype has a `mass_GeV::Float64` field; use [`mass_GeV`](@ref) to read it
generically and [`flavor_symbol`](@ref) for the label a [`Meson`](@ref) stores.
"""
abstract type AbstractQuark end

"""
    LightQuark(mass_GeV)

An up/down constituent quark. The GI parameter set assigns `u` and `d` the *same*
mass (`m_ud_avg`), so the model cannot resolve them: nothing in the dynamics
branches on up-vs-down.

Consequently [`charge`](@ref) is **deliberately not defined** for `LightQuark` —
calling it raises `MethodError` rather than returning a wrong number. A light
quark has no charge of its own here: charge belongs to the meson's flavor
wavefunction, which distinguishes states the constituent masses cannot
(`π⁺ = u d̄` carries charge while `π⁰ = (u ū − d d̄)/√2` does not, at identical
mass). Callers that need light-quark charges supply them per meson row.
"""
struct LightQuark <: AbstractQuark
    mass_GeV::Float64
    function LightQuark(mass::Real)
        mass > 0 || throw(ArgumentError("quark mass must be positive, got $mass"))
        new(Float64(mass))
    end
end

"""
    StrangeQuark(mass_GeV)

The strange constituent quark. It gets its own type because it is the model's
SU(3)-breaking axis: `m_s ≠ m_ud` splits the strange sectors, and the
nonstrange/strange (`ns`/`ss`) distinction drives the isoscalar annihilation
basis. Unlike [`LightQuark`](@ref) it has no mass-degenerate isospin partner, so
its charge is unambiguous: `-1//3`.
"""
struct StrangeQuark <: AbstractQuark
    mass_GeV::Float64
    function StrangeQuark(mass::Real)
        mass > 0 || throw(ArgumentError("quark mass must be positive, got $mass"))
        new(Float64(mass))
    end
end

"""
    HeavyQuark{T}(mass_GeV, name)

A heavy constituent quark, where `T` is the weak-isospin class — `:up` for
charm/top (`+2//3`) or `:down` for bottom (`-1//3`). Two charge classes genuinely
exist among the heavy flavors, so the tag is what distinguishes them; the `name`
(`:c`, `:b`, `:t`, …) is bookkeeping and is dynamically inert.

`mass_GeV` is a **free external dial**: heavy quarks decouple from the light
sector (changing `m_c` leaves the ρ bit-identical), so it can be varied without
retuning the medium. Heavy-light observables then organize into heavy-quark
symmetry — `Λ̄ = M̄ − m_Q` and `ΔM_hf · m_Q` both flatten as `m_Q` grows.

    HeavyQuark{:up}(1.628, :c)     # charm
    HeavyQuark{:down}(4.977, :b)   # bottom

!!! warning
    The FD radial grid and the fixed harmonic-oscillator β grid are tuned for
    masses up to bottomonium; quarkonium much heavier than `m_b` becomes too
    compact to resolve on the defaults.
"""
struct HeavyQuark{T} <: AbstractQuark
    mass_GeV::Float64
    name::Symbol
    function HeavyQuark{T}(mass::Real, name::Symbol) where {T}
        T in (:up, :down) || throw(ArgumentError(
            "HeavyQuark tag must be :up or :down (weak-isospin class), got `$T`",
        ))
        mass > 0 || throw(ArgumentError("quark mass must be positive, got $mass"))
        new{T}(Float64(mass), name)
    end
end

"""
    mass_GeV(q::AbstractQuark) -> Float64

Constituent mass in GeV — the only dynamical input the GI medium responds to.
"""
mass_GeV(q::AbstractQuark) = q.mass_GeV

"""
    flavor_symbol(q::AbstractQuark) -> Symbol

The flavor label a [`Meson`](@ref) stores for this quark: `:q` for light (the
isospin-averaged stand-in), `:s` for strange, and the given `name` for heavy.
Equality of these symbols is what [`is_equal_flavor`](@ref) tests.
"""
flavor_symbol(::LightQuark) = :q
flavor_symbol(::StrangeQuark) = :s
flavor_symbol(q::HeavyQuark) = q.name

"""
    charge(q::AbstractQuark) -> Rational{Int}

Electric charge in units of `e`. Defined for [`StrangeQuark`](@ref) (`-1//3`) and
[`HeavyQuark`](@ref) (`+2//3` for `:up`, `-1//3` for `:down`).

**Not defined for [`LightQuark`](@ref)** — see its docstring; the model cannot
resolve `u` from `d`, so a light quark has no charge of its own and the call
raises `MethodError` instead of guessing.
"""
charge(::StrangeQuark) = -1 // 3
charge(::HeavyQuark{:up}) = 2 // 3
charge(::HeavyQuark{:down}) = -1 // 3
