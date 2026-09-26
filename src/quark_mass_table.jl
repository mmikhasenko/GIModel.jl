# Table-II-style quark masses from the parameters TOML (`[masses]` block).
#
# Public API (exported from GIModel.jl):
#   QuarkMassTable, load_quark_masses, load_parameters_and_quark_masses

"""
    QuarkMassTable

Flavor-keyed constituent masses in GeV. This is an alias for
`Dict{String,Float64}`, with keys `"u"`, `"d"`, `"q"` (light average),
`"s"`, `"c"`, and `"b"`. A manually created table only needs the flavors you use.

## Example

```julia
using GIModel
path = joinpath(pkgdir(GIModel), "data", "parameters.provisional.toml")
mq = load_quark_masses(path)
mq["c"]                           # charm constituent mass in GeV
```

Or create a table and pass it to [`Meson`](@ref):

```julia
custom = QuarkMassTable("c" => 1.628, "b" => 4.977)
Meson(custom, :c, :b)
```

## Related

- [`load_quark_masses`](@ref): masses from TOML.
- [`load_parameters_and_quark_masses`](@ref): parameters and masses.
- [`Meson`](@ref): resolve a flavor pair.
- [`ConstituentMasses`](@ref): explicit pair of masses.
"""
const QuarkMassTable = Dict{String,Float64}

function quark_masses_from_raw(raw)::QuarkMassTable
    m = validate_parameter_section(raw, "masses", ("m_ud_avg_MeV", "m_s_MeV", "m_c_MeV", "m_b_MeV"))
    all(>(0), values(m)) || throw(ArgumentError("constituent masses must be positive"))
    return Dict{String,Float64}(
        "u" => m["m_ud_avg_MeV"] / 1000,
        "d" => m["m_ud_avg_MeV"] / 1000,
        "q" => m["m_ud_avg_MeV"] / 1000,
        "s" => m["m_s_MeV"] / 1000,
        "c" => m["m_c_MeV"] / 1000,
        "b" => m["m_b_MeV"] / 1000,
    )
end

"""
    load_quark_masses(path) -> QuarkMassTable

Read the `[masses]` block of a parameter TOML file, converting MeV to GeV.
The light-average mass supplies all three keys `"u"`, `"d"`, and `"q"`.

## Example

```julia
using GIModel
path = joinpath(pkgdir(GIModel), "data", "parameters.provisional.toml")
mq = load_quark_masses(path)
Meson(mq, :c, :b)
```

## Related

[`QuarkMassTable`](@ref), [`Meson`](@ref), [`load_parameters_and_quark_masses`](@ref).
"""
function load_quark_masses(path::AbstractString)::QuarkMassTable
    return quark_masses_from_raw(TOML.parsefile(path))
end

"""
    load_parameters_and_quark_masses(path) -> (GIParameters, QuarkMassTable)

Read model parameters and constituent masses from one TOML file.
Destructure the returned tuple as `params, mq`.

## Example

```julia
using GIModel
path = joinpath(pkgdir(GIModel), "data", "parameters.provisional.toml")
params, mq = load_parameters_and_quark_masses(path)
meson = Meson(mq, :c, :b)
```

## Related

- [`load_parameters`](@ref): only [`GIParameters`](@ref).
- [`load_quark_masses`](@ref): only the [`QuarkMassTable`](@ref).
- [`Meson`](@ref): choose a channel.
- [`compute_spectrum`](@ref): calculate its levels.
"""
function load_parameters_and_quark_masses(path::AbstractString)
    raw = TOML.parsefile(path)
    return gi_parameters_from_raw(raw), quark_masses_from_raw(raw)
end
