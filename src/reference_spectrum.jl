# Experimental / catalog rows and comparison-side wrappers (CSV → structs).
# Not part of the spin-independent Hamiltonian or channel solver core.

"""
    ReferenceState

One row from a reference spectrum CSV (`sector`, `quark_content`, radial/spin labels,
experimental `mass_GeV`, etc.). Used with [`load_reference_spectrum`](@ref), then optionally
[`attach_constituent_masses`](@ref), before [`compute_sector`](@ref) / [`compare`](@ref).
"""
struct ReferenceState
    sector::String
    quark_content::String
    composition::String
    n::Int
    multiplicity::Int
    L::String
    J::Int
    mass_GeV::Float64
    confidence::String
end

"""
    ReferenceStateWithMasses(state, constituent_masses)

Reference CSV row plus [`ConstituentMasses`](@ref) resolved **once** at setup.

Use [`attach_constituent_masses`](@ref) / [`resolve_constituent_masses`](@ref) with a
[`QuarkMassTable`](@ref); pass the vector to [`compute_sector`](@ref) / [`compare`](@ref).

For scans with arbitrary masses, construct manually:
`ReferenceStateWithMasses(state, ConstituentMasses(m1, m2))`.
"""
struct ReferenceStateWithMasses
    state::ReferenceState
    constituent_masses::ConstituentMasses
end

"""
    FineStructureMultiplet(state::ReferenceState)
    FineStructureMultiplet(row::ReferenceStateWithMasses)

Copy `(L, multiplicity, J)` from a reference row into [`FineStructureMultiplet`](@ref).
"""
FineStructureMultiplet(state::ReferenceState) =
    FineStructureMultiplet(state.L, state.multiplicity, state.J)

FineStructureMultiplet(row::ReferenceStateWithMasses) = FineStructureMultiplet(row.state)

function load_reference_spectrum(path::AbstractString)
    states = ReferenceState[]
    for row in CSV.File(path)
        push!(
            states,
            ReferenceState(
                String(row.sector),
                String(row.quark_content),
                String(row.composition_raw),
                Int(row.n),
                Int(row.multiplicity),
                String(row.L),
                Int(row.J),
                Float64(row.mass_GeV),
                String(row.confidence),
            ),
        )
    end
    states
end
