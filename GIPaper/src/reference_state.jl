# Experimental/catalog row types for reference spectra + CSV loader (`reference_spectrum_*.csv`).
#
# Public API (exported from GIPaper.jl): ReferenceState, load_reference_spectrum

import GIModel: FineStructureMultiplet

"""
    ReferenceState

One row from a reference spectrum CSV (`sector`, `quark_content`, radial/spin labels,
experimental `mass_GeV`, etc.). Loaded via [`load_reference_spectrum`](@ref), mapped to a
[`GIModel.Meson`](@ref) by [`reference_meson`](@ref), then compared via [`compare_reference`](@ref).
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
    FineStructureMultiplet(state::ReferenceState)

Copy `(L, multiplicity, J)` from a reference row into [`GIModel.FineStructureMultiplet`](@ref).
"""
FineStructureMultiplet(state::ReferenceState) =
    FineStructureMultiplet(state.L, state.multiplicity, state.J)

"""
    load_reference_spectrum(path::AbstractString) -> Vector{ReferenceState}

Read a GI-style reference spectrum CSV. Required columns match
the data-invariant tests under `GIPaper/test/`.

See [`ReferenceState`](@ref); use [`reference_spectrum_path`](@ref) for the bundled catalogs.
"""
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
