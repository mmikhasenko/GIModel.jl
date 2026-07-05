# Core layout types: constituent masses, spin multiplet labels, radial FD samples.
# Experimental/catalog row types live in the GIPaper comparison package.

"""
    ConstituentMasses(m1_GeV, m2_GeV)

Constituent quark masses in GeV for a meson **sector** / spin-independent radial channel.

Components are rounded to 12 significant digits so values stay aligned with
[`RadialChannelKey`](@ref) cache keys and redundant FD solves collapse as intended.

This is plain physics input (not a subset of [`GIParameters`](@ref)); masses appear in
spin-dependent operators, smearing widths, and kinetic factors throughout the package.
"""
struct ConstituentMasses
    m1_GeV::Float64
    m2_GeV::Float64
    function ConstituentMasses(m1::Real, m2::Real)
        new(
            round(Float64(m1); sigdigits = 12),
            round(Float64(m2); sigdigits = 12),
        )
    end
end

"""
    reduced_mass(masses::ConstituentMasses)

Reduced mass ``μ = m_1 m_2 / (m_1 + m_2)`` in GeV (same units as the constituents).
"""
function reduced_mass(m::ConstituentMasses)
    m.m1_GeV * m.m2_GeV / (m.m1_GeV + m.m2_GeV)
end

"""
    FineStructureMultiplet(L_label, multiplicity, J)

Spectroscopic spin/orbital labels shared by spin-dependent corrections:

  - [`fine_structure_components`](@ref) uses `L_label`, `multiplicity`, and `J`.
  - [`contact_hyperfine_shift`](@ref) uses only `L_label` and `multiplicity` (`J` is ignored).

The comparison layer (GIPaper) adds an overload for its reference-catalog rows.
"""
struct FineStructureMultiplet
    L_label::String
    multiplicity::Int
    J::Int
    function FineStructureMultiplet(
        L_label::AbstractString,
        multiplicity::Integer,
        J::Integer,
    )
        new(String(L_label), Int(multiplicity), Int(J))
    end
end

"""
    RadialWaveOnUniformMesh(u, r, h)
    RadialWaveOnUniformMesh(u, r)

Reduced radial wavefunction ``u(r)`` on a **uniform** interior grid: samples `uᵢ` and radii `rᵢ`
with spacing `h` (for `length(r) ≥ 2`, the two-argument form sets `h = r[2] - r[1]`).

This bundles the data [`fine_structure_components`](@ref), [`contact_hyperfine_shift`](@ref),
and [`physical_u_norm`](@ref) rely on. It is **one radial eigenlevel** on the mesh — not the
full multi-level output of [`channel_solution`](@ref).

For the cached workflow object [`ChannelRadialSolution`](@ref), use the constructor
`RadialWaveOnUniformMesh(solution, radial_level)` defined in `sector_solver.jl`: it takes
column `radial_level` of `solution.eigenvectors` together with `solution.r`.

The explicit `h` argument must agree with the uniform spacing implied by `r` (guardrail).
"""
struct RadialWaveOnUniformMesh
    u::Vector{Float64}
    r::Vector{Float64}
    h::Float64
    function RadialWaveOnUniformMesh(
        u::AbstractVector{<:Real},
        r::AbstractVector{<:Real},
        h::Real,
    )
        length(u) == length(r) ||
            throw(ArgumentError("RadialWaveOnUniformMesh: length(u) != length(r)"))
        length(r) >= 1 || throw(ArgumentError("RadialWaveOnUniformMesh: empty r"))
        hf = float(h)
        isfinite(hf) && hf > 0 ||
            throw(ArgumentError("RadialWaveOnUniformMesh: invalid mesh spacing h=$h"))
        if length(r) >= 2
            hinfer = float(r[2]) - float(r[1])
            isapprox(hinfer, hf; rtol = 1e-10, atol = 1e-12) ||
                throw(
                    ArgumentError(
                        "RadialWaveOnUniformMesh: h=$hf inconsistent with r spacing $hinfer",
                    ),
                )
        end
        return new(collect(Float64, u), collect(Float64, r), hf)
    end
end

function RadialWaveOnUniformMesh(u::AbstractVector{<:Real}, r::AbstractVector{<:Real})
    length(r) >= 2 ||
        throw(ArgumentError("RadialWaveOnUniformMesh(u,r): need length(r) ≥ 2 to infer h"))
    return RadialWaveOnUniformMesh(u, r, r[2] - r[1])
end
