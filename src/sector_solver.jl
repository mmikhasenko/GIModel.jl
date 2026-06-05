# Spin-independent radial sector: cache keys, typed FD solutions, diagnostic equal-mass sweeps.
#
# Public API (exported from GIModel.jl):
#   RadialChannelKey, ChannelRadialSolution, SectorComputation, solve_sector
#   RadialWaveOnUniformMesh(solution::ChannelRadialSolution, radial_level)
#
# Batching reference rows into [`SectorComputation`](@ref) is [`compute_sector`](@ref) in
# `sector_comparison.jl`.

"""
    RadialChannelKey(masses, L_label)
    RadialChannelKey(m1_GeV, m2_GeV, L_label)

Dict key for one **spin-independent radial channel**: [`ConstituentMasses`](@ref) and orbital
label (`S`, `P`, …), same convention as `ReferenceState.L`. Mass rounding lives in
[`ConstituentMasses`](@ref) so distinct finite-difference solves that share the same physics
collapse to one cache entry.

The three-argument form builds [`ConstituentMasses`](@ref)`(m1_GeV, m2_GeV)` for convenience.
"""
struct RadialChannelKey
    masses::ConstituentMasses
    L_label::String
    function RadialChannelKey(masses::ConstituentMasses, L_label::AbstractString)
        new(masses, String(L_label))
    end
end

RadialChannelKey(m1::Real, m2::Real, L_label::AbstractString) =
    RadialChannelKey(ConstituentMasses(m1, m2), L_label)

"""
    ChannelRadialSolution(eigenvalues_GeV, eigenvectors, r)

Output of one `channel_solution` call, stored in `SectorComputation.channel_cache`:

  - `eigenvalues_GeV`: lowest radial eigenvalues (GeV) of the central Hamiltonian on the mesh.
  - `eigenvectors`: columns are reduced radial functions ``u_n(r)`` for each level.
  - `r`: uniform interior radial grid (same spacing as in the FD builder).

Spin-dependent expectations ([`fine_structure_components`](@ref), contact hyperfine) need a
**single level** ``u_n`` on this mesh — see [`RadialWaveOnUniformMesh`](@ref)`(solution, n)` below,
which wraps column `n` with the correct spacing `h`.
"""
struct ChannelRadialSolution
    eigenvalues_GeV::Vector{Float64}
    eigenvectors::Matrix{Float64}
    r::Vector{Float64}
end

"""
    RadialWaveOnUniformMesh(solution::ChannelRadialSolution, radial_level::Integer)

Build [`RadialWaveOnUniformMesh`](@ref) for eigenvector column `radial_level` of `solution`
(shared mesh `solution.r`, spacing ``h = r_2 - r_1``). Use this in [`compare`](@ref) /
tooling when you already hold cached [`ChannelRadialSolution`](@ref) data instead of raw `(u, r)` vectors.
"""
function RadialWaveOnUniformMesh(sol::ChannelRadialSolution, radial_level::Integer)
    radial_level >= 1 ||
        throw(ArgumentError("radial_level must be ≥ 1, got $radial_level"))
    radial_level <= size(sol.eigenvectors, 2) ||
        throw(ArgumentError(
            "radial_level=$radial_level exceeds number of stored eigenvectors $(size(sol.eigenvectors, 2))",
        ))
    length(sol.r) >= 2 ||
        throw(ArgumentError("ChannelRadialSolution.r must have length ≥ 2"))
    ucol = view(sol.eigenvectors, :, radial_level)
    return RadialWaveOnUniformMesh(ucol, sol.r)
end

"""
    SectorComputation(params, channel_cache)

Container filled by [`compute_sector`](@ref) (`sector_comparison.jl`): precomputed radial FD solves per distinct channel.

  - `params`: [`GIParameters`](@ref) used to build each central Hamiltonian.
  - `channel_cache`: map `RadialChannelKey` → `ChannelRadialSolution`.

[`compare`](@ref) expects the same [`ReferenceStateWithMasses`](@ref) rows used to build the cache.
"""
struct SectorComputation
    params::GIParameters
    channel_cache::Dict{RadialChannelKey,ChannelRadialSolution}
    ho_wave_cache::Dict{RadialChannelKey,ChannelRadialSolution}
end

SectorComputation(params, cache) =
    SectorComputation(params, cache, Dict{RadialChannelKey,ChannelRadialSolution}())

"""
    solve_sector(params, equal_mass_GeV; …)

Equal-mass diagnostic sweep over every orbital letter in `L_SYMBOLS`, using `equal_mass_GeV`
for both constituents (same reduced dynamics as charmonium/bottomonium with that mass).
"""
function solve_sector(
    params::GIParameters,
    equal_mass_GeV::Real;
    maxn::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
)
    m = Float64(equal_mass_GeV)
    mm = ConstituentMasses(m, m)
    results = Dict{Tuple{Int,String},Float64}()
    for (symbol, L) in L_SYMBOLS
        levels = solve_channel(
            params,
            mm,
            L;
            nlevels = maxn,
            ngrid = ngrid,
            rmax = rmax,
            kinetic = kinetic,
            eigensolver = eigensolver,
        )
        for n = 1:length(levels)
            results[(n, symbol)] = levels[n]
        end
    end
    return results
end
