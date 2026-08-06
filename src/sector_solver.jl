# Spin-independent radial sector: cache keys, typed FD solutions, diagnostic equal-mass sweeps.
#
# Public API (exported from GIModel.jl):
#   RadialChannelKey, ChannelRadialSolution, SectorComputation, solve_sector
#   RadialWaveOnUniformMesh(solution::ChannelRadialSolution, radial_level)
#
# [`compute_spectrum`](@ref) (`spectrum.jl`) fills [`SectorComputation`](@ref) with one
# solve per distinct orbital channel.

"""
    RadialChannelKey(masses, L_label)
    RadialChannelKey(m1_GeV, m2_GeV, L_label)

Dict key for one **spin-independent radial channel**: [`ConstituentMasses`](@ref) and orbital
label (`S`, `P`, …), same convention as `SpectrumState.L`. The key normalizes
masses to 12 significant digits so numerically identical channels collapse to
one cache entry without rounding the physics values in [`ConstituentMasses`](@ref).

The three-argument form builds [`ConstituentMasses`](@ref)`(m1_GeV, m2_GeV)` for convenience.
"""
struct RadialChannelKey
    m1_GeV::Float64
    m2_GeV::Float64
    L_label::String
    function RadialChannelKey(masses::ConstituentMasses, L_label::AbstractString)
        new(
            round(masses.m1_GeV; sigdigits = 12),
            round(masses.m2_GeV; sigdigits = 12),
            String(L_label),
        )
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
(shared mesh `solution.r`, spacing ``h = r_2 - r_1``). Use this in spectrum/comparison
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
    SectorComputation(params, solver, channel_cache)

Container filled by [`central_spectrum`](@ref) (`spectrum.jl`): precomputed radial
solves per distinct channel, plus how they were produced.

  - `params`: [`GIParameters`](@ref) used to build each central Hamiltonian.
  - `solver`: the [`RadialSolver`](@ref) that produced them. Later stages resolve
    their own solves (the non-perturbative contact term) and must use the same
    method — otherwise stage 1 and stage 2 of one spectrum would disagree about
    which calculation this is.
  - `channel_cache`: map `RadialChannelKey` → `ChannelRadialSolution`.

[`Spectrum`](@ref) keeps this alive so later stages and two-meson flavor mixing
can reuse the cached solves.
"""
struct SectorComputation{P<:GIParameters,S<:RadialSolver}
    params::P
    solver::S
    channel_cache::Dict{RadialChannelKey,ChannelRadialSolution}
end

"""
    solve_sector(params, equal_mass_GeV; …)

Equal-mass diagnostic sweep over every orbital letter in `L_SYMBOLS`, using `equal_mass_GeV`
for both constituents (same reduced dynamics as charmonium/bottomonium with that mass).
"""
function solve_sector(
    params::GIParameters,
    equal_mass_GeV::Real;
    maxn::Integer = 6,
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    m = Float64(equal_mass_GeV)
    mm = ConstituentMasses(m, m)
    results = Dict{Tuple{Int,String},Float64}()
    for (symbol, L) in L_SYMBOLS
        levels = solve_channel(params, mm, L; nlevels = maxn, solver = solver)
        for n = 1:length(levels)
            results[(n, symbol)] = levels[n]
        end
    end
    return results
end
