# Spin-independent radial sector: cache keys, representation-independent
# solutions, and diagnostic equal-mass sweeps.
#
# Public API (exported from GIModel.jl):
#   RadialChannelKey, ChannelRadialSolution, SectorComputation, solve_sector
#
# [`compute_spectrum`](@ref) (`spectrum.jl`) fills [`SectorComputation`](@ref) with one
# solve per distinct orbital channel.

"""
    RadialChannelKey(masses, L_label[, multiplicity, J])

Dict key for a radial channel: [`ConstituentMasses`](@ref), orbital label, and
optionally spin multiplicity/J. `(multiplicity,J) == (0,0)` identifies a central
spin-independent channel; physical fixed sectors use their actual values. The key normalizes
masses to 12 significant digits so numerically identical channels collapse to
one cache entry without rounding the physics values in [`ConstituentMasses`](@ref).

The three-argument form builds [`ConstituentMasses`](@ref)`(m1_GeV, m2_GeV)` for convenience.
"""
struct RadialChannelKey
    m1_GeV::Float64
    m2_GeV::Float64
    L_label::String
    multiplicity::Int
    J::Int
    function RadialChannelKey(
        masses::ConstituentMasses,
        L_label::AbstractString,
        multiplicity::Integer = 0,
        J::Integer = 0,
    )
        new(
            round(masses.m1_GeV; sigdigits = 12),
            round(masses.m2_GeV; sigdigits = 12),
            String(L_label),
            Int(multiplicity),
            Int(J),
        )
    end
end

RadialChannelKey(m1::Real, m2::Real, L_label::AbstractString) =
    RadialChannelKey(ConstituentMasses(m1, m2), L_label)

RadialChannelKey(
    m1::Real,
    m2::Real,
    L_label::AbstractString,
    multiplicity::Integer,
    J::Integer,
) = RadialChannelKey(ConstituentMasses(m1, m2), L_label, multiplicity, J)

"""
    ChannelRadialSolution(eigenvalues_GeV, waves)
    ChannelRadialSolution(eigenvalues_GeV, mesh_eigenvectors, r)

Output of one `channel_solution` call, stored in `SectorComputation.channel_cache`:

  - `eigenvalues_GeV`: lowest radial eigenvalues in GeV;
  - `waves`: one native [`RadialWave`](@ref) per eigenvalue (`MeshWave` for FD,
    `OscillatorWave` for native HO).

Use [`radial_wave`](@ref)`(solution, n)` to retrieve it. The three-argument
constructor is a convenience for numerical methods that naturally produce a
matrix of mesh samples; the matrix and grid are immediately folded into
`MeshWave` objects and are not retained separately.
"""
struct ChannelRadialSolution{W<:RadialWave}
    eigenvalues_GeV::Vector{Float64}
    waves::Vector{W}
    function ChannelRadialSolution(
        eigenvalues_GeV::AbstractVector{<:Real},
        waves::AbstractVector{W},
    ) where {W<:RadialWave}
        length(eigenvalues_GeV) == length(waves) || throw(ArgumentError(
            "ChannelRadialSolution: eigenvalue/wave counts differ",
        ))
        phased = [fix_outer_phase(wave) for wave in waves]
        return new{W}(collect(Float64, eigenvalues_GeV), collect(W, phased))
    end
end

function ChannelRadialSolution(
    eigenvalues_GeV::AbstractVector{<:Real},
    eigenvectors::AbstractMatrix{<:Real},
    r::AbstractVector{<:Real},
)
    size(eigenvectors, 2) == length(eigenvalues_GeV) || throw(ArgumentError(
        "ChannelRadialSolution: eigenvalue/eigenvector counts differ",
    ))
    waves = [MeshWave(view(eigenvectors, :, n), r) for n in axes(eigenvectors, 2)]
    return ChannelRadialSolution(eigenvalues_GeV, waves)
end

"""Return one native radial eigenlevel from a channel solution."""
function radial_wave(sol::ChannelRadialSolution, radial_level::Integer)
    1 <= radial_level <= length(sol.waves) || throw(ArgumentError(
        "radial_level=$radial_level outside stored range 1:$(length(sol.waves))",
    ))
    return sol.waves[radial_level]
end

"""
    SectorComputation(params, solver, channel_cache)

Container filled by one spectrum calculation: its native radial solves plus how
they were produced. A production `fixed_spectrum` stores fixed `(L,S,J)` keys;
an independent `central_spectrum` stores central `L` keys. The two calculations
are not combined implicitly.

  - `params`: [`GIParameters`](@ref) used to build each central Hamiltonian.
  - `solver`: the [`RadialSolver`](@ref) that produced them and is reused for
    every fixed-sector solve.
  - `channel_cache`: map `RadialChannelKey` → `ChannelRadialSolution`.

[`Spectrum`](@ref) keeps this alive so later stages and two-meson flavor mixing
can reuse the relevant solved states.
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
