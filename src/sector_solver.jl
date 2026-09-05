# Spin-independent radial sector: cache keys, representation-independent
# solutions, and diagnostic equal-mass sweeps.
#
# Public API (exported from GIModel.jl):
#   RadialChannelKey, OscillatorConvergence, ChannelRadialSolution,
#   SectorComputation, solve_sector
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
    OscillatorConvergence

Numerical certificate attached to an HO [`ChannelRadialSolution`](@ref).
`status` is `:converged` for an adaptive solve or `:unchecked` when the caller
explicitly requested `OscillatorSolver(converge=false)`. `energy_delta_GeV` is
the largest requested-eigenvalue change on the final basis refinement;
`max_wave_overlap_defect` is `max(1-|<u_N|u_previous>|)` over those levels.
The latter is recorded for wave-sensitive follow-up audits but is not used as
an energy-convergence substitute.
"""
struct OscillatorConvergence
    status::Symbol
    beta_GeV::Float64
    nbasis::Int
    energy_delta_GeV::Union{Nothing,Float64}
    max_wave_overlap_defect::Union{Nothing,Float64}
    tolerance_GeV::Float64
    refinements::Int
    function OscillatorConvergence(
        status::Symbol,
        beta_GeV::Real,
        nbasis::Integer,
        energy_delta_GeV::Union{Nothing,Real},
        max_wave_overlap_defect::Union{Nothing,Real},
        tolerance_GeV::Real,
        refinements::Integer,
    )
        status in (:converged, :unchecked) || throw(ArgumentError(
            "OscillatorConvergence: status must be :converged or :unchecked",
        ))
        beta_GeV > 0 || throw(ArgumentError(
            "OscillatorConvergence: beta must be positive",
        ))
        nbasis >= 1 || throw(ArgumentError(
            "OscillatorConvergence: nbasis must be positive",
        ))
        tolerance_GeV > 0 || throw(ArgumentError(
            "OscillatorConvergence: tolerance must be positive",
        ))
        refinements >= 0 || throw(ArgumentError(
            "OscillatorConvergence: refinements must be non-negative",
        ))
        delta = isnothing(energy_delta_GeV) ? nothing : float(energy_delta_GeV)
        defect = isnothing(max_wave_overlap_defect) ? nothing :
                 float(max_wave_overlap_defect)
        status == :converged && (isnothing(delta) || isnothing(defect)) &&
            throw(ArgumentError(
                "OscillatorConvergence: converged status requires final deltas",
            ))
        !isnothing(delta) && (!isfinite(delta) || delta < 0) && throw(ArgumentError(
            "OscillatorConvergence: energy delta must be finite and non-negative",
        ))
        !isnothing(defect) && (!isfinite(defect) || defect < 0) && throw(ArgumentError(
            "OscillatorConvergence: overlap defect must be finite and non-negative",
        ))
        return new(
            status,
            float(beta_GeV),
            Int(nbasis),
            delta,
            defect,
            float(tolerance_GeV),
            Int(refinements),
        )
    end
end

"""
    ChannelRadialSolution(eigenvalues_GeV, waves)
    ChannelRadialSolution(eigenvalues_GeV, mesh_eigenvectors, r)

Output of one `channel_solution` call, stored in `SectorComputation.channel_cache`:

  - `eigenvalues_GeV`: lowest radial eigenvalues in GeV;
  - `waves`: one native [`RadialWave`](@ref) per eigenvalue (`MeshWave` for FD,
    `OscillatorWave` for native HO);
  - `convergence`: the [`OscillatorConvergence`](@ref) certificate for HO, or
    `nothing` for a method without an automatic certificate.

Use [`radial_wave`](@ref)`(solution, n)` to retrieve it. The three-argument
constructor is a convenience for numerical methods that naturally produce a
matrix of mesh samples; the matrix and grid are immediately folded into
`MeshWave` objects and are not retained separately.
"""
struct ChannelRadialSolution{W<:RadialWave}
    eigenvalues_GeV::Vector{Float64}
    waves::Vector{W}
    convergence::Union{Nothing,OscillatorConvergence}
    function ChannelRadialSolution(
        eigenvalues_GeV::AbstractVector{<:Real},
        waves::AbstractVector{W},
        ;
        convergence::Union{Nothing,OscillatorConvergence} = nothing,
    ) where {W<:RadialWave}
        length(eigenvalues_GeV) == length(waves) || throw(ArgumentError(
            "ChannelRadialSolution: eigenvalue/wave counts differ",
        ))
        phased = [fix_outer_phase(wave) for wave in waves]
        return new{W}(
            collect(Float64, eigenvalues_GeV), collect(W, phased), convergence,
        )
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
