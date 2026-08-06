#!/usr/bin/env julia
# Campaign-style diagnostics for the Table III P1/P2 reproduction mismatch.

using Pkg
# Archived script: runs in the GIPaper environment at the repository root.
Pkg.activate(joinpath(@__DIR__, "..", "..", "GIPaper"))

using CSV
using LinearAlgebra
using Printf

root = joinpath(dirname(dirname(@__DIR__)), "GIPaper")
using GIModel

const TARGET_BASIS = ["1 ns", "1 ss", "1 cc", "1 bb", "2 ns", "2 ss", "2 cc"]
const REPORT = joinpath(@__DIR__, "table_iii_suspect_investigation.md")
const PERT_COEFF = 2π / 3 * (log(2) - 1)

struct BasisProbe
    label::String
    flavor::String
    n::Int
    mass::Float64
    diagonal::Float64
    wave::RadialWaveOnUniformMesh
end

function reduced_p2(wave)
    u = wave.u
    h = wave.h
    norm = sum(abs2, u) * h
    prev = 0.0
    accum = 0.0
    for ui in u
        accum += (ui - prev)^2 / h
        prev = ui
    end
    accum += prev^2 / h
    return accum / norm
end

function origin_R(wave)
    return wave.u[1] / wave.r[1]
end

function alpha_mass(m)
    return GIModel.alpha_s_q(max(m, 1e-9))
end

function load_basis(params, mq; ngrid = 220, rmax = 22.0)
    channels = Dict{String,Any}()
    for (label, mass_key) in [("ns", "q"), ("ss", "s"), ("cc", "c"), ("bb", "b")]
        masses = ConstituentMasses(mq[mass_key], mq[mass_key])
        ev, vecs, r = channel_solution(
            params, masses, 0;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax),
        )
        sol = GIModel.ChannelRadialSolution(ev, vecs, r)
        levels = GIModel.contact_hyperfine_nonperturbative_levels(params, masses, "S", 1, r, 2)
        channels[label] = (mass = mq[mass_key], solution = sol, levels = levels)
    end
    out = BasisProbe[]
    for label in TARGET_BASIS
        n = parse(Int, first(split(label)))
        flavor = last(split(label))
        channel = channels[flavor]
        push!(
            out,
            BasisProbe(
                label,
                flavor,
                n,
                channel.mass,
                channel.levels[n],
                RadialWaveOnUniformMesh(channel.solution, n),
            ),
        )
    end
    return out
end

function table_targets(model)
    rows = collect(CSV.File(joinpath(root, "data", "clean", "mixings.csv")))
    dict = Dict{String,Vector{Float64}}()
    shifts = Dict{String,Union{Missing,Float64}}()
    for row in rows
        ismissing(row.model) && continue
        String(row.model) == model || continue
        state = String(row.state_name)
        vec = get!(dict, state, zeros(Float64, length(TARGET_BASIS)))
        idx = findfirst(==(String(row.basis)), TARGET_BASIS)
        isnothing(idx) && continue
        vec[idx] = Float64(row.amplitude)
        if !ismissing(row.mass_shift_MeV)
            shifts[state] = Float64(row.mass_shift_MeV)
        elseif !haskey(shifts, state)
            shifts[state] = missing
        end
    end
    ordered = sort(collect(keys(dict)))
    return ordered, [dict[state] for state in ordered], shifts
end

function smearing_factor(b::BasisProbe, variant::Symbol)
    input = pseudoscalar_annihilation_basis_input(
        b.label,
        b.mass,
        b.diagonal,
        b.wave,
    )
    R0 = origin_R(b.wave)
    p2 = reduced_p2(b.wave)
    rel = b.mass / sqrt(b.mass^2 + p2)
    if variant == :current_abs
        return abs(R0) * rel / sqrt(4π)
    elseif variant == :signed
        return R0 * rel / sqrt(4π)
    elseif variant == :no_rel
        return abs(R0) / sqrt(4π)
    elseif variant == :no_sqrt4pi
        return abs(R0) * rel
    elseif variant == :origin_abs
        return abs(R0)
    elseif variant == :momentum_integral
        return GIModel._s0_smearing_factor(GIModel.FDMomentumIntegralSmearing(), input)
    elseif variant == :origin_p2_coherent
        return GIModel._flavor_coherence_factor(input) * abs(R0) * rel / sqrt(4π)
    elseif variant == :momentum_integral_coherent
        return GIModel._annihilation_overlap_factor(GIModel.FDMomentumIntegralSmearing(), input)
    else
        error("unknown smearing variant $variant")
    end
end

function p1_bracket(params, left::BasisProbe, right::BasisProbe)
    nonpert = params.annihilation_p1_A_np *
              exp(-(left.mass^2 + right.mass^2) / params.annihilation_p1_m_eta^2)
    pert = PERT_COEFF * alpha_mass(left.diagonal) * alpha_mass(right.diagonal) / π^2
    return nonpert + pert
end

function p1_matrix(params, basis; variant = :current_abs, strength = 1.0)
    diag = [b.diagonal for b in basis]
    f = [smearing_factor(b, variant) for b in basis]
    A = zeros(Float64, length(basis), length(basis))
    for j in eachindex(basis), i in eachindex(basis)
        A[j, i] = 4π * p1_bracket(params, basis[j], basis[i]) * f[j] * f[i] /
                  (basis[j].mass * basis[i].mass)
    end
    return Matrix(Diagonal(diag)) + strength * A, A
end

function align_score(vectors, targets; basis_phase_search = false)
    best = Inf
    best_cols = Int[]
    best_signs = Float64[]
    phase_patterns = basis_phase_search ? Iterators.product(fill((-1.0, 1.0), length(TARGET_BASIS))...) : ((ones(Float64, length(TARGET_BASIS))...,),)
    for pattern_tuple in phase_patterns
        phases = collect(Float64, pattern_tuple)
        transformed = Diagonal(phases) * vectors
        used = Set{Int}()
        scores = Float64[]
        cols = Int[]
        signs = Float64[]
        for target in targets
            col = 0
            dotmax = -Inf
            for candidate in axes(transformed, 2)
                candidate in used && continue
                d = dot(transformed[:, candidate], target)
                if abs(d) > dotmax
                    dotmax = abs(d)
                    col = candidate
                end
            end
            push!(used, col)
            sgn = dot(transformed[:, col], target) < 0 ? -1.0 : 1.0
            diff = sgn .* transformed[:, col] .- target
            push!(scores, sqrt(sum(abs2, diff) / length(diff)))
            push!(cols, col)
            push!(signs, sgn)
        end
        score = sum(scores) / length(scores)
        if score < best
            best = score
            best_cols = cols
            best_signs = signs
        end
    end
    return best, best_cols, best_signs
end

function score_variant(params, basis, targets; variant = :current_abs, strength = 1.0, phase_search = false)
    matrix, A = p1_matrix(params, basis; variant = variant, strength = strength)
    fact = eigen(Symmetric(matrix))
    score, cols, _ = align_score(fact.vectors, targets; basis_phase_search = phase_search)
    return (score = score, masses = fact.values, cols = cols, matrix = matrix, annihilation = A)
end

function best_strength(params, basis, targets; variant = :current_abs)
    grid = 10 .^ range(-2, 4; length = 241)
    best = (score = Inf, strength = NaN, masses = Float64[])
    for strength in grid
        scored = score_variant(params, basis, targets; variant = variant, strength = strength)
        if scored.score < best.score
            best = (score = scored.score, strength = strength, masses = scored.masses)
        end
    end
    return best
end

function closest_strength_for_first_mass(params, basis, targets, target_mass; variant = :current_abs)
    grid = 10 .^ range(-2, 4; length = 241)
    best = (distance = Inf, strength = NaN, score = Inf, mass = NaN)
    for strength in grid
        scored = score_variant(params, basis, targets; variant = variant, strength = strength)
        mass = scored.masses[1]
        distance = abs(mass - target_mass)
        if distance < best.distance
            best = (distance = distance, strength = strength, score = scored.score, mass = mass)
        end
    end
    return best
end

function main()
    params, mq = load_parameters_and_quark_masses(joinpath(dirname(root), "data", "parameters.provisional.toml"))
    basis = load_basis(params, mq)
    states, p1_targets, _ = table_targets("P1")

    current = score_variant(params, basis, p1_targets; variant = :momentum_integral_coherent)
    phase_only = score_variant(params, basis, p1_targets; phase_search = true)
    variants = [
        :current_abs,
        :signed,
        :momentum_integral,
        :origin_p2_coherent,
        :momentum_integral_coherent,
        :no_rel,
        :no_sqrt4pi,
        :origin_abs,
    ]
    variant_scores = [(variant, score_variant(params, basis, p1_targets; variant = variant).score) for variant in variants]
    strength = best_strength(params, basis, p1_targets; variant = :momentum_integral_coherent)
    eta_strength = closest_strength_for_first_mass(
        params,
        basis,
        p1_targets,
        0.520;
        variant = :momentum_integral_coherent,
    )

    open(REPORT, "w") do io
        println(io, "# Table III Suspect Investigation")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/investigate_table_iii_suspects.jl`.")
        println(io)
        println(io, "This campaign probes whether the Table III mismatch behaves like a phase convention issue, an overall annihilation-strength error, a coordinate/momentum `alpha_s` convention issue, an Eq. (17) wavefunction-factor problem, or a missing `ns` flavor-coherence factor.")
        println(io)
        println(io, "## Basis Diagnostics")
        println(io)
        println(io, "| basis | diagonal GeV | R(0) proxy | <p^2> GeV^2 | m/sqrt(m^2+<p^2>) | legacy S0 factor |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for b in basis
            p2 = reduced_p2(b.wave)
            rel = b.mass / sqrt(b.mass^2 + p2)
            println(
                io,
                @sprintf(
                    "| `%s` | %.3f | %.3f | %.3f | %.3f | %.4f |",
                    b.label,
                    b.diagonal,
                    origin_R(b.wave),
                    p2,
                    rel,
                    smearing_factor(b, :current_abs),
                ),
            )
        end
        println(io)
        println(io, "## Suspect Ranking")
        println(io)
        println(io, "| suspect | probe | result | interpretation |")
        println(io, "|---|---|---:|---|")
        println(io, @sprintf("| basis phase convention | exhaustive per-basis phase search | %.3f RMS | phase choices do not rescue the amplitudes |", phase_only.score))
        println(io, @sprintf("| current implementation | momentum Eq. (17) + coherent `ns` FD P1 | %.3f RMS | baseline mismatch |", current.score))
        println(io, @sprintf("| overall strength | best scalar multiplier %.3g | %.3f RMS | scalar normalization helps, but still not a clean reproduction |", strength.strength, strength.score))
        println(io, @sprintf("| mass-only retune | closest first pole %.3f GeV at strength %.3g | %.3f RMS | scalar strength cannot cleanly target the eta mass or composition |", eta_strength.mass, eta_strength.strength, eta_strength.score))
        best_variant = first(sort(variant_scores; by = x -> x[2]))
        println(io, @sprintf("| Eq. (17) / flavor variants | best local variant `%s` | %.3f RMS | local factors alone do not reproduce Table III |", String(best_variant[1]), best_variant[2]))
        println(io)
        println(io, "## Simple Eq. (17) Variants")
        println(io)
        println(io, "| variant | mean P1 vector RMS |")
        println(io, "|---|---:|")
        for (variant, score) in variant_scores
            println(io, @sprintf("| `%s` | %.3f |", String(variant), score))
        end
        println(io)
        println(io, "## Strength Scan")
        println(io)
        println(io, @sprintf("Best scalar multiplier on the current P1 annihilation matrix: `%.4g`, with mean vector RMS `%.3f`.", strength.strength, strength.score))
        println(io, @sprintf("Closest scalar multiplier in the scan to the Table III/GI `eta(548)` mass is `%.4g`, giving first pole `%.3f GeV` and mean vector RMS `%.3f`.", eta_strength.strength, eta_strength.mass, eta_strength.score))
        println(io)
        println(io, "## Current Finding")
        println(io)
        println(io, "This does not look like a column-ordering, phase, one-line sign, coordinate/momentum `alpha_s`, Eq. (17) FD integral, or `ns` flavor-coherence bug. The exact FD momentum integral and coherent `ns` factor are now implemented in the paper P1/P2 path, but they still under-mix Table III. The remaining focused suspect is the paper-order basis entering the annihilation matrix: HO-diagonalized wavefunctions/diagonal masses and the phenomenological pseudoscalar prescription, rather than a local prefactor retune.")
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
