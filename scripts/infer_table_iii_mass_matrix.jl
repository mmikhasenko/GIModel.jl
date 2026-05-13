#!/usr/bin/env julia
# Invert the visible light-pseudoscalar Table III amplitudes into an effective
# mass matrix and compare it with the current FD Eq. (16)-(18) construction.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using LinearAlgebra
using Printf
using Statistics

root = dirname(@__DIR__)
using GIModel

const LIGHT_BASIS = ["1 ns", "1 ss", "2 ns", "2 ss"]
const STATE_ORDER = ["eta(548)", "eta_prime(958)", "eta_r(?)", "eta_r_prime(?)"]
const TARGET_MASSES = Dict(
    "P1" => [0.520, 0.960, 1.440, 1.630],
    "P2" => [0.520, 0.960, 1.270, 1.550],
)
const REPORT = joinpath(root, "docs", "residual_reports", "table_iii_inverse_matrix.md")

function light_s_wave_basis(params, mq; ngrid = 220, rmax = 22.0)
    channels = Dict{String,Any}()
    for (label, mass_key) in [("ns", "q"), ("ss", "s")]
        masses = ConstituentMasses(mq[mass_key], mq[mass_key])
        ev, vecs, r = channel_solution(
            params,
            masses,
            0;
            nlevels = 2,
            ngrid = ngrid,
            rmax = rmax,
            kinetic = :relativistic,
        )
        sol = GIModel.ChannelRadialSolution(ev, vecs, r)
        levels = GIModel.contact_hyperfine_nonperturbative_levels(
            params,
            masses,
            "S",
            1,
            r,
            2,
        )
        channels[label] = (mass = mq[mass_key], solution = sol, levels = levels)
    end

    basis = GIModel.PseudoscalarAnnihilationBasisInput[]
    for label in LIGHT_BASIS
        n = parse(Int, first(split(label)))
        flavor = last(split(label))
        channel = channels[flavor]
        push!(
            basis,
            pseudoscalar_annihilation_basis_input(
                label,
                channel.mass,
                channel.levels[n],
                RadialWaveOnUniformMesh(channel.solution, n),
            ),
        )
    end
    return basis
end

function table_vectors(model)
    rows = collect(CSV.File(joinpath(root, "data", "clean", "mixings.csv")))
    vectors = Dict{String,Vector{Float64}}(
        state => zeros(Float64, length(LIGHT_BASIS)) for state in STATE_ORDER
    )
    for row in rows
        ismissing(row.model) && continue
        String(row.model) == model || continue
        state = String(row.state_name)
        haskey(vectors, state) || continue
        idx = findfirst(==(String(row.basis)), LIGHT_BASIS)
        isnothing(idx) && continue
        vectors[state][idx] = Float64(row.amplitude)
    end
    V = hcat([vectors[state] ./ norm(vectors[state]) for state in STATE_ORDER]...)
    return V
end

function implied_matrix(model)
    V = table_vectors(model)
    masses = TARGET_MASSES[model]
    raw = V * Diagonal(masses) / V
    return Symmetric(0.5 .* (raw .+ raw')), V
end

function paper_matrix(model, params, basis)
    solution = if model == "P1"
        isoscalar_pseudoscalar_annihilation_solution(PaperP1Annihilation(), params, basis)
    else
        isoscalar_pseudoscalar_annihilation_solution(PaperP2Annihilation(), params, basis)
    end
    return solution.block.matrix
end

function matrix_table(io, matrix; scale = 1000.0)
    println(io, "| row \\ col | ", join(["`$b`" for b in LIGHT_BASIS], " | "), " |")
    println(io, "|---", repeat("|---:", length(LIGHT_BASIS)), "|")
    for (i, label) in enumerate(LIGHT_BASIS)
        vals = [@sprintf("%.1f", scale * matrix[i, j]) for j in axes(matrix, 2)]
        println(io, "| `$label` | ", join(vals, " | "), " |")
    end
end

function offdiag_vector(matrix)
    out = Float64[]
    for j in axes(matrix, 2), i in 1:(j - 1)
        push!(out, matrix[i, j])
    end
    return out
end

function scalar_fit(required, current)
    r = offdiag_vector(required)
    c = offdiag_vector(current)
    denom = dot(c, c)
    denom <= eps(Float64) && return NaN
    return dot(r, c) / denom
end

function phase_optimized_fit(required, current)
    best = (rms = Inf, fit = NaN, phases = ones(Float64, size(required, 1)))
    for tuple in Iterators.product(fill((-1.0, 1.0), size(required, 1))...)
        phases = collect(Float64, tuple)
        transformed = Diagonal(phases) * current * Diagonal(phases)
        fit = scalar_fit(required, transformed)
        fitted = fit .* transformed
        rms = sqrt(mean(abs2, offdiag_vector(required - fitted))) * 1000
        if rms < best.rms
            best = (rms = rms, fit = fit, phases = phases)
        end
    end
    return best
end

function write_model(io, model, params, basis)
    implied, V = implied_matrix(model)
    current = paper_matrix(model, params, basis)
    diag = Diagonal([state.diagonal_GeV for state in basis])
    required_annihilation = Matrix(implied - diag)
    current_annihilation = Matrix(current - diag)
    fit = scalar_fit(required_annihilation, current_annihilation)
    fitted = fit .* current_annihilation
    offdiag_rms = sqrt(mean(abs2, offdiag_vector(required_annihilation - fitted))) * 1000
    phase_fit = phase_optimized_fit(required_annihilation, current_annihilation)

    println(io, "## `$model`")
    println(io)
    println(io, "Target masses used: ", join([@sprintf("%.3f", m) for m in TARGET_MASSES[model]], ", "), " GeV.")
    println(io, @sprintf("Table-vector Gram max off-diagonal: `%.3f`.", maximum(abs.(V' * V - I))))
    println(io, @sprintf("Best scalar multiplier from current off-diagonal annihilation to the Table-implied off-diagonal block: `%.3g`; residual off-diagonal RMS: `%.1f MeV`.", fit, offdiag_rms))
    println(io, @sprintf("After exhaustive basis-phase flips: best scalar multiplier `%.3g`; residual off-diagonal RMS `%.1f MeV`; phases `%s`.", phase_fit.fit, phase_fit.rms, join([@sprintf("%+.0f", p) for p in phase_fit.phases], ", ")))
    println(io)
    println(io, "### Current FD Diagonal")
    println(io)
    println(io, "| basis | diagonal GeV |")
    println(io, "|---|---:|")
    for state in basis
        println(io, @sprintf("| `%s` | %.3f |", state.label, state.diagonal_GeV))
    end
    println(io)
    println(io, "### Table-Implied Annihilation Block (MeV)")
    println(io)
    matrix_table(io, required_annihilation)
    println(io)
    println(io, "### Current Eq. (16)-(18) Annihilation Block (MeV)")
    println(io)
    matrix_table(io, current_annihilation)
    println(io)
end

function main()
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    basis = light_s_wave_basis(params, mq)
    open(REPORT, "w") do io
        println(io, "# Table III Inverse Mass-Matrix Probe")
        println(io)
        println(io, "Generated by `julia scripts/infer_table_iii_mass_matrix.jl`.")
        println(io)
        println(io, "This report asks a scientist's inverse question: if the visible light pseudoscalar Table III amplitudes are treated as eigenvectors, what effective mass matrix do they imply, and does it resemble the current FD Eq. (16)-(18) matrix?")
        println(io)
        println(io, "The inversion uses only the light block `1 ns`, `1 ss`, `2 ns`, `2 ss`. Charm and bottom components are omitted because they are tiny in the visible light pseudoscalar rows; the remaining vectors are renormalized. The P2 reconstruction is only a diagnostic, since the paper notes that P2 poles are mass-dependent and need not be orthogonal.")
        println(io)
        write_model(io, "P1", params, basis)
        write_model(io, "P2", params, basis)
        println(io, "## Interpretation")
        println(io)
        println(io, "If the Table-implied block differs from the FD Eq. (16)-(18) block by more than a mostly uniform scale, the mismatch is unlikely to be a local sign or normalization bug. It points instead to the basis entering the annihilation problem: paper-order HO eigenvectors/diagonal masses, omitted coupled-channel/chiral pseudoscalar physics, or the phenomenological construction behind P1/P2.")
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
