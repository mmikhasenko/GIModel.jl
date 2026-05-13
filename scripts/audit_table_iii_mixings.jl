#!/usr/bin/env julia
# Compare the implemented pseudoscalar P1/P2 annihilation modes with the
# promoted Table III amplitude rows.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using LinearAlgebra
using Printf

root = dirname(@__DIR__)
using GIModel

const TARGET_BASIS = ["1 ns", "1 ss", "1 cc", "1 bb", "2 ns", "2 ss", "2 cc"]
const REPORT = joinpath(root, "docs", "residual_reports", "table_iii_mixing_audit.md")

function s_wave_basis(params, mq; ngrid = 220, rmax = 22.0)
    channels = Dict{String,Any}()
    for (label, mass_key) in [("ns", "q"), ("ss", "s"), ("cc", "c"), ("bb", "b")]
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
    for label in TARGET_BASIS
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

function table_iii_targets()
    rows = collect(CSV.File(joinpath(root, "data", "clean", "mixings.csv")))
    targets = Dict{Tuple{String,String},Dict{String,Float64}}()
    shifts = Dict{Tuple{String,String},Union{Missing,Float64}}()
    for row in rows
        model = ismissing(row.model) ? "" : String(row.model)
        model in ("P1", "P2") || continue
        key = (String(row.state_name), model)
        dict = get!(targets, key, Dict{String,Float64}())
        dict[String(row.basis)] = Float64(row.amplitude)
        if !ismissing(row.mass_shift_MeV)
            shifts[key] = Float64(row.mass_shift_MeV)
        elseif !haskey(shifts, key)
            shifts[key] = missing
        end
    end
    return targets, shifts
end

function aligned_columns(vectors::AbstractMatrix, target_order)
    used = Set{Int}()
    cols = Int[]
    signs = Float64[]
    for target in target_order
        best = 0
        bestdot = -Inf
        for col in axes(vectors, 2)
            col in used && continue
            d = dot(vectors[:, col], target)
            if abs(d) > bestdot
                bestdot = abs(d)
                best = col
            end
        end
        push!(used, best)
        push!(cols, best)
        push!(signs, dot(vectors[:, best], target) < 0 ? -1.0 : 1.0)
    end
    return cols, signs
end

function target_vector(dict)
    return [get(dict, label, 0.0) for label in TARGET_BASIS]
end

function model_symbol(::PaperP1Annihilation)
    return "P1"
end

function model_symbol(::PaperP2Annihilation)
    return "P2"
end

function write_model_section(io, title, solution, targets, shifts)
    model = title
    target_keys = sort([key for key in keys(targets) if key[2] == model]; by = key -> key[1])
    tvecs = [target_vector(targets[key]) for key in target_keys]
    cols, signs = aligned_columns(solution.vectors, tvecs)

    println(io, "## `$model`")
    println(io)
    println(io, "| state | pole / column | mass GeV | Table Δ MeV | vector RMS | max | largest model components |")
    println(io, "|---|---:|---:|---:|---:|---:|---|")
    rms_values = Float64[]
    for (i, key) in enumerate(target_keys)
        col = cols[i]
        pred = signs[i] .* solution.vectors[:, col]
        target = tvecs[i]
        diffs = pred .- target
        rms = sqrt(sum(abs2, diffs) / length(diffs))
        push!(rms_values, rms)
        maxerr = maximum(abs.(diffs))
        largest = sort(
            collect(zip(TARGET_BASIS, pred));
            by = item -> -abs(item[2]),
        )[1:3]
        largest_text = join([@sprintf("`%s`=%+.3f", label, amp) for (label, amp) in largest], ", ")
        shift = get(shifts, key, missing)
        shift_text = ismissing(shift) ? "n/a" : @sprintf("%.0f", shift)
        println(
            io,
            @sprintf(
                "| `%s` | %d | %.3f | %s | %.3f | %.3f | %s |",
                key[1],
                col,
                solution.masses[col],
                shift_text,
                rms,
                maxerr,
                largest_text,
            ),
        )
    end
    println(io)
    println(
        io,
        @sprintf(
            "Mean vector RMS for `%s`: `%.3f`; max vector RMS: `%.3f`.",
            model,
            sum(rms_values) / length(rms_values),
            maximum(rms_values),
        ),
    )
    println(io)
end

function main()
    params, mq = load_parameters_and_quark_masses(joinpath(root, "data", "parameters.provisional.toml"))
    basis = s_wave_basis(params, mq)
    targets, shifts = table_iii_targets()

    p1 = isoscalar_pseudoscalar_annihilation_solution(PaperP1Annihilation(), params, basis)
    p2 = isoscalar_pseudoscalar_annihilation_solution(PaperP2Annihilation(), params, basis)

    open(REPORT, "w") do io
        println(io, "# Table III Mixing Audit")
        println(io)
        println(io, "Generated by `julia scripts/audit_table_iii_mixings.jl`.")
        println(io)
        println(io, "This audit compares the implemented literal FD pseudoscalar P1/P2 modes with the promoted Table III amplitude rows in `data/clean/mixings.csv`. The basis is ordered as:")
        println(io)
        println(io, join(["`$label`" for label in TARGET_BASIS], ", "))
        println(io)
        println(io, "## Current Conclusion")
        println(io)
        println(io, "The current literal FD Eq. (18a,b) modes do **not** reproduce the Table III amplitudes. They now use the paper-facing momentum `alpha_s(Q^2)`, a direct FD S-wave Eq. (17) momentum integral, and coherent `ns` flavor normalization, but still generate small annihilation-driven admixtures while Table III's light pseudoscalars have large `ns`/`ss` and radial components. This confirms that the calibrated mass-control path is not yet a Table-III eigenvector reproduction, and that the next technical step is paper-order basis/matrix-element fidelity rather than retuning comparison residuals.")
        println(io)
        write_model_section(io, "P1", p1, targets, shifts)
        write_model_section(io, "P2", p2, targets, shifts)
        println(io, "## Non-Pseudoscalar Rows")
        println(io)
        println(io, "The visible `1^3S_1` and `1^3P_2` rows are now promoted into `data/clean/mixings.csv` from the page-011 markdown. They use `A(^3S_1)=+2.5` and `A(^3P_2)=-0.8`, but this audit does not score them yet because the codebase does not yet have a general non-pseudoscalar Eq. (16) annihilation model.")
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
