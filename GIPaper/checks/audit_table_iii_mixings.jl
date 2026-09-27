#!/usr/bin/env julia
# Compare the implemented pseudoscalar P1/P2 annihilation modes with the
# promoted Table III amplitude rows.

mkpath(joinpath(dirname(@__DIR__), "reports"))

using Pkg
Pkg.activate(dirname(@__DIR__))

using CSV
using LinearAlgebra
using Printf

root = dirname(@__DIR__)
using GIModel

const TARGET_BASIS = ["1 ns", "1 ss", "1 cc", "1 bb", "2 ns", "2 ss", "2 cc"]
const REPORT = joinpath(root, "reports", "table_iii_mixing_audit.md")

# Paper mass targets for the literal modes. P1 masses are the digitized Fig. 5
# isoscalar pseudoscalar labels (and Fig. 6 for the eta_c). P2 masses follow
# the Table III Δm column relative to the isovector pi(0.15)/pi'(1.30) for the
# eta/eta', and the Sec. VA text (1.27 and 1.55 GeV) for the third and fourth
# poles.
const MASS_TARGETS_GEV = Dict(
    "P1" => Dict(
        "eta(548)" => 0.52,
        "eta_prime(958)" => 0.96,
        "eta_r(?)" => 1.44,
        "eta_r_prime(?)" => 1.63,
        "eta_c(2980)" => 2.97,
    ),
    "P2" => Dict(
        "eta(548)" => 0.49,
        "eta_prime(958)" => 0.93,
        "eta_r(?)" => 1.27,
        "eta_r_prime(?)" => 1.55,
        "eta_c(2980)" => 2.97,
    ),
)

function s_wave_basis(params, mq)
    levels = [BasisState(n, "S", 1, 0) for n in 1:2]
    # Use the same native full fixed-channel HO masses and waves as the decay
    # drivers; do not mix FD diagonal masses with spin-independent HO waves.
    channels = Dict(flavor => fixed_spectrum(params, Meson(mq, f, f);
        levels, solver=OscillatorSolver()) for
        (flavor, f) in (("ns", :q), ("ss", :s), ("cc", :c), ("bb", :b)))
    return [annihilation_basis_input(channels[last(split(label))],
        levels[parse(Int, first(split(label)))]) for label in TARGET_BASIS]
end

function table_iii_targets()
    rows = collect(CSV.File(joinpath(root, "data", "transcription_mixings.csv")))
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
    println(io, "| state | pole / column | mass GeV | target GeV | mass Δ MeV | Table Δ MeV | vector RMS | max | largest model components |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---|")
    rms_values = Float64[]
    mass_residuals = Float64[]
    mass_targets = get(MASS_TARGETS_GEV, model, Dict{String,Float64}())
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
        mass_target = get(mass_targets, key[1], missing)
        mass_target_text = ismissing(mass_target) ? "n/a" : @sprintf("%.3f", mass_target)
        mass_delta_text = "n/a"
        if !ismissing(mass_target)
            delta = 1000 * (solution.masses[col] - mass_target)
            push!(mass_residuals, abs(delta))
            mass_delta_text = @sprintf("%+.0f", delta)
        end
        println(
            io,
            @sprintf(
                "| `%s` | %d | %.3f | %s | %s | %s | %.3f | %.3f | %s |",
                key[1],
                col,
                solution.masses[col],
                mass_target_text,
                mass_delta_text,
                shift_text,
                rms,
                maxerr,
                largest_text,
            ),
        )
    end
    println(io)
    println(io, "### Complete signed components")
    println(io)
    println(io, "| state | " * join(["$b model | $b paper" for b in TARGET_BASIS], " | ") * " |")
    println(io, "|---|" * repeat("---:|", 2length(TARGET_BASIS)))
    for (i, key) in enumerate(target_keys)
        pred = signs[i] .* solution.vectors[:, cols[i]]
        @assert isapprox(sum(abs2, pred), 1; atol=1e-10)
        cells = [@sprintf("%+.6g | %+.6g", pred[j], tvecs[i][j]) for j in eachindex(pred)]
        println(io, "| `$(key[1])` | " * join(cells, " | ") * " |")
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
    if !isempty(mass_residuals)
        println(
            io,
            @sprintf(
                "Mean abs mass residual for `%s`: `%.0f MeV` against the paper-model targets.",
                model,
                sum(mass_residuals) / length(mass_residuals),
            ),
        )
    end
    println(io)
end

function write_general_sections(io, params, mq)
    rows = collect(CSV.File(joinpath(root, "data", "transcription_mixings.csv")))
    for (group, L, J) in [("1 3S1", "S", 1), ("1 3P2", "P", 2)]
        level = BasisState(1, L, 3, J)
        basis = [annihilation_basis_input(fixed_spectrum(params, Meson(mq, f, f);
            levels=[level], solver=OscillatorSolver()), level) for f in (:q, :s, :c, :b)]
        solution = isoscalar_general_annihilation_solution(params, basis;
            amplitude_A=L == "S" ? params.annihilation.s1_A : params.annihilation.a_3p2,
            L=L == "S" ? 0 : 1, multiplicity=3, J=J)
        labels = TARGET_BASIS[1:4]
        selected = filter(r -> r.state_group == group, rows)
        names = unique(String(r.state_name) for r in selected)
        targets = [[only(Float64(r.amplitude) for r in selected if r.state_name == name && r.basis == b)
            for b in labels] for name in names]
        cols, signs = aligned_columns(solution.vectors, targets)
        println(io, "### $group — Eq. (16), four flavors")
        println(io)
        println(io, "Radial n=2 entries are not part of these original Table III blocks; they are not zero predictions.")
        println(io)
        println(io, "| state | mass GeV | " * join(["$b model | $b paper" for b in labels], " | ") * " |")
        println(io, "|---|---:|" * repeat("---:|", 2length(labels)))
        for i in eachindex(names)
            pred = signs[i] .* solution.vectors[:, cols[i]]
            @assert isapprox(sum(abs2, pred), 1; atol=1e-10)
            cells = [@sprintf("%+.6g | %+.6g", pred[j], targets[i][j]) for j in eachindex(pred)]
            println(io, "| `$(names[i])` | $(@sprintf("%.6f", solution.masses[cols[i]])) | " * join(cells, " | ") * " |")
        end
        group == "1 3S1" && println(io)
    end
end

function main()
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    basis = s_wave_basis(params, mq)
    targets, shifts = table_iii_targets()

    p1 = isoscalar_pseudoscalar_annihilation_solution(PaperP1Annihilation(), params, basis)
    p2 = isoscalar_pseudoscalar_annihilation_solution(PaperP2Annihilation(), params, basis)

    open(REPORT, "w") do io
        println(io, "# Table III Mixing Audit")
        println(io)
        println(io, "Generated by `julia GIPaper/checks/audit_table_iii_mixings.jl`.")
        println(io, numerics_provenance(OscillatorSolver()))
        println(io)
        println(io, "This audit compares the implemented literal pseudoscalar P1/P2 modes (Eq. 18a,b on HO radial wavefunctions) with the promoted Table III amplitude rows in `data/transcription_mixings.csv`. The basis is ordered as:")
        println(io)
        println(io, join(["`$label`" for label in TARGET_BASIS], ", "))
        println(io)
        println(io, "## Current Conclusion")
        println(io)
        println(io, "All diagonal masses and Eq. (17) overlaps use the same native full fixed-channel oscillator wavefunctions. P1/P2 use seven flavor/radial entries; vector and tensor channels use the original four-flavor ground-state blocks. Every signed component is retained below, including small charm and bottom admixtures. State assignment and global signs maximize overlap with the reference vector; no component is adjusted to its paper value.")
        println(io)
        write_model_section(io, "P1", p1, targets, shifts)
        write_model_section(io, "P2", p2, targets, shifts)
        println(io, "## Non-Pseudoscalar Rows")
        println(io)
        write_general_sections(io, params, mq)
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
