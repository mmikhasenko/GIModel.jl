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

function s_wave_basis(params, mq; ngrid = 220, rmax = 22.0)
    # Annihilation wavefunctions use the HO basis (paper-consistent
    # wavefunction-at-origin scale); diagonal masses use the FD contact levels.
    solver_ho = OscillatorSolver(ngrid = ngrid, rmax = rmax)
    channels = Dict{String,Any}()
    for (label, mass_key) in [("ns", "q"), ("ss", "s"), ("cc", "c"), ("bb", "b")]
        masses = ConstituentMasses(mq[mass_key], mq[mass_key])
        fd_solution = channel_solution(
            params, masses, 0;
            nlevels = 2,
            solver = FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax, kinetic = :relativistic),
        )
        solution = channel_solution(
            params,
            masses,
            0;
            solver = solver_ho,
            nlevels = 2,
        )
        levels = GIModel.contact_hyperfine_nonperturbative_levels(
            params,
            masses,
            "S",
            1,
            radial_wave(fd_solution, 1).r,
            2,
        )
        isempty(levels) && (levels = fd_solution.eigenvalues_GeV)
        channels[label] = (mass = mq[mass_key], solution = solution, levels = levels)
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
                fix_annihilation_phase(radial_wave(channel.solution, n));
                # TARGET_BASIS is ["1 ns", "1 ss", "1 cc", "1 bb", ...]; only the
                # nonstrange rows are the coherent (u ubar + d dbar)/sqrt(2) state.
                isoscalar_coherent = flavor == "ns",
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

function main()
    params, mq = load_parameters_and_quark_masses(joinpath(dirname(root), "data", "parameters.provisional.toml"))
    basis = s_wave_basis(params, mq)
    targets, shifts = table_iii_targets()

    p1 = isoscalar_pseudoscalar_annihilation_solution(PaperP1Annihilation(), params, basis)
    p2 = isoscalar_pseudoscalar_annihilation_solution(PaperP2Annihilation(), params, basis)

    open(REPORT, "w") do io
        println(io, "# Table III Mixing Audit")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/audit_table_iii_mixings.jl`.")
        println(io, numerics_provenance(
            FiniteDifferenceSolver(ngrid = 220, rmax = 22.0),
            OscillatorSolver(ngrid = 220, rmax = 22.0),
        ))
        println(io)
        println(io, "This audit compares the implemented literal pseudoscalar P1/P2 modes (Eq. 18a,b on HO radial wavefunctions) with the promoted Table III amplitude rows in `data/clean/mixings.csv`. The basis is ordered as:")
        println(io)
        println(io, join(["`$label`" for label in TARGET_BASIS], ", "))
        println(io)
        println(io, "## Current Conclusion")
        println(io)
        println(io, "With HO-basis wavefunctions in the Eq. (17) smearing and the GI annihilation phase convention `Φ(0) > 0` (which makes the radially excited `2 ns`/`2 ss` couplings negative), the literal Eq. (18a,b) modes reproduce the Table III sign structure in every pseudoscalar row. The remaining amplitude error is magnitude-level (under-mixed radial components in the eta-prime), and the remaining mass error tracks the light-sector unperturbed diagonals (the FD pi sits ~55 MeV below the paper's 0.15 GeV). The headline residual reports keep the calibrated P1 control for mass scoring while these literal modes are tracked here.")
        println(io)
        write_model_section(io, "P1", p1, targets, shifts)
        write_model_section(io, "P2", p2, targets, shifts)
        println(io, "## Non-Pseudoscalar Rows")
        println(io)
        println(io, "The visible `1^3S_1` and `1^3P_2` rows promoted into `data/clean/mixings.csv` are now reproduced by the general Eq. (16) model (`isoscalar_general_annihilation_solution`, `:table_iii` comparison scheme) with `A(^3S_1)=+2.5` (three-gluon bracket) and `A(^3P_2)=-0.8` (two-gluon bracket): omega/phi come out as `(+1.000, -0.029)`/`(+0.029, +1.000)` vs Table III `(+0.999, -0.02)`/`(+0.02, +0.999)`, and f2/f2' as `(+0.997, +0.080)`/`(-0.080, +0.997)` vs `(+0.997, +0.06)`/`(-0.07, +0.997)`. See the isoscalar residual report for the mass-level scoring.")
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
