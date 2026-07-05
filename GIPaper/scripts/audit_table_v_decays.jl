#!/usr/bin/env julia
# Table V strong-decay audit: compute the light 1S+1P harmonic-oscillator
# amplitudes from the two-parameter (A, S0) model and compare against the
# paper's numeric column. See docs/observable_ledger.md for conventions.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using Printf
using GIModel

root = dirname(@__DIR__)
const TABLE = joinpath(
    root, "data", "raw", "digitized_tables", "table_v_strong_decays",
    "table_v_light_1s_1p.provisional.csv",
)
const REPORT = joinpath(root, "docs", "residual_reports", "table_v_light_decays.md")

# Parent masses (GeV): 1984-era experimental values for established states,
# GI model masses (Figs. 3-5) for states the paper itself predicted
# (H, H', delta2, eps, eps', kappa as parents).
const PARENT_MASS = Dict(
    "rho" => 0.769,
    "phi" => 1.0195,
    "Kstar" => 0.8921,
    "A2" => 1.318,
    "f" => 1.273,
    "fprime" => 1.525,
    "Kstar2" => 1.434,
    "A1" => 1.27,
    "D1285" => 1.283,
    "E1420" => 1.42,
    "Q1" => 1.28,
    "Q2" => 1.40,
    "B" => 1.231,
    "H" => 1.22,        # GI model 1^1P_1 isoscalar
    "Hprime" => 1.47,   # GI model 1^1P_1 isoscalar (ss)
    "delta2" => 1.09,   # GI model 1^3P_0 isovector
    "eps" => 1.09,      # GI model 1^3P_0 isoscalar
    "epsprime" => 1.36, # GI model 1^3P_0 isoscalar (ss)
    "kappa" => 1.24,    # GI model 1^3P_0 strange
)

# Daughter masses (GeV). Quasi-two-body subchannels ((eta pi)_delta2,
# (pi pi)_eps, (K pi)_kappa) carry an effective mass; values below are the
# 1984 resonance masses with the open conventions flagged in the table notes.
const DAUGHTER_MASS = Dict(
    "pi" => 0.138,
    "K" => 0.4957,
    "eta" => 0.5488,
    "etaprime" => 0.9575,
    "rho" => 0.769,
    "omega" => 0.7826,
    "Kstar" => 0.8921,
    "delta2" => 0.98,
    "eps" => 1.00,
    "kappa" => 0.98,
)

# Rows whose paper numerics depend on machinery beyond the two-parameter
# equal-mass convention; they are computed but excluded from the headline
# score, each with its reason.
function exclusion_reason(row)
    parent = String(row.parent)
    label = String(row.decay_label)
    notes = ismissing(row.notes) ? "" : String(row.notes)
    if parent in ("Q1", "Q2")
        # The numeric column uses the model's strange-axial (K1) mixing angle
        # on top of the unmixed formula coefficients.
        return "model K1 mixing angle"
    end
    occursin("kinematics convention open", notes) && return "quasi-two-body daughter mass"
    occursin("image audit pending", notes) && return "image audit pending"
    label == "D -> [(K pi)_K* Kbar]_S" && return "below nominal K* K threshold (paper integrates lineshape)"
    label == "D -> [(K pi)_K* Kbar]_D" && return "below nominal K* K threshold (paper integrates lineshape)"
    label == "D -> (eta pi)_delta2 pi" && return "image audit pending (coefficient suspect)"
    return nothing
end

function main()
    rows = collect(CSV.File(TABLE))

    fit_kinematics = Dict{String,Float64}()
    for row in rows
        label = String(row.decay_label)
        if label == "rho -> pi pi" || label == "B -> [omega pi]_S"
            M = PARENT_MASS[String(row.parent)]
            m1 = DAUGHTER_MASS[String(row.daughter1)]
            m2 = DAUGHTER_MASS[String(row.daughter2)]
            fit_kinematics[label] = decay_momentum(M, m1, m2)
        end
    end
    model = calibrate_strong_decay_model(
        fit_kinematics["rho -> pi pi"],
        fit_kinematics["B -> [omega pi]_S"],
    )

    computed_rows = NamedTuple[]
    for row in rows
        class = Symbol(String(row.amplitude_class))
        label = String(row.decay_label)
        paper_text = strip(String(row.paper_MeV12))
        M = PARENT_MASS[String(row.parent)]
        m1 = DAUGHTER_MASS[String(row.daughter1)]
        m2 = DAUGHTER_MASS[String(row.daughter2)]
        q = decay_momentum(M, m1, m2)
        class_map = Dict(
            :A => :A, :A0 => :A0, :Aprime => :Aprime, :Adoubleprime => :Adoubleprime,
            :S => :S, :D => :D, :P => :P,
        )
        if class === :mixing_only
            push!(computed_rows, (
                label = label, section = String(row.section), class = "mixing",
                q_GeV = q, computed = missing, paper = paper_text,
                ratio = missing, flagged = false,
                excluded = "Table III mixing only", small = false,
                notes = ismissing(row.notes) ? "" : String(row.notes),
            ))
            continue
        end
        amp = strong_decay_amplitude(
            model,
            Float64(row.coefficient_value),
            class_map[class],
            Int(row.qbar_power),
            q,
        )
        paper_value = tryparse(Float64, replace(paper_text, "+" => ""))
        below = paper_text == "below_threshold"
        ratio = (isnothing(paper_value) || paper_value == 0) ? missing : amp / paper_value
        excluded = exclusion_reason(row)
        # Small paper amplitudes are rounded to one digit; score them on the
        # absolute difference instead of the ratio.
        small = !isnothing(paper_value) && abs(paper_value) <= 1.0
        flagged = isnothing(excluded) && !below && !ismissing(ratio) &&
                  (small ? abs(amp - paper_value) > 0.3 : abs(ratio - 1) > 0.20)
        push!(computed_rows, (
            label = label, section = String(row.section), class = String(row.amplitude_class),
            q_GeV = q, computed = amp,
            paper = below ? "below threshold" : paper_text,
            ratio = ratio, flagged = flagged, excluded = excluded, small = small,
            notes = ismissing(row.notes) ? "" : String(row.notes),
        ))
    end

    scored = [r for r in computed_rows if !ismissing(r.ratio) && isnothing(r.excluded) && !r.small]
    deviations = [abs(r.ratio - 1) for r in scored]
    n_flagged = count(r -> r.flagged, computed_rows)
    n_excluded = count(r -> hasproperty(r, :excluded) && !isnothing(r.excluded), computed_rows)

    open(REPORT, "w") do io
        println(io, "# Table V Light 1S+1P Strong-Decay Audit")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/audit_table_v_decays.jl`.")
        println(io)
        println(io, "Two-parameter Table IV/V model calibrated on the paper's fit rows:")
        println(io, @sprintf("- `A = %.3f` from `rho -> pi pi = +12.4 MeV^(1/2)` (q = %.0f MeV)", model.A, 1000fit_kinematics["rho -> pi pi"]))
        println(io, @sprintf("- `S0 = %.3f` from `B -> [omega pi]_S = -11 MeV^(1/2)` (q = %.0f MeV)", model.S0, 1000fit_kinematics["B -> [omega pi]_S"]))
        println(io, "- `beta = 0.40 GeV`; amplitude conventions in `docs/observable_ledger.md`.")
        println(io)
        println(io, @sprintf(
            "Headline rows (clean two-parameter convention, |paper| > 1): %d; median abs deviation %.0f%%; flagged (>20%%): %d. Excluded rows (computed but not scored): %d.",
            length(scored), 100 * sort(deviations)[max(1, (length(deviations) + 1) ÷ 2)], n_flagged, n_excluded,
        ))
        println(io)
        println(io, "| decay | class | q MeV | computed | paper | ratio | status |")
        println(io, "|---|---|---:|---:|---:|---:|---|")
        for r in computed_rows
            computed_text = ismissing(r.computed) ? "n/a" : @sprintf("%+.2f", r.computed)
            ratio_text = ismissing(r.ratio) ? "" : @sprintf("%.2f", r.ratio)
            status = if !isnothing(r.excluded)
                "excluded: $(r.excluded)"
            elseif r.flagged
                "**flag**"
            elseif r.small
                "small-amplitude (abs)"
            else
                "ok"
            end
            println(io, @sprintf(
                "| `%s` | %s | %.0f | %s | %s | %s | %s |",
                r.label, r.class, 1000r.q_GeV, computed_text, r.paper, ratio_text, status,
            ))
        end
        println(io)
        println(io, "## Open Conventions")
        println(io)
        println(io, "- `Q1`/`Q2` numeric amplitudes fold in the model's strange-axial (`K1`) mixing angle on top of the unmixed formula coefficients; reproducing them needs the `1^3P_1`/`1^1P_1` strange mixing block (machinery exists in the spectrum layer).")
        println(io, "- Quasi-two-body subchannel daughters (`(pi pi)_eps`, `(K pi)_kappa`, `(eta pi)_delta2`) and sub-threshold modes (`D -> K* Kbar`) depend on lineshape conventions the paper does not state; nominal-mass kinematics cannot reproduce them.")
        println(io, "- Strange-parent normalization: `K*2 -> K pi` computes ~sqrt(3) high while `K*(892) -> K pi` is exact, pointing at unequal-mass recoil factors in the Appendix B amplitudes not yet encoded.")
    end
    println("wrote ", REPORT)
    @printf("A=%.3f S0=%.3f scored=%d flagged=%d excluded=%d\n", model.A, model.S0, length(scored), n_flagged, n_excluded)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
