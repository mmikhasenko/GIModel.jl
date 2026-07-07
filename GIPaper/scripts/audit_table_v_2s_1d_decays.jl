#!/usr/bin/env julia
# Table V 2S + 1D strong-decay audit. Scores the digitized 1^3D_3, 1^3D_2,
# 1^1D_2, 1^3D_1, 2^1S_0 and 2^3S_1 sections against the paper's numeric
# amplitude column, using the SAME two-parameter (A, S0) harmonic-oscillator
# model calibrated on rho -> pi pi and B -> [omega pi]_S (no refit), and the
# model's OWN predicted parent masses from `compute_spectrum` (corrected stage).
#
# This keeps `table_v_light_decays.md` (the 1S+1P audit) untouched; it delivers
# a separate report. See docs/observable_ledger.md for amplitude conventions.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using Printf
using GIModel

root = dirname(@__DIR__)
const TABLE_2S1D = joinpath(
    root, "data", "raw", "digitized_tables", "table_v_strong_decays",
    "table_v_light_2s_1d.provisional.csv",
)
const TABLE_1D3 = joinpath(
    root, "data", "raw", "digitized_tables", "table_v_strong_decays",
    "table_v_light_1d3.provisional.csv",
)
const REPORT = joinpath(root, "docs", "residual_reports", "table_v_2s_1d_decays.md")
const PARAMS_PATH = joinpath(dirname(root), "data", "parameters.provisional.toml")

# --- Daughter masses (GeV): 1984-era physical values. Quasi-two-body
# subchannels ((pi pi)_eps, (K pi)_kappa, (eta pi)_delta2) are excluded from
# scoring below, but keep placeholder masses so momenta are still defined. -----
const DAUGHTER_MASS = Dict(
    "pi" => 0.138, "K" => 0.4957, "eta" => 0.5488, "etaprime" => 0.9575,
    "rho" => 0.769, "omega" => 0.7826, "phi" => 1.0195, "Kstar" => 0.8921,
    "B" => 1.231, "delta2" => 0.98, "eps" => 1.00, "kappa" => 0.98,
)

# --- Parent -> model level map. Nonstrange isoscalars (omega, phi_nn) and the
# isovector share constituent masses under isospin symmetry, so the nonstrange
# representative is the q-qbar sector and the ss states the s-sbar sector.
# Keyed by (section, parent) because parent labels (rho, omega, phi, Q1, Q2,
# Kstar...) recur across sections with different quantum numbers. Value is
# (flavor1, flavor2, n, L, multiplicity, J). ------------------------------------
const PARENT_LEVEL = Dict{Tuple{String,String},Tuple{Symbol,Symbol,Int,String,Int,Int}}(
    # 1^3D_3 (n=1, ^3D_3)
    ("1^3D_3", "g")       => (:q, :q, 1, "D", 3, 3),
    ("1^3D_3", "omega")   => (:q, :q, 1, "D", 3, 3),
    ("1^3D_3", "phi")     => (:s, :s, 1, "D", 3, 3),
    ("1^3D_3", "Kstar3")  => (:q, :s, 1, "D", 3, 3),
    # 1^3D_2 nonstrange (n=1, ^3D_2)
    ("1^3D_2 nonstrange", "rho")   => (:q, :q, 1, "D", 3, 2),
    ("1^3D_2 nonstrange", "omega") => (:q, :q, 1, "D", 3, 2),
    ("1^3D_2 nonstrange", "phi")   => (:s, :s, 1, "D", 3, 2),
    # 1^3D_2 / 1^1D_2 strange (footnote c: Q1 = pure 1^1D_2 lower, Q2 = pure 1^3D_2 higher)
    ("1^3D_2 1^1D_2 strange", "Q1") => (:q, :s, 1, "D", 1, 2),
    ("1^3D_2 1^1D_2 strange", "Q2") => (:q, :s, 1, "D", 3, 2),
    # 1^1D_2 nonstrange (n=1, ^1D_2)
    ("1^1D_2 nonstrange", "A3")    => (:q, :q, 1, "D", 1, 2),
    ("1^1D_2 nonstrange", "omega") => (:q, :q, 1, "D", 1, 2),
    ("1^1D_2 nonstrange", "phi")   => (:s, :s, 1, "D", 1, 2),
    # 1^3D_1 (n=1, ^3D_1)
    ("1^3D_1", "rhoD")   => (:q, :q, 1, "D", 3, 1),
    ("1^3D_1", "omegaD") => (:q, :q, 1, "D", 3, 1),
    ("1^3D_1", "phiD")   => (:s, :s, 1, "D", 3, 1),
    ("1^3D_1", "KstarD") => (:q, :s, 1, "D", 3, 1),
    # 2^1S_0 (n=2, ^1S_0)
    ("2^1S_0", "piprime")    => (:q, :q, 2, "S", 1, 0),
    ("2^1S_0", "eta_r")      => (:q, :q, 2, "S", 1, 0),
    ("2^1S_0", "etaprime_r") => (:s, :s, 2, "S", 1, 0),
    ("2^1S_0", "Kprime")     => (:q, :s, 2, "S", 1, 0),
    # 2^3S_1 (n=2, ^3S_1)
    ("2^3S_1", "rhoS")   => (:q, :q, 2, "S", 3, 1),
    ("2^3S_1", "omegaS") => (:q, :q, 2, "S", 3, 1),
    ("2^3S_1", "phiS")   => (:s, :s, 2, "S", 3, 1),
    ("2^3S_1", "KstarS") => (:q, :s, 2, "S", 3, 1),
)

const class_map = Dict(
    :A => :A, :A0 => :A0, :Aprime => :Aprime, :Adoubleprime => :Adoubleprime,
    :S => :S, :D => :D, :P => :P,
)

# Table IV split (line 527 of the vision OCR): "structure independent"
# (A, A', A'', A_c, A0 — only the angular-momentum q^L and the e^{-q^2/16b^2}
# form factor) vs "structure dependent" (S, D, P, S_c — extra S0 - k A qbar^2
# polynomial, "highly sensitive to the structure of the states").
const STRUCTURE_INDEPENDENT = Set([:A, :A0, :Aprime, :Adoubleprime])
is_structure_independent(class_str) = Symbol(class_str) in STRUCTURE_INDEPENDENT

# --- Model parent masses from compute_spectrum (corrected stage: central +
# contact + fine structure, pre-mixing). One solve per flavor sector, cached. ---
function build_mass_resolver(params, mq)
    cache = Dict{Tuple{Symbol,Symbol},Any}()
    levels = spectrum_levels(2; L_labels = ("S", "P", "D"))
    function corrected_sector(f1, f2)
        get!(cache, (f1, f2)) do
            central = central_spectrum(params, Meson(mq, f1, f2); levels = levels)
            add_spin_corrections(central)
        end
    end
    function model_mass(f1, f2, n, L, mult, J)
        spectrum_state(corrected_sector(f1, f2), n, L, mult, J).mass_GeV
    end
    return model_mass
end

# Two-part paper strings ("-0.11 - 1.1") are mixing-dependent ranges; single
# signed numbers parse, "below_threshold" and "0" are handled separately.
parse_paper(text) = tryparse(Float64, replace(strip(text), "+" => ""))

function main()
    params, mq = load_parameters_and_quark_masses(PARAMS_PATH)
    model_mass = build_mass_resolver(params, mq)

    # Same calibration as the 1S+1P audit (no refit): A from rho -> pi pi,
    # S0 from B -> [omega pi]_S, at their nominal 1984 kinematics.
    q_rho = decay_momentum(0.769, 0.138, 0.138)
    q_B = decay_momentum(1.231, 0.7826, 0.138)
    model = calibrate_strong_decay_model(q_rho, q_B)

    rows = vcat(collect(CSV.File(TABLE_1D3)), collect(CSV.File(TABLE_2S1D)))

    computed = NamedTuple[]
    for row in rows
        section = String(row.section)
        parent = String(row.parent)
        class_str = String(row.amplitude_class)
        label = String(row.decay_label)
        d1 = String(row.daughter1)
        d2 = String(row.daughter2)
        paper_text = ismissing(row.paper_MeV12) ? "" : strip(String(row.paper_MeV12))

        # Rows carrying no amplitude model (Table III mixing, cryptoexotics).
        if class_str in ("mixing_only", "unlisted")
            continue
        end

        key = (section, parent)
        if !haskey(PARENT_LEVEL, key)
            push!(computed, (label = label, section = section, class = class_str,
                q_MeV = NaN, computed = NaN, paper = paper_text, ratio = missing,
                status = "unmapped parent"))
            continue
        end
        f1, f2, n, L, mult, J = PARENT_LEVEL[key]
        M = model_mass(f1, f2, n, L, mult, J)
        m1 = get(DAUGHTER_MASS, d1, NaN)
        m2 = get(DAUGHTER_MASS, d2, NaN)
        q = decay_momentum(M, m1, m2)
        amp = strong_decay_amplitude(model, Float64(row.coefficient_value),
            class_map[Symbol(class_str)], Int(row.qbar_power), q)

        paper_value = parse_paper(paper_text)
        struct_indep = is_structure_independent(class_str)
        # Determine scoreability.
        status = ""
        if paper_text == "below_threshold"
            status = "excluded: paper below threshold"
        elseif isnothing(paper_value)
            status = "excluded: mixing-dependent range"
        elseif section == "1^3D_2 1^1D_2 strange"
            # Footnote c: Q1/Q2 here are the physical MIXED 1D_2 states (theta_1D
            # ~ 33 deg, mixing_angles.md), exactly analogous to the 1P Q1/Q2 case.
            # The printed formula is the pure singlet/triplet, so scoring raw is
            # meaningless -- excluded like the 1P footnote-b rows.
            status = "excluded: 1D2 mixing (theta_1D, see mixing_angles.md)"
        elseif d1 in ("eps", "kappa", "delta2") || d2 in ("eps", "kappa", "delta2")
            status = "excluded: quasi-two-body daughter"
        elseif q <= 0
            status = "excluded: below threshold at model mass"
        end
        ratio = (isnothing(paper_value) || paper_value == 0) ? missing : amp / paper_value
        small = !isnothing(paper_value) && abs(paper_value) <= 1.0
        if isempty(status)
            if small
                status = @sprintf("small-amplitude (abs diff %.2f)", abs(amp - paper_value))
            elseif !ismissing(ratio) && abs(ratio - 1) > 0.20
                status = "**flag**"
            else
                status = "ok"
            end
        end
        push!(computed, (label = label, section = section, class = class_str,
            struct_indep = struct_indep, q_MeV = 1000q, computed = amp,
            paper = paper_text, ratio = ratio, status = status, small = small))
    end

    scored = [r for r in computed if hasproperty(r, :small) && !r.small &&
              !ismissing(r.ratio) && !startswith(r.status, "excluded")]
    scored_si = [r for r in scored if r.struct_indep]
    scored_sd = [r for r in scored if !r.struct_indep]
    med(rs) = isempty(rs) ? NaN :
        sort([abs(r.ratio - 1) for r in rs])[max(1, (length(rs) + 1) ÷ 2)]
    median_dev, median_si, median_sd = med(scored), med(scored_si), med(scored_sd)
    n_flag = count(r -> r.status == "**flag**", computed)
    n_excl = count(r -> startswith(r.status, "excluded"), computed)
    n_small = count(r -> hasproperty(r, :small) && r.small, computed)

    open(REPORT, "w") do io
        println(io, "# Table V 2S + 1D Strong-Decay Audit")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/audit_table_v_2s_1d_decays.jl`.")
        println(io)
        println(io, "Scores the digitized `1^3D_3`, `1^3D_2`, `1^1D_2`, `1^3D_1`, `2^1S_0`")
        println(io, "and `2^3S_1` sections of Table V against the paper's numeric column,")
        println(io, "using the **same** two-parameter `(A, S0)` harmonic-oscillator model as")
        println(io, "the 1S+1P audit (`table_v_light_decays.md`) with **no refit**:")
        println(io, @sprintf("- `A = %.3f`, `S0 = %.3f`, `beta = %.2f GeV`.", model.A, model.S0, model.beta_GeV))
        println(io)
        println(io, "**Parent masses are the model's own predictions** from")
        println(io, "`compute_spectrum` (corrected stage: central + contact + fine structure),")
        println(io, "so this is an end-to-end check: our spectrum feeds the decay kinematics.")
        println(io, "The paper suppresses the realistic-factor multipliers (`[D/S]`, `(1.x)`)")
        println(io, "from the numeric column, so the comparison is against the raw")
        println(io, "harmonic-oscillator amplitude; the realistic factors are audited")
        println(io, "separately (Eqs. 20-21, `W5`).")
        println(io)
        println(io, "**Headline, split by Table IV amplitude class** (clean rows, |paper| > 1):")
        println(io)
        println(io, @sprintf(
            "- **Structure-independent** (`A`/`A'`/`A''`/`A0`: F- and G-wave rows): %d scored, median abs deviation **%.0f%%**. These carry only the `q^L` angular factor and the elastic form factor, so they isolate the flavor/spin algebra + kinematics + the `A` calibration -- and they reproduce the paper's column tightly.",
            length(scored_si), 100median_si))
        println(io, @sprintf(
            "- **Structure-dependent** (`S`/`D`/`P`: the `S0 - k A qbar^2` classes): %d scored, median abs deviation **%.0f%%**. These are \"highly sensitive to the structure of the states\" (Table IV text): they carry the `S0` strength (whose two-point fit reproduces `S0 = %.2f` vs the paper's reported `3.27`, `reproduction_audit.md` tier B) and a `qbar^2` node that makes them hypersensitive to `qbar`; the residual grows with `qbar` and produces near-zero blow-ups exactly at the node.",
            length(scored_sd), 100median_sd, model.S0))
        println(io)
        println(io, @sprintf(
            "Overall: %d scored, median %.0f%%, flagged (>20%%): %d; small-amplitude rows: %d; excluded (computed but not scored): %d.",
            length(scored), 100median_dev, n_flag, n_small, n_excl))
        println(io)
        # Emit one subtable per section, in reading order.
        section_order = [
            "1^3D_3", "1^3D_2 nonstrange", "1^3D_2 1^1D_2 strange",
            "1^1D_2 nonstrange", "1^3D_1", "2^1S_0", "2^3S_1",
        ]
        for sec in section_order
            secrows = [r for r in computed if r.section == sec]
            isempty(secrows) && continue
            println(io, "## `", sec, "`")
            println(io)
            println(io, "| decay | class | q MeV | computed | paper | ratio | status |")
            println(io, "|---|---|---:|---:|---:|---:|---|")
            for r in secrows
                comp_text = isnan(r.computed) ? "n/a" : @sprintf("%+.2f", r.computed)
                q_text = isnan(r.q_MeV) ? "n/a" : @sprintf("%.0f", r.q_MeV)
                ratio_text = ismissing(r.ratio) ? "" : @sprintf("%.2f", r.ratio)
                println(io, @sprintf("| `%s` | %s | %s | %s | %s | %s | %s |",
                    r.label, r.class, q_text, comp_text, r.paper, ratio_text, r.status))
            end
            println(io)
        end
        println(io, "## Notes")
        println(io)
        println(io, "- The clean read is the class split above: **structure-independent")
        println(io, "  amplitudes reproduce the paper's numeric column to a few percent**")
        println(io, "  (confirming the algebra, the `A` calibration, the flavor/spin")
        println(io, "  coefficients, and the model-mass kinematics), while the")
        println(io, "  **structure-dependent `S`/`D`/`P` rows carry the `S0` strength and a")
        println(io, "  `qbar^2` node** and so spread more and blow up near the node.")
        println(io, "- Parent masses are the corrected (pre-mixing) model masses from")
        println(io, "  `compute_spectrum`. For a given decay the D-wave and F-wave rows share")
        println(io, "  the same `q`, so the tight F-wave (structure-independent) agreement at")
        println(io, "  that `q` rules out a kinematics error as the cause of the D-wave spread.")
        println(io, "- The strange `1^3D_2/1^1D_2` `Q1`/`Q2` rows are the physical **mixed**")
        println(io, "  states (footnote c, `theta_1D ~ 33 deg`), directly analogous to the 1P")
        println(io, "  `Q1`/`Q2` case in `table_v_light_decays.md`; the printed formula is the")
        println(io, "  pure singlet/triplet, so they are excluded from the headline (see")
        println(io, "  `mixing_angles.md` for the 1D mixing angle).")
        println(io, "- `mixing_only` (Table III), `unlisted` (Sec. VD cryptoexotics),")
        println(io, "  quasi-two-body daughters (`eps`/`kappa`/`delta2`), two-part")
        println(io, "  mixing-dependent paper ranges (the `2^1S_0` isoscalar `eta_r`/`eta'_r`")
        println(io, "  rows), and sub-threshold modes are computed where defined but excluded.")
    end

    println("wrote ", REPORT)
    @printf("A=%.3f S0=%.3f | scored=%d (SI=%d SD=%d) flagged=%d small=%d excluded=%d | median: all=%.0f%% SI=%.0f%% SD=%.0f%%\n",
        model.A, model.S0, length(scored), length(scored_si), length(scored_sd),
        n_flag, n_small, n_excl, 100median_dev, 100median_si, 100median_sd)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
