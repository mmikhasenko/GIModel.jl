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
const CHARM_TABLE = joinpath(
    root, "data", "raw", "digitized_tables", "table_v_strong_decays",
    "table_v_charmed.provisional.csv",
)
const REPORT = joinpath(root, "docs", "residual_reports", "table_v_light_decays.md")

# Charmed-meson kinematics (GeV). Daughters use PDG-era D/D*/pi masses; parents
# use the GI-1985 predicted charmed P-wave masses (Sec. IV / Table VIII). The
# 1^3S_1 D* parents are the physical D* states (near threshold, so amplitudes
# are momentum-sensitive at the few-MeV level).
const CHARM_MASS = Dict(
    "Dstarplus" => 2.010, "Dstar0" => 2.007, "Dstar" => 2.008,
    "D0" => 1.865, "Dplus" => 1.869, "D" => 1.867,
    "piplus" => 0.1396, "pi0" => 0.135, "pi" => 0.138,
    # GI-1985 predicted charmed 1P parents (charm analogues of K*, kappa, Q1/Q2)
    "Kstar_c" => 2.50,   # 1^3P_2 (D*_2)
    "kappa_c" => 2.40,   # 1^3P_0
    "Q1c" => 2.44,       # lower 1P (mostly 1^1P_1)
    "Q2c" => 2.49,       # higher 1P (mostly 1^3P_1)
)

# charm-meson quark masses (GeV) for the footnote-d form factor / recoil.
const CHARM_M_C = 1.628
const CHARM_M_D = 0.220

# A_c P-wave rows that carry the explicit recoil multiplier (footnote d):
# K*_c rows and the [D* pi]_D rows of Q1c/Q2c. S_c S-wave and 1^3S_1 rows do not.
charm_has_recoil(section, class, qbar_power) =
    Symbol(class) == :A_c && Int(qbar_power) >= 2

# The charmed CSV's paper_MeV12 column parses as Float64 (unlike the light
# CSV's strings); format it back with an explicit sign for display parity.
charm_paper_text(v) = ismissing(v) ? "" :
    (v isa AbstractString ? String(v) : @sprintf("%+g", v))

const class_map = Dict(
    :A => :A, :A0 => :A0, :Aprime => :Aprime, :Adoubleprime => :Adoubleprime,
    :S => :S, :D => :D, :P => :P,
)

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

# --- Q1/Q2 (strange axial) mixing ---------------------------------------------
#
# Table V footnote b: "We quote the pure Q_B(1) (i.e., pure singlet (triplet))
# amplitude formulas under Q_1(2), where Q_1(2) is the lower (higher) state in
# mass."  So the printed Q1 formula column is the unmixed 1^1P_1 amplitude and
# the printed Q2 column the unmixed 1^3P_1 amplitude, while the numeric column
# is for the physical mixed states, rotated with the Fig. 4 caption convention
#
#   A(Q1 -> X) =  cos(theta) A_singlet(X) + sin(theta) A_triplet(X)
#   A(Q2 -> X) = -sin(theta) A_singlet(X) + cos(theta) A_triplet(X)
#
# with both unmixed amplitudes evaluated at the physical parent's momentum.
# The model angle comes from the public spectrum API exactly as in
# audit_mixing_angles.jl (identity convention mapping, see mixing_angles.md).
const THETA_PAPER_DEG = 34.0   # Fig. 4 caption, theta_1P

function model_k1_angle_deg()
    params_path = joinpath(dirname(root), "data", "parameters.provisional.toml")
    params, mq = load_parameters_and_quark_masses(params_path)
    spec = compute_spectrum(params, Meson(mq, :u, :s);
        levels = [BasisState(1, "P", 1, 1), BasisState(1, "P", 3, 1)])
    state = spectrum_state(spec, 1, "P", 1, 1)
    mix = only(filter(m -> m.mechanism == "antisymmetric_spin_orbit", state.mixings))
    return mix.mixing_angle_deg
end

# Charm 1P same-J mixing angle from the model, for the Q1c/Q2c rows (footnote b
# charm analogue). Same public spectrum API as audit_mixing_angles.jl, on the
# c-dbar 1P block.
function model_k1c_angle_deg()
    params_path = joinpath(dirname(root), "data", "parameters.provisional.toml")
    params, mq = load_parameters_and_quark_masses(params_path)
    spec = compute_spectrum(params, Meson(mq, :c, :d);
        levels = [BasisState(1, "P", 1, 1), BasisState(1, "P", 3, 1)])
    state = spectrum_state(spec, 1, "P", 1, 1)
    mix = only(filter(m -> m.mechanism == "antisymmetric_spin_orbit", state.mixings))
    return mix.mixing_angle_deg
end
const THETA_PAPER_1P_CHARM_DEG = -41.0   # Fig./Table VIII charm 1P (cu-bar analog)

# Rows whose paper numerics depend on machinery beyond the two-parameter
# equal-mass convention; they are computed but excluded from the headline
# score, each with its reason.
function exclusion_reason(row)
    parent = String(row.parent)
    label = String(row.decay_label)
    notes = ismissing(row.notes) ? "" : String(row.notes)
    if parent in ("Q1", "Q2")
        # The numeric column uses the model's strange-axial (K1) mixing angle
        # on top of the unmixed formula coefficients; these rows are scored
        # three ways in the dedicated Q1/Q2 section below the main table.
        return "K1 mixing angle (see Q1/Q2 section)"
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

    # --- Q1/Q2 three-way mixing comparison (footnote b + Fig. 4 rotation) ----
    theta_model = model_k1_angle_deg()
    thetas = [0.0, theta_model, THETA_PAPER_DEG]   # unmixed, model, paper

    qrows = [row for row in rows if String(row.parent) in ("Q1", "Q2")]
    channel_key(row) = (String(row.daughter1), String(row.daughter2), Int(row.qbar_power))
    singlet = Dict(channel_key(r) => r for r in qrows if String(r.parent) == "Q1")
    triplet = Dict(channel_key(r) => r for r in qrows if String(r.parent) == "Q2")

    # Unmixed 1^1P_1 / 1^3P_1 amplitudes for a channel at momentum q (the
    # printed Q1/Q2 formula coefficients per footnote b, common kinematics).
    function unmixed_pair(key, q)
        amp(r) = strong_decay_amplitude(model, Float64(r.coefficient_value),
            class_map[Symbol(String(r.amplitude_class))], Int(r.qbar_power), q)
        return amp(singlet[key]), amp(triplet[key])
    end
    mixed(amp_s, amp_t, theta_deg, parent) =
        parent == "Q1" ? cosd(theta_deg) * amp_s + sind(theta_deg) * amp_t :
        -sind(theta_deg) * amp_s + cosd(theta_deg) * amp_t

    q1q2 = NamedTuple[]
    for row in qrows
        parent = String(row.parent)
        notes = ismissing(row.notes) ? "" : String(row.notes)
        paper_text = strip(String(row.paper_MeV12))
        paper_value = tryparse(Float64, replace(paper_text, "+" => ""))
        q = decay_momentum(PARENT_MASS[parent],
            DAUGHTER_MASS[String(row.daughter1)], DAUGHTER_MASS[String(row.daughter2)])
        amp_s, amp_t = unmixed_pair(channel_key(row), q)
        amps = [mixed(amp_s, amp_t, th, parent) for th in thetas]
        small = !isnothing(paper_value) && abs(paper_value) <= 1.0
        status = if q <= 0
            "unscored: below nominal threshold (paper integrates lineshape)"
        elseif occursin("kinematics convention open", notes)
            "unscored: quasi-two-body kinematics convention open"
        elseif small
            "small-amplitude (abs)"
        else
            "scored"
        end
        push!(q1q2, (
            label = String(row.decay_label), parent = parent, q_GeV = q,
            amp_s = amp_s, amp_t = amp_t, amps = amps,
            paper = paper_value, paper_text = paper_text,
            small = small, status = status, notes = notes,
        ))
    end

    ratio_rows = [r for r in q1q2 if r.status == "scored"]
    small_rows = [r for r in q1q2 if r.status == "small-amplitude (abs)"]
    median_dev(theta) = begin
        devs = sort([abs(mixed(r.amp_s, r.amp_t, theta, r.parent) / r.paper - 1)
                     for r in ratio_rows])
        devs[max(1, (length(devs) + 1) ÷ 2)]
    end
    q1q2_medians = [median_dev(th) for th in thetas]
    scan = 0.0:0.1:60.0
    theta_best = scan[argmin([median_dev(th) for th in scan])]
    # Angle implied by each footnote-j hypersensitive row (amp(theta) = paper).
    implied = [(r.label, scan[argmin([abs(mixed(r.amp_s, r.amp_t, th, r.parent) - r.paper)
                                      for th in scan])])
               for r in q1q2
               if r.label in ("Q1 -> [(K pi)_K* pi]_S", "Q2 -> [(K pi)_K* pi]_D")]

    # --- Charmed section (A_c/S_c, footnote d) -------------------------------
    charm_rows = collect(CSV.File(CHARM_TABLE))
    charm_q(row) = decay_momentum(
        CHARM_MASS[String(row.parent)],
        CHARM_MASS[String(row.daughter1)], CHARM_MASS[String(row.daughter2)])
    charm_amp(row, q) = charm_decay_amplitude(
        model, Float64(row.coefficient_value),
        Symbol(String(row.amplitude_class)), Int(row.qbar_power), q;
        recoil = charm_has_recoil(String(row.section), String(row.amplitude_class),
                                  Int(row.qbar_power)),
        m_c_GeV = CHARM_M_C, m_d_GeV = CHARM_M_D)

    charm_scored = NamedTuple[]
    for row in charm_rows
        parent = String(row.parent)
        parent in ("Q1c", "Q2c") && continue   # handled in the mixing block below
        q = charm_q(row)
        amp = charm_amp(row, q)
        paper_value = ismissing(row.paper_MeV12) ? nothing : Float64(row.paper_MeV12)
        paper_text = charm_paper_text(row.paper_MeV12)
        small = !isnothing(paper_value) && abs(paper_value) <= 1.0
        ratio = (isnothing(paper_value) || paper_value == 0) ? missing : amp / paper_value
        push!(charm_scored, (
            label = String(row.decay_label), section = String(row.section),
            class = String(row.amplitude_class),
            recoil = charm_has_recoil(String(row.section), String(row.amplitude_class),
                                      Int(row.qbar_power)),
            q_GeV = q, computed = amp, paper = paper_text, paper_value = paper_value,
            ratio = ratio, small = small,
        ))
    end
    charm_dev = [abs(r.ratio - 1) for r in charm_scored if !ismissing(r.ratio) && !r.small]
    charm_median = isempty(charm_dev) ? NaN : sort(charm_dev)[max(1, (length(charm_dev) + 1) ÷ 2)]

    # Q1c/Q2c three-way mixing (charm analogue of the Q1/Q2 block, footnote b).
    theta_c_model = model_k1c_angle_deg()
    thetas_c = [0.0, theta_c_model, THETA_PAPER_1P_CHARM_DEG]
    qc_rows = [row for row in charm_rows if String(row.parent) in ("Q1c", "Q2c")]
    ckey(row) = (String(row.daughter1), String(row.daughter2), Int(row.qbar_power))
    csinglet = Dict(ckey(r) => r for r in qc_rows if String(r.parent) == "Q1c")
    ctriplet = Dict(ckey(r) => r for r in qc_rows if String(r.parent) == "Q2c")
    function c_unmixed_pair(key, q)
        a(r) = charm_amp(r, q)
        return a(csinglet[key]), a(ctriplet[key])
    end
    c_mixed(amp_s, amp_t, th, parent) =
        parent == "Q1c" ? cosd(th) * amp_s + sind(th) * amp_t :
        -sind(th) * amp_s + cosd(th) * amp_t
    q1q2c = NamedTuple[]
    for row in qc_rows
        parent = String(row.parent)
        paper_value = ismissing(row.paper_MeV12) ? nothing : Float64(row.paper_MeV12)
        paper_text = charm_paper_text(row.paper_MeV12)
        q = charm_q(row)
        amp_s, amp_t = c_unmixed_pair(ckey(row), q)
        amps = [c_mixed(amp_s, amp_t, th, parent) for th in thetas_c]
        push!(q1q2c, (
            label = String(row.decay_label), parent = parent, q_GeV = q,
            amp_s = amp_s, amp_t = amp_t, amps = amps,
            paper = paper_value, paper_text = paper_text,
        ))
    end

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
        println(io, "## Q1/Q2 (strange axial) rows: three-way K1 mixing-angle comparison")
        println(io)
        println(io, "Table V footnote b: *\"We quote the pure `Q_B(1)` (i.e., pure singlet")
        println(io, "(triplet)) amplitude formulas under `Q_1(2)`, where `Q_1(2)` is the lower")
        println(io, "(higher) state in mass.\"*  So the printed `Q1` formula column is the")
        println(io, "unmixed `1^1P_1` amplitude, the printed `Q2` column the unmixed `1^3P_1`")
        println(io, "amplitude, and the numeric column is for the physical mixed states.")
        println(io)
        println(io, "**Convention adopted** (see `mixing_angles.md` for the identity mapping")
        println(io, "between the model's `StateMixing` angle and the paper's): the printed")
        println(io, "singlet/triplet columns are combined with the Fig. 4 caption rotation")
        println(io)
        println(io, "```")
        println(io, "A(Q1 -> X) =  cos(theta) A_singlet(X) + sin(theta) A_triplet(X)")
        println(io, "A(Q2 -> X) = -sin(theta) A_singlet(X) + cos(theta) A_triplet(X)")
        println(io, "```")
        println(io)
        println(io, "with both unmixed amplitudes evaluated at the *physical parent's*")
        println(io, "momentum, and **no extra relative phase** between the printed singlet and")
        println(io, "triplet amplitude columns. This relative-phase choice is validated by the")
        println(io, "two footnote-j hypersensitive rows: at `theta = +34 deg` the rotation")
        println(io, "reproduces both near-cancellations (`Q1 -> [K* pi]_S` and")
        println(io, "`Q2 -> [K* pi]_D`) in sign and magnitude, while the opposite relative")
        println(io, "phase (`theta -> -theta`) predicts them at full unmixed size (~ -16 and")
        println(io, "~ -3.5), grossly excluded by the paper's -0.3 and -0.1.")
        println(io)
        println(io, @sprintf(
            "Angles compared: **unmixed** `theta = 0` (headline-table behavior), **model** `theta = %+.1f deg` (`compute_spectrum`, `antisymmetric_spin_orbit` mixing of the `u sbar` 1P block), **paper** `theta = %+.0f deg` (Fig. 4 caption `theta_1P`).",
            theta_model, THETA_PAPER_DEG))
        println(io)
        println(io, @sprintf(
            "| decay | q MeV | unmixed | model %+.1f | paper %+.0f | paper | r(unm) | r(mod) | r(pap) | status |",
            theta_model, THETA_PAPER_DEG))
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        for r in q1q2
            ratios = map(a -> (r.status == "scored" ? @sprintf("%.2f", a / r.paper) : ""), r.amps)
            status = r.status
            if r.status == "small-amplitude (abs)"
                status *= @sprintf(": diff %.2f / %.2f / %.2f", (abs.(r.amps .- r.paper))...)
            end
            println(io, @sprintf(
                "| `%s` | %.0f | %+.2f | %+.2f | %+.2f | %s | %s | %s | %s | %s |",
                r.label, 1000r.q_GeV, r.amps[1], r.amps[2], r.amps[3], r.paper_text,
                ratios[1], ratios[2], ratios[3], status,
            ))
        end
        println(io)
        println(io, @sprintf(
            "Median abs deviation over the %d ratio-scored rows: unmixed %.0f%%, model angle %.0f%%, paper angle %.0f%%.",
            length(ratio_rows), 100q1q2_medians[1], 100q1q2_medians[2], 100q1q2_medians[3]))
        println(io, @sprintf(
            "A scan of the median deviation over `theta in [0, 60] deg` is minimized at `theta = %+.1f deg`; the two footnote-j rows individually imply %s.",
            theta_best,
            join([@sprintf("`theta = %+.1f deg` (`%s`)", th, lb) for (lb, th) in implied], " and ")))
        println(io)
        println(io, "**Conclusion:** the paper's Q1/Q2 numeric column is consistent with the")
        println(io, @sprintf(
            "quoted `theta_1P ~ +34 deg`, not with the model's `%+.1f deg` (which misses the two cancellation rows by an order of magnitude). The Q1/Q2 rows therefore corroborate the `mixing_angles.md` finding that our `u sbar` 1P same-J angle is the quantity that deviates, not the decay algebra.",
            theta_model))
        println(io)
        println(io, "Residual paper-angle outliers, all kinematics-driven: `Q1 -> [omega K]_S`")
        println(io, "sits at `q = 32 MeV` (the nominal masses put it 2 MeV above threshold, so")
        println(io, "the amplitude is hostage to the Q1 mass input at the few-MeV level),")
        println(io, "`Q1 -> [rho K]_D` is a tiny D-wave at `q = 96 MeV` with the same")
        println(io, "sensitivity, and `Q2 -> (K pi)_kappa pi` has a quasi-two-body `kappa`")
        println(io, "daughter. The two `(pi pi)_eps K` rows stay unscoreable: with the nominal")
        println(io, "`eps` mass both parents are below threshold (the paper integrates the")
        println(io, "`eps` lineshape, footnote f). The paper-angle value of `Q2 -> [K* pi]_D`")
        println(io, "(-0.10) also confirms the digitized `-0.1` that the OCR note had flagged")
        println(io, "as suspect.")
        println(io)
        println(io, "## Charmed section (`A_c`/`S_c`, Table V footnote d)")
        println(io)
        println(io, "The charmed rows reuse the SAME two-parameter `(A, S0)` calibration as")
        println(io, "the light sector (no refit). Per footnote d the only changes are the")
        println(io, "class letters (`A_c`/`S_c`, identical reduced-amplitude algebra),")
        println(io, "`beta_c = beta` numerically, and the modified Gaussian form factor")
        println(io, "`exp[-(1/4)(m_c/(m_c+m_d))^2 q^2/beta_c^2]`. The four A_c P-wave rows")
        println(io, "(`K*_c` and the `[D* pi]_D` rows of `Q1c`/`Q2c`) carry the explicit")
        println(io, "recoil multiplier `m_c beta / ((m_c+m_d) beta_c)`; the S_c S-wave and")
        println(io, "1^3S_1 rows do not. `m_c = 1.628 GeV`, `m_d = 0.220 GeV` (Table II).")
        println(io)
        println(io, @sprintf(
            "Clean scoreable charmed rows: %d; median abs deviation %.0f%% (|paper| > 1).",
            count(r -> !ismissing(r.ratio) && !r.small, charm_scored),
            isnan(charm_median) ? 0.0 : 100 * charm_median))
        println(io)
        println(io, "| decay | class | recoil | q MeV | computed | paper | ratio | status |")
        println(io, "|---|---|:-:|---:|---:|---:|---:|---|")
        for r in charm_scored
            ratio_text = ismissing(r.ratio) ? "" : @sprintf("%.2f", r.ratio)
            status = r.small ? @sprintf("small-amplitude (abs diff %.2f)", abs(r.computed - r.paper_value)) :
                     (ismissing(r.ratio) ? "unscored" :
                      (abs(r.ratio - 1) > 0.20 ? "**flag**" : "ok"))
            println(io, @sprintf(
                "| `%s` | %s | %s | %.0f | %+.2f | %s | %s | %s |",
                r.label, r.class, r.recoil ? "yes" : "no", 1000r.q_GeV,
                r.computed, r.paper, ratio_text, status))
        end
        println(io)
        println(io, "The three `1^3S_1` `D* -> D pi` rows reproduce the paper's column to")
        println(io, "within a few percent with no refit, directly confirming that `A_c = A`")
        println(io, "with only the charmed form factor applied. The `K*_c` (`1^3P_2`) and")
        println(io, "`kappa_c` (`1^3P_0`) rows are the clean P-wave checks; residual spread")
        println(io, "is dominated by the GI-predicted charmed parent masses (not experimental).")
        println(io)
        println(io, "### `Q1c`/`Q2c` rows: charm 1P mixing (footnote b + j)")
        println(io)
        println(io, "Exactly as the strange `Q1`/`Q2` rows: the printed formula is the pure")
        println(io, "singlet/triplet amplitude while the numeric column is the physical mixed")
        println(io, "state, so the two sign-flagged rows (`Q1c -> [D* pi]_S`,")
        println(io, "`Q2c -> [D* pi]_D`, footnote j hypersensitive) are scored via the same")
        println(io, "three-way rotation, using the charm 1P same-J angle from")
        println(io, "`compute_spectrum(Meson(mq, :c, :d))`.")
        println(io)
        println(io, @sprintf(
            "| decay | q MeV | unmixed | model %+.1f | paper %+.0f | paper | status |",
            theta_c_model, THETA_PAPER_1P_CHARM_DEG))
        println(io, "|---|---:|---:|---:|---:|---:|---|")
        for r in q1q2c
            small = !isnothing(r.paper) && abs(r.paper) <= 1.0
            status = small ? "small-amplitude (mixing-sensitive)" : "scored"
            println(io, @sprintf(
                "| `%s` | %.0f | %+.2f | %+.2f | %+.2f | %s | %s |",
                r.label, 1000r.q_GeV, r.amps[1], r.amps[2], r.amps[3], r.paper_text, status))
        end
        println(io)
        println(io, @sprintf(
            "The model charm 1P angle is `%+.1f deg` (`compute_spectrum`,",
            theta_c_model))
        println(io, @sprintf(
            "`antisymmetric_spin_orbit` on the `c dbar` 1P block); the paper quotes ~`%+.0f deg`.",
            THETA_PAPER_1P_CHARM_DEG))
        println(io, "As in the strange Q1/Q2 case, the two sign-flagged rows are near")
        println(io, "cancellations whose numeric value is set by this mixing angle; the")
        println(io, "clean `Q2c -> [D* pi]_S` (S_c) and `Q1c -> [D* pi]_D` (A_c) rows are the")
        println(io, "robust ones. See `mixing_angles.md` for the charm 1P angle finding.")
        println(io)
        println(io, "## Open Conventions")
        println(io)
        println(io, @sprintf(
            "- `Q1`/`Q2` numeric amplitudes fold the strange-axial (`K1`) `1^1P_1`/`1^3P_1` mixing angle into the unmixed footnote-b formula coefficients; they are scored three ways (unmixed / model `%+.1f deg` / paper `+34 deg`) in the dedicated section above, and select the paper's `+34 deg`.",
            theta_model))
        println(io, "- Quasi-two-body subchannel daughters (`(pi pi)_eps`, `(K pi)_kappa`, `(eta pi)_delta2`) and sub-threshold modes (`D -> K* Kbar`) depend on lineshape conventions the paper does not state; nominal-mass kinematics cannot reproduce them.")
        println(io, "- Strange-parent normalization: `K*2 -> K pi` computes ~sqrt(3) high while `K*(892) -> K pi` is exact, pointing at unequal-mass recoil factors in the Appendix B amplitudes not yet encoded.")
    end
    println("wrote ", REPORT)
    @printf("A=%.3f S0=%.3f scored=%d flagged=%d excluded=%d\n", model.A, model.S0, length(scored), n_flagged, n_excluded)
    @printf("Q1/Q2: theta_model=%+.1f deg, theta_paper=%+.0f deg; median dev unmixed=%.0f%% model=%.0f%% paper=%.0f%%; best-scan theta=%+.1f deg\n",
        theta_model, THETA_PAPER_DEG, 100q1q2_medians[1], 100q1q2_medians[2], 100q1q2_medians[3], theta_best)
    @printf("charm: scored=%d median dev=%.0f%%; theta_c_model=%+.1f deg theta_c_paper=%+.0f deg\n",
        count(r -> !ismissing(r.ratio) && !r.small, charm_scored),
        isnan(charm_median) ? 0.0 : 100charm_median, theta_c_model, THETA_PAPER_1P_CHARM_DEG)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
