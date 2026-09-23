#!/usr/bin/env julia
# =============================================================================
# Table V reproduction harness
# =============================================================================
# Loads EVERY row of the canonical `table_v_strong_decays.csv`, computes its
# amplitude with the row-oriented transition API (leading-S0 convention), and
# compares against the paper's tabulated MeV^(1/2) column. Writes a per-row
# report with an "N matched / M scoreable" headline and classifies every
# non-match by the specific paper convention it depends on.
#
# Kinematic masses: pinned PDG experimental inputs only. Unassigned states
# remain unavailable; wavefunctions and decay operators retain the GI prescription.

using Pkg
Pkg.activate(@__DIR__)

using CSV
using Printf
using GIModel
using QuarkModelTransitions
using GIPaper: load_table_v, quark_for, experimental_mass, historical_mass

const ROOT = dirname(@__DIR__)
const TABLE = joinpath(ROOT, "data", "raw", "digitized_tables", "table_v_strong_decays.csv")
const REPORT = joinpath(ROOT, "docs", "residual_reports", "table_v_reproduction.md")
const PARAMS_PATH = default_parameters_path()

# --- Flavor content from the parent name -------------------------------------
const NN = Set([
    "rho", "rho2", "rho3", "rhoD", "rhoS", "piprime", "A1", "A2", "A3", "B",
    "delta", "delta2", "g", "D1285", "eps", "f", "H", "omega", "omega2",
    "omega3", "omegaD", "omegaS", "omega_1D2", "eta_r", "h",  # h = f4 (nn): decays to pi pi, K Kbar
])
const SS = Set([
    "phi", "phi2", "phi3", "phiD", "phiS", "phi_1D2", "Hprime", "fprime",
    "epsprime", "hprime", "etaprime_r", "E1420",  # h' = f4' (ss): decays to K Kbar, eta eta
])
const STRANGE = Set([
    "Kstar", "Kstar2", "Kstar3", "Kstar4", "KstarD", "KstarS", "Kprime",
    "kappa", "Q1", "Q2", "Q1_1D2", "Q2_1D2",
])
const CHARMED = Set(["Dstar0", "Dstarplus", "Kstar_c", "kappa_c", "Q1c", "Q2c"])

function parent_flavor(parent)
    parent in NN && return (:q, :q)
    parent in SS && return (:s, :s)
    parent in STRANGE && return (:q, :s)
    parent in CHARMED && return (:c, :d)
    error("unknown parent flavor for `$parent`")
end

# Heavy-quark share of the constituent mass, r = m_Q/(m_Q + m_q), which selects
# the decay form factor (see GIModel.spatial_overlap). Read from the same TOML
# as everything else — the value is NOT frozen in src.
#
# Note what this makes visible: the paper applies its footnote-d unequal-mass
# form factor to the CHARMED rows only. Strange rows are also unequal-mass
# (m_s = 419 vs m_ud = 220 MeV, r = 0.66) yet keep the equal-mass SHO Gaussian.
# So 0.5 here encodes the paper's convention, not the kinematics.
const _QUARK_MASSES = load_quark_masses(PARAMS_PATH)
const CHARM_FRACTION = _QUARK_MASSES["c"] / (_QUARK_MASSES["c"] + _QUARK_MASSES["d"])
heavy_fraction_of(parent) = parent in CHARMED ? CHARM_FRACTION : 0.5

# --- (n, L, multiplicity, J) from the section header (+ parent for mixed Q's) -
function parent_level(section, parent)
    s = section
    if occursin("1P_1", s)   # strange/charmed 1P: Q1 = pure singlet, Q2 = pure triplet
        parent in ("Q1", "Q1c") && return (1, "P", 1, 1)
        parent in ("Q2", "Q2c") && return (1, "P", 3, 1)
    end
    if s == "1^3D_2 1^1D_2 strange"
        parent == "Q1_1D2" && return (1, "D", 1, 2)
        parent == "Q2_1D2" && return (1, "D", 3, 2)
    end
    s in ("1^3S_1", "1^3S_1 charmed") && return (1, "S", 3, 1)
    s in ("1^3P_2", "1^3P_2 charmed") && return (1, "P", 3, 2)
    s == "1^3P_1 nonstrange" && return (1, "P", 3, 1)
    s == "1^1P_1 nonstrange" && return (1, "P", 1, 1)
    s in ("1^3P_0", "1^3P_0 charmed") && return (1, "P", 3, 0)
    s == "1^3D_3" && return (1, "D", 3, 3)
    s == "1^3D_2 nonstrange" && return (1, "D", 3, 2)
    s == "1^1D_2 nonstrange" && return (1, "D", 1, 2)
    s == "1^3D_1" && return (1, "D", 3, 1)
    s == "1^3F_4" && return (1, "F", 3, 4)
    s == "2^1S_0" && return (2, "S", 1, 0)
    s == "2^3S_1" && return (2, "S", 3, 1)
    error("unknown level for section `$section`, parent `$parent`")
end

# All charge conventions and unavailable assignments live in the shared registry.
const DAUGHTER_MASS = Dict(label => something(experimental_mass("V", label), NaN)
    for label in unique(String(r.daughter1) for r in CSV.File(TABLE) if !ismissing(r.daughter1)))
for r in CSV.File(TABLE)
    ismissing(r.daughter2) && continue
    label = String(r.daughter2)
    DAUGHTER_MASS[label] = something(experimental_mass("V", label), NaN)
end
kinematic_mass(label) = something(experimental_mass("V", label), NaN)
modern_momentum(M,m1,m2) = all(isfinite,(M,m1,m2)) ? decay_momentum(M,m1,m2) : NaN

function historical_q(parent,d1,d2)
    masses = [historical_mass("V",label) for label in (parent,d1,d2)]
    any(isnothing,masses) && return NaN
    return 1000decay_momentum(masses...)
end

const QUASI_TWO_BODY = Set(["delta2", "eps", "kappa"])

# Same-J mixing sections: the printed Q1/Q2 formulas are the PURE singlet/triplet
# amplitudes (footnotes b, c); the numeric column is for the physical mixed
# states. Each entry: (paper angle deg, singlet parent, triplet parent,
# flavor sector, singlet level, triplet level).
struct MixSpec
    theta::Float64
    q1::String; q2::String
    flavor::Tuple{Symbol,Symbol}
    lvl1::Tuple{Int,String,Int,Int}; lvl2::Tuple{Int,String,Int,Int}
end
const MIX = Dict(
    # theta_1P ~ +34 deg (Fig. 4); K1(1270)/K1(1400) physical masses.
    "1P_1 strange" => MixSpec(34.0, "Q1", "Q2", (:q, :s), (1, "P", 1, 1), (1, "P", 3, 1)),
    # theta_1D ~ +33 deg; explicit experimental assignments in the registry.
    "1^3D_2 1^1D_2 strange" => MixSpec(33.0, "Q1_1D2", "Q2_1D2", (:q, :s), (1, "D", 1, 2), (1, "D", 3, 2)),
    # charm 1P ~ -41 deg; unresolved experimental correspondence is deferred.
    "1P_1 charmed" => MixSpec(-41.0, "Q1c", "Q2c", (:c, :d), (1, "P", 1, 1), (1, "P", 3, 1)),
)
const MIXING_SECTIONS = Set(keys(MIX))

# Parse the paper's amp_MeV string -> (value, kind).
function parse_amp(text)
    t = strip(String(text))
    t == "below_threshold" && return (nothing, :below)
    (t == "~0" || t == "0") && return (0.0, :zero)
    isempty(t) && return (nothing, :blank)
    body = replace(t, r"^[-+]" => "")     # strip a leading sign before scanning
    occursin(r"[0-9]\s*[-+]\s*[0-9]", body) && return (nothing, :range)  # two-part range
    v = tryparse(Float64, replace(t, "+" => ""))
    return isnothing(v) ? (nothing, :range) : (v, :number)
end

function main()
    params, mq = load_parameters_and_quark_masses(PARAMS_PATH)
    q_rho = modern_momentum(kinematic_mass("rho"), kinematic_mass("pi"), kinematic_mass("pi"))
    q_B = modern_momentum(kinematic_mass("B"), kinematic_mass("omega"), kinematic_mass("pi"))
    model = calibrate_strong_decay_model(q_rho, q_B; convention = LeadingS0())

    parent_mass(section, parent) = (parent == "delta" && section == "1^3F_4" ?
        something(experimental_mass("V:1^3F_4", parent), NaN) : kinematic_mass(parent), "PDG 2026")

    rows = NamedTuple[]
    for r in CSV.File(TABLE)
        loaded = load_table_v((r,); heavy_fraction_for = row -> heavy_fraction_of(String(row.parent)))
        if isempty(loaded)
            # `mixing_only` and Sec.-V-D context rows are part of the canonical
            # 220-row table even though they do not define a pure-state decay
            # operator.  They must remain visible in a coverage audit rather
            # than silently disappearing from the report.
            section = String(r.section); parent = String(r.parent)
            d1 = String(r.daughter1); d2 = String(r.daughter2)
            M, msrc = parent_mass(section, parent)
            m1 = get(DAUGHTER_MASS, d1, NaN); m2 = get(DAUGHTER_MASS, d2, NaN)
            q = modern_momentum(M, m1, m2)
            pv, kind = parse_amp(r.amp_MeV)
            class_str = String(r.amp_class)
            status = class_str == "mixing_only" ?
                "convention: mixing-only row requires physical isoscalar composition" :
                "no paper value"
            push!(rows, (label = String(r.decay), section = section, class = class_str,
                q_MeV = 1000q, q_historical_MeV = historical_q(parent,d1,d2),
                msrc = msrc, M_GeV = M, implied_M = NaN, computed = NaN,
                paper_text = strip(String(r.amp_MeV)), paper = pv, kind = kind,
                ratio = missing, status = status,
                notes = ismissing(r.notes) ? "" : String(r.notes)))
            continue
        end
        ch = only(loaded)
        class_str = String(r.amp_class)
        section = String(r.section); parent = String(r.parent)
        d1 = String(r.daughter1); d2 = String(r.daughter2)

        M, msrc = parent_mass(section, parent)
        m1 = get(DAUGHTER_MASS, d1, NaN); m2 = get(DAUGHTER_MASS, d2, NaN)
        q = modern_momentum(M, m1, m2)
        a = isfinite(q) ? decay_amplitude(model, ch, q; convention = LeadingS0()) : (total=NaN,)

        pv, kind = parse_amp(r.amp_MeV)
        conf = ismissing(r.confidence) ? "" : String(r.confidence)
        notes = ismissing(r.notes) ? "" : String(r.notes)

        # Classify.
        status, ratio = classify(a.total, q, pv, kind, section, d1, d2, conf, notes)
        # For off rows on a structure-independent class, report the parent mass
        # that would reproduce the paper (proves a MISS is a kinematics effect).
        implied = NaN
        if (startswith(status, "MISS") || startswith(status, "near")) &&
           ch.class in STRUCT_INDEP && !isnothing(pv) && !isnan(m1) && !isnan(m2)
            implied = implied_parent_mass(model, ch, m1, m2, pv)
        end
        push!(rows, (label = String(r.decay), section = section, class = class_str,
            q_MeV = 1000q, q_historical_MeV = historical_q(parent,d1,d2), msrc = msrc, M_GeV = M, implied_M = implied,
            computed = a.total, paper_text = strip(String(r.amp_MeV)),
            paper = pv, kind = kind, ratio = ratio, status = status, notes = notes))
    end

    # Score the same-J mixing sections via the singlet/triplet rotation, in
    # place of the deferred placeholders from the main loop.
    filter!(r -> !(r.section in MIXING_SECTIONS), rows)
    append!(rows, mixing_pass(model))

    write_report(model, rows, q_rho, q_B)
    summarize(model, rows)
end

# Physical mixed-state amplitude via the Fig.-4 rotation:
#   A(Q1) =  cos(t) A_singlet + sin(t) A_triplet   (at the Q1 mass)
#   A(Q2) = -sin(t) A_singlet + cos(t) A_triplet   (at the Q2 mass)
# The singlet/triplet configuration amplitudes reuse each row's own coefficient,
# class and qbar power (footnotes b, c), evaluated at the physical parent's q.
function mixing_pass(model)
    out = NamedTuple[]
    for (section, spec) in MIX
        M1 = kinematic_mass(spec.q1)
        M2 = kinematic_mass(spec.q2)
        singlet = Dict{Tuple{String,String,Int},Any}()
        triplet = Dict{Tuple{String,String,Int},Any}()
        for r in CSV.File(TABLE)
            String(r.section) == section || continue
            isempty(load_table_v((r,); heavy_fraction_for = row -> heavy_fraction_of(String(row.parent)))) && continue
            key = (String(r.daughter1), String(r.daughter2), Int(r.qbar_power))
            String(r.parent) == spec.q1 && (singlet[key] = r)
            String(r.parent) == spec.q2 && (triplet[key] = r)
        end
        # configuration amplitude of a row at a given parent mass
        cfg_amp(r, M) = begin
            m1 = get(DAUGHTER_MASS, String(r.daughter1), NaN)
            m2 = get(DAUGHTER_MASS, String(r.daughter2), NaN)
            q = modern_momentum(M, m1, m2)
            ch = only(load_table_v(
                (r,);
                heavy_fraction_for = row -> heavy_fraction_of(String(row.parent)),
            ))
            (isfinite(q) ? decay_amplitude(model, ch, q; convention = LeadingS0()).total : NaN, q,
             String(r.daughter1), String(r.daughter2))
        end
        t = spec.theta
        for (parent, isQ1) in ((spec.q1, true), (spec.q2, false))
            src = isQ1 ? singlet : triplet
            for (key, r) in src
                srow = isQ1 ? r : get(singlet, key, nothing)
                trow = isQ1 ? get(triplet, key, nothing) : r
                M = isQ1 ? M1 : M2
                as = isnothing(srow) ? 0.0 : cfg_amp(srow, M)[1]
                at = isnothing(trow) ? 0.0 : cfg_amp(trow, M)[1]
                amp = isQ1 ? cosd(t) * as + sind(t) * at : -sind(t) * as + cosd(t) * at
                d1, d2 = String(r.daughter1), String(r.daughter2)
                pv, kind = parse_amp(r.amp_MeV)
                q = modern_momentum(M, get(DAUGHTER_MASS, d1, NaN), get(DAUGHTER_MASS, d2, NaN))
                status, ratio = classify_mix(amp, q, pv, kind, d1, d2)
                push!(out, (label = String(r.decay), section = section,
                    class = "mix($(round(Int,t)))", q_MeV = 1000q, q_historical_MeV = historical_q(parent,d1,d2), msrc = "PDG 2026",
                    M_GeV = M, implied_M = NaN, computed = amp,
                    paper_text = strip(String(r.amp_MeV)), paper = pv, kind = kind,
                    ratio = ratio, status = status, notes = ""))
            end
        end
    end
    return out
end

function classify_mix(amp, q, pv, kind, d1, d2)
    isfinite(q) || return ("convention: experimental mass unavailable", missing)
    kind == :below && return ("convention: below threshold (paper integrates lineshape)", missing)
    kind == :range && return ("convention: mixing-dependent range", missing)
    kind == :blank && return ("no paper value", missing)
    (d1 in QUASI_TWO_BODY || d2 in QUASI_TWO_BODY) &&
        return ("convention: quasi-two-body daughter mass", missing)
    q <= 0 && return ("convention: below threshold at chosen mass", missing)
    ratio = pv == 0 ? missing : amp / pv
    if abs(pv) <= 1.0
        d = abs(amp - pv)
        return (d <= 0.4 ? "match (abs $(round(d,digits=2)))" :
                "MISS (abs $(round(d,digits=2)))", ratio)
    end
    dev = abs(ratio - 1)
    return (dev <= 0.15 ? "match" : dev <= 0.25 ? "near ($(round(100dev))%)" :
            "MISS ($(round(100dev))%)", ratio)
end

# Return (status::String, ratio::Union{Missing,Float64}).
function classify(computed, q, pv, kind, section, d1, d2, conf, notes)
    isfinite(q) || return ("convention: experimental mass unavailable", missing)
    kind == :below && return ("convention: below threshold (paper integrates lineshape)", missing)
    kind == :range && return ("convention: mixing-dependent range", missing)
    kind == :blank && return ("no paper value", missing)
    section in MIXING_SECTIONS && return ("convention: same-J mixing angle (Q1/Q2)", missing)
    (d1 in QUASI_TWO_BODY || d2 in QUASI_TWO_BODY) &&
        return ("convention: quasi-two-body daughter mass", missing)
    q <= 0 && return ("convention: below threshold at chosen mass", missing)
    if kind == :zero
        return (abs(computed) <= 0.5 ? "match (~0)" : "MISS (paper ~0)", missing)
    end
    # Documented on-page sign mismatches (formula sign disagrees with printed
    # value): score on magnitude, the sign is the paper's to resolve.
    sign_flag = occursin("sign mismatch", notes)
    cmp = sign_flag ? abs(computed) : computed
    ref = sign_flag ? abs(pv) : pv
    tag = sign_flag ? "sign-flag(page) " : ""
    ratio = pv == 0 ? missing : computed / pv
    small = abs(pv) <= 1.0
    if small
        d = abs(cmp - ref)
        return (d <= 0.3 ? "match $(tag)(abs $(round(d,digits=2)))" :
                "MISS $(tag)(abs $(round(d,digits=2)))", ratio)
    end
    dev = abs(cmp / ref - 1)
    return (dev <= 0.10 ? "match $(tag)"|>strip :
            dev <= 0.20 ? "near $(tag)($(round(100dev))%)" :
            "MISS $(tag)($(round(100dev))%)", ratio)
end

# Structure-independent (A-class) amplitude is a pure function of the parent
# mass M (coefficient, A, qbar^L, form factor). Invert it: scan M for the value
# that reproduces the paper's number, to show a MISS is a parent-mass effect.
const STRUCT_INDEP = Set([:A, :A0, :Aprime, :Adoubleprime, :A_c])

# Physics interpretation belongs beside the authoritative residual calculation,
# not in downstream publications. Keys include the section because historical
# labels such as `K*` are reused for several multiplets.
const MISS_DIAGNOSES = Dict(
    ("1^3P_1 nonstrange", "E -> [(K pi)_K* Kbar]_D") =>
        "D-wave threshold sensitivity; the paper amplitude requires a larger parent mass.",
    ("1^1P_1 nonstrange", "H -> [rho pi]_D") =>
        "The q² factor makes the result sensitive to the historical parent-mass assignment.",
    ("1^1P_1 nonstrange", "H' -> [(K pi)_K* Kbar]_D") =>
        "Near-threshold q² suppression; the modern parent assignment gives too little phase space.",
    ("1^1D_2 nonstrange", "phi -> [K* Kbar]_F") =>
        "The q³ F-wave factor amplifies the parent-mass difference.",
    ("2^3S_1", "rhoS -> (K pi)_K* Kbar") =>
        "A small radial-transition amplitude is dominated by node cancellation and the modern state assignment.",
    ("2^3S_1", "omegaS -> omega eta") =>
        "Radial-node cancellation magnifies the changed breakup momentum.",
    ("2^3S_1", "K*S -> K eta") =>
        "The 2S overlap lies near a node and is sensitive to the modern parent and daughter masses.",
    ("2^3S_1", "K*S -> rho K") =>
        "The 2S overlap lies near a node and is sensitive to the modern parent and daughter masses.",
    ("2^3S_1", "K*S -> omega K") =>
        "The same node-sensitive spatial overlap as rho K appears with a different flavor coefficient.",
    ("2^3S_1", "K*S -> K* pi") =>
        "The 2S overlap and changed breakup momentum jointly shift the amplitude.",
    ("1^3F_4", "K* -> rho K") =>
        "The structure-independent amplitude is reproduced by a higher parent mass than the modern assignment.",
    ("1P_1 strange", "Q1 -> [(K pi)_K* pi]_D") =>
        "Cancellation between rotated singlet/triplet amplitudes is compounded by the changed breakup momentum.",
    ("1P_1 strange", "Q2 -> [(pi pi)_rho K]_S") =>
        "The same-J rotation produces a cancellation-sensitive amplitude; historical and modern momenta also differ.",
)

function implied_parent_mass(model, ch, m1, m2, pv)
    best_M, best_err = NaN, Inf
    for M in 0.9:0.001:2.6
        q = modern_momentum(M, m1, m2)
        q <= 0 && continue
        amp = decay_amplitude(model, ch, q; convention = LeadingS0()).total
        err = abs(amp - pv)
        err < best_err && ((best_M, best_err) = (M, err))
    end
    return best_M
end

function summarize(model, rows)
    scoreable = [r for r in rows if !startswith(r.status, "convention") && r.status != "no paper value"]
    matched = count(r -> startswith(r.status, "match"), scoreable)
    near = count(r -> startswith(r.status, "near"), scoreable)
    miss = count(r -> startswith(r.status, "MISS"), scoreable)
    conv = count(r -> startswith(r.status, "convention"), rows)
    @printf("A=%.3f S0=%.3f (leading)\n", model.A, model.S0)
    @printf("scoreable=%d  matched=%d  near=%d  MISS=%d  | convention-deferred=%d  total-amp-rows=%d\n",
        length(scoreable), matched, near, miss, conv, length(rows))
    for r in rows
        startswith(r.status, "MISS") || continue
        imp = isnan(r.implied_M) ? "" :
            @sprintf("  [M used %.3f -> paper implies %.3f GeV]", r.M_GeV, r.implied_M)
        @printf("  MISS  %-34s computed %+7.2f  paper %-8s  %s%s\n",
            r.label, r.computed, r.paper_text, r.status, imp)
    end
end

function write_report(model, rows, q_rho, q_B)
    scoreable = [r for r in rows if !startswith(r.status, "convention") && r.status != "no paper value"]
    matched = count(r -> startswith(r.status, "match"), scoreable)
    near = count(r -> startswith(r.status, "near"), scoreable)
    miss = count(r -> startswith(r.status, "MISS"), scoreable)
    conv = count(r -> startswith(r.status, "convention"), rows)
    open(REPORT, "w") do io
        println(io, "# Table V Strong-Decay Reproduction")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/reproduce_table_v.jl`. Every")
        println(io, "canonical `table_v_strong_decays.csv` amplitude row is accounted for. Rows with")
        println(io, "a pure-state operator are computed with the")
        println(io, "row-oriented `decay_amplitude` API. Each amplitude factorizes as")
        println(io, "`c * X(qbar) * spatial_overlap` (matrix element x numerical overlap); see")
        println(io, "`docs/observable_ledger.md` for the [PAPER]/[DERIVED] provenance of each factor.")
        println(io)
        println(io, @sprintf("- `A = %.3f` from `rho -> pi pi = +12.4` (q = %.0f MeV)", model.A, 1000q_rho))
        println(io, @sprintf("- `S0 = %.3f` from `B -> [omega pi]_S = -11` (q = %.0f MeV), **leading-S0** convention", model.S0, 1000q_B))
        println(io, "- Kinematics use the shared PDG 2026 mass registry, including stated charge averages.")
        println(io, "- NaN means a required experimental assignment is unavailable; no GI/model mass fallback is used.")
        println(io, "- The paper's two amplitude calibration anchors are retained and evaluated at modern kinematics.")
        println(io, "- Fixed beta and mixing angles retain the GI operator prescription. Broad effective daughters are deferred.")
        println(io, "\n## Comparison")
        println(io, "Differences from GI include updated experimental kinematics as well as the overlap/operator calculation.")
        secs = unique(r.section for r in rows)
        for sec in secs
            secrows = [r for r in rows if r.section == sec]
            println(io, "## `", sec, "`")
            println(io)
            println(io, "| decay | class | q MeV | q reference MeV | src | computed | paper | ratio | status |")
            println(io, "|---|---|---:|---:|:-:|---:|---:|---:|---|")
            for r in secrows
                ratio_text = ismissing(r.ratio) ? "" : @sprintf("%.2f", r.ratio)
                st = r.status
                isnan(r.implied_M) || (st *= @sprintf(" [implies M=%.3f]", r.implied_M))
                println(io, @sprintf("| `%s` | %s | %.0f | %s | %s | %+.2f | %s | %s | %s |",
                    r.label, r.class, r.q_MeV, isnan(r.q_historical_MeV) ? "—" : @sprintf("%.0f",r.q_historical_MeV), r.msrc, r.computed, r.paper_text, ratio_text, st))
            end
            println(io)
        end
        println(io, "## Investigated misses")
        println(io)
        println(io, "This is the canonical interpretation ledger for explicit tolerance misses; downstream publications link to this section rather than maintaining a second diagnosis table.")
        println(io)
        println(io, "| decay | computed / paper | kinematic evidence | diagnosis |")
        println(io, "|---|---:|---|---|")
        for r in rows
            startswith(r.status, "MISS") || continue
            evidence = isnan(r.implied_M) ?
                @sprintf("q = %.0f MeV", r.q_MeV) :
                @sprintf("M used %.3f GeV; paper value implies %.3f GeV", r.M_GeV, r.implied_M)
            diagnosis = get(MISS_DIAGNOSES, (r.section, r.label),
                "Residual is retained without a unique attribution.")
            println(io, @sprintf("| `%s` | %+.2f / %s | %s | %s |",
                r.label, r.computed, r.paper_text, evidence, diagnosis))
        end
        @assert count(r -> startswith(r.status, "MISS"), rows) == length(MISS_DIAGNOSES)
        println(io)
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
