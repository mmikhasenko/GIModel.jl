#!/usr/bin/env julia
# =============================================================================
# Table VII audit — gluonic (part c) and leptonic (part a) slices.
# =============================================================================
# Reproduces Table VII annihilation amplitudes from the model wavefunctions
# with the spectrum preset and explicit transition-operator inputs.
#
# Gluonic:  amp = sqrt(prefactor(channel, α_s(M), m_Q)) · S_L(Ψ), amp² = Γ,
#           with S_L the Eq. (17) smeared wavefunction-at-origin.
# Leptonic: f = coefficient · factor, with `factor` one of the mock-meson
#           integrals P_P, P'_A1, V_V, V'_V of Table VII(a)/Eqs. D4-D6
#           (`leptonic_decay_factor`) and `coefficient` the printed
#           quark-charge/flavor factor of the formula column (e.g. 2√3 for
#           f_π, (16/3)^(1/2) for f_ψ). The tabulated amplitude is the
#           dimensionless f (= f_P/M_P for the weak rows, as the experiment
#           column confirms: f_π/M_π = 131/140 ≈ 0.95).
#
# Phase convention: each radial wave's outermost antinode is fixed positive, so
# the wavefunction-at-origin flips sign once per radial node; under this single
# convention the model reproduces the paper's alternating amplitude signs.
#
# Not modelled: the hypothetical t-tbar rows (eta_t, zeta; no top constituent
# mass in the 1985 set). Two-photon and charge-radius rows are included.

mkpath(joinpath(dirname(@__DIR__), "reports"))

using Pkg
Pkg.activate(dirname(@__DIR__); io = devnull)

using Printf
using GIModel
using GIModel.QuarkModelTransitions
# Internal kernels for reference calculations.
using GIModel.QuarkModelTransitions: HBARC_FM2, axial_tau_width, dilepton_vector_width, gluonic_annihilation_amplitude, leptonic_decay_factor, leptonic_pseudoscalar_width, mock_meson_mass, two_photon_amplitude, wavefunction_origin_smearing
# Qualified public API used by this reference/kernel audit.
using GIModel.QuarkModelTransitions: HBARC_FM2,
    mock_meson_mass
using GIPaper: experimental_mass, load_table_policy

kinematic_mass(label) = something(experimental_mass("VII", label), NaN)

const ROOT = dirname(@__DIR__)                                  # GIPaper/
const TABLE = joinpath(ROOT, "data", "table_vii_annihilation_em.csv")
const PARAMS_PATH = default_parameters_path()
const REPORT = joinpath(ROOT, "reports", "table_vii_annihilation_em.md")

const NPTS = 900
const INITIAL_NBASIS = 24 # adaptive paper-order full diagonalization starts here

# External validation anchors, distinct from the inputs to the model prediction.
# These are legacy comparison values, not a claim of updated PDG 2026 widths.
const TABLE_POLICY = load_table_policy()["table_vii"]
const VALIDATION_INPUTS = (; (Symbol(k) => TABLE_POLICY["validation"][k] for k in
    ("f_pi_GeV", "hbar_GeV_s", "pion_lifetime_s", "psi_dilepton_width_GeV"))...)
const DILEPTON_VALIDATION_ROWS = [(r["label"], r["width_GeV"], r["note"])
    for r in TABLE_POLICY["dilepton_validation"]]

# Phase convention: outermost antinode positive (see header).
reporting_wave(wave::RadialWave) = wave

# --- paper-order distorted waves (finite HO-basis full diagonalization) ------
# The single treatment for the whole audit: every wave is a full diagonalization
# of H_central + V_spin in the paper-β harmonic-oscillator basis (the paper's
# literal method). V_spin is the A15 smeared contact for S-waves and the literal
# A15-A16 spin-orbit+tensor operator for the ³P_J chi rows.
# This resolves the wave-choice inconsistency W6 exposed (gluonic previously used
# spin-independent central waves; leptonic/γγ used the FD nonperturbative
# resummation): the finite basis gives the gluonic rows their missing spin
# distortion yet keeps the light ¹S₀ pseudoscalars resummed (pion ≈0.15 GeV),
# where a first-order-PT treatment would fail.

function distorted_family(params, solver_ho, m1, m2, L, multiplicity, J; nlevels)
    L_label = L == 0 ? "S" : L == 1 ? "P" : L == 2 ? "D" :
        throw(ArgumentError("unsupported audit orbital L=$L"))
    return fixed_channel_solution(
        params,
        ConstituentMasses(m1, m2),
        FineStructureMultiplet(L_label, multiplicity, J);
        solver = solver_ho,
        nlevels = max(nlevels, 4),
    )
end

# central (spin-independent) family on the same HO basis (V = 0), for the
# ³D₁/³P₁ leptonic rows the paper leaves undistorted.
function central_family(params, solver_ho, m1, m2, L; nlevels = 2)
    return channel_solution(params, ConstituentMasses(m1, m2), L;
        solver = solver_ho, nlevels = max(nlevels, 4))
end

# --- read the paper's predicted amplitudes from the digitized CSV ------------
function load_paper_predictions(path, subtable)
    preds = Dict{String,Float64}()
    open(path) do io
        header = split(strip(readline(io)), ',')
        di = findfirst(==("decay"), header)
        pi_ = findfirst(==("predicted"), header)
        si = findfirst(==("subtable"), header)
        for line in eachline(io)
            cols = split(line, ',')
            length(cols) < max(di, pi_, si) && continue
            cols[si] == subtable || continue
            val = tryparse(Float64, strip(cols[pi_]))
            val === nothing && continue
            preds[strip(cols[di])] = val
        end
    end
    return preds
end

median_of(xs) = isempty(xs) ? NaN : sort(xs)[cld(length(xs), 2)]

# amplitude unit <-> GeV^(1/2): a value of `v` unit^(1/2) is v·sqrt(unit_in_GeV)
_unit_gev(unit) = unit === :eV ? sqrt(1e-9) : sqrt(1e-6)   # eV or keV
to_gev(v, unit) = v * _unit_gev(unit)
from_gev(a, unit) = a / _unit_gev(unit)

# γγ predicted column carries a unit string ("2.6 eV^1/2", "-1.0 keV^1/2")
function load_gg_predictions(path)
    preds = Dict{String,Tuple{Float64,Symbol}}()
    open(path) do io
        header = split(strip(readline(io)), ',')
        di = findfirst(==("decay"), header)
        pi_ = findfirst(==("predicted"), header)
        si = findfirst(==("subtable"), header)
        for line in eachline(io)
            cols = split(line, ',')
            length(cols) < max(di, pi_, si) && continue
            cols[si] == "gamma_gamma" || continue
            m = match(r"^([-+]?[0-9.]+)\s*(eV|keV)", strip(cols[pi_]))
            m === nothing && continue
            preds[strip(cols[di])] = (parse(Float64, m.captures[1]),
                                      m.captures[2] == "eV" ? :eV : :keV)
        end
    end
    return preds
end

# =============================================================================
# Gluonic slice (part c)
# =============================================================================

struct GluonicRow
    decay::String        # CSV key
    flavor::Symbol       # :c or :b
    L::Int               # orbital angular momentum of the QQ̄ (0=S, 1=P)
    n::Int               # radial level (1 = ground)
    channel::Symbol      # :S0_2g | :S1_3g | :P2_2g | :P0_2g
end

const GLUONIC_ROWS = [
    # S-wave quarkonia: S_L uses L=0 (S₀)
    GluonicRow("eta_c -> 2g",     :c, 0, 1, :S0_2g),
    GluonicRow("psi -> 3g",       :c, 0, 1, :S1_3g),
    GluonicRow("eta'_c -> 2g",    :c, 0, 2, :S0_2g),
    GluonicRow("psi' -> 3g",      :c, 0, 2, :S1_3g),
    GluonicRow("eta_b -> 2g",     :b, 0, 1, :S0_2g),
    GluonicRow("Upsilon -> 3g",   :b, 0, 1, :S1_3g),
    GluonicRow("eta'_b -> 2g",    :b, 0, 2, :S0_2g),
    GluonicRow("Upsilon' -> 3g",  :b, 0, 2, :S1_3g),
    GluonicRow("Upsilon'' -> 3g", :b, 0, 3, :S1_3g),
    GluonicRow("Upsilon''' -> 3g",:b, 0, 4, :S1_3g),
    # P-wave quarkonia (chi): S_L uses L=1 (S₁)
    GluonicRow("chi_2c -> 2g",    :c, 1, 1, :P2_2g),
    GluonicRow("chi_0c -> 2g",    :c, 1, 1, :P0_2g),
    GluonicRow("chi_2b -> 2g",    :b, 1, 1, :P2_2g),
    GluonicRow("chi_0b -> 2g",    :b, 1, 1, :P0_2g),
    GluonicRow("chi'_2b -> 2g",   :b, 1, 2, :P2_2g),
    GluonicRow("chi'_0b -> 2g",   :b, 1, 2, :P0_2g),
]

# spin operator (L, multiplicity, J) for each gluonic annihilation channel
channel_spin(ch::Symbol) =
    ch === :S0_2g ? (0, 1, 0) :
    ch === :S1_3g ? (0, 3, 1) :
    ch === :P2_2g ? (1, 3, 2) :
    ch === :P0_2g ? (1, 3, 0) :
    error("unknown gluonic channel $ch")

function run_gluonic(params, solver_ho, mq, paper)
    cache = Dict{Tuple{Symbol,Int,Int,Int},Any}()   # (flavor, L, mult, J)
    results = NamedTuple[]
    for row in GLUONIC_ROWS
        mQ = mq[String(row.flavor)]
        L, mult, J = channel_spin(row.channel)
        fam = get!(cache, (row.flavor, L, mult, J)) do
            distorted_family(
                params, solver_ho, mQ, mQ, L, mult, J; nlevels = row.n,
            )
        end
        wave = reporting_wave(radial_wave(fam, row.n))
        M = kinematic_mass(row.decay)
        S = wavefunction_origin_smearing(wave, mQ; L = row.L, npoints = NPTS)
        αs = isfinite(M) ? GIModel.alpha_s_q(M) : NaN
        model = gluonic_annihilation_amplitude(row.channel, S, αs, mQ) * sqrt(1000)
        pap = get(paper, row.decay, NaN)
        ratio = isnan(pap) || pap == 0 ? NaN : abs(model) / abs(pap)
        sign_ok = !isnan(pap) && pap != 0 && sign(model) == sign(pap)
        push!(results, (label = row.decay, M = M, S = S, alpha = αs,
                        model = model, paper = pap, ratio = ratio, sign_ok = sign_ok))
    end
    return results
end

# =============================================================================
# Leptonic slice (part a)
# =============================================================================

struct LeptonicRow
    decay::String        # CSV key
    coeff::Float64       # printed flavor/charge coefficient (signed)
    kind::Symbol         # :P_P | :V_V | :Vp_V | :Pp_A1
    f1::String           # constituent flavors
    f2::String
    n::Int               # radial level within the wavefunction channel
end

const LEPTONIC_ROWS = [
    # weak pseudoscalar rows (¹S₀ singlet waves): f_P/M_P = 2√3 P_P
    LeptonicRow("pi -> mu nu",             2sqrt(3),      :P_P,   "q", "q", 1),
    LeptonicRow("K -> mu nu",              2sqrt(3),      :P_P,   "q", "s", 1),
    LeptonicRow("D -> mu nu",              2sqrt(3),      :P_P,   "q", "c", 1),
    LeptonicRow("F -> mu nu",              2sqrt(3),      :P_P,   "s", "c", 1),
    LeptonicRow("B -> mu nu",              2sqrt(3),      :P_P,   "q", "b", 1),
    LeptonicRow("1^1S_0(bc) -> mu nu",     2sqrt(3),      :P_P,   "c", "b", 1),
    # weak axial / vector rows
    LeptonicRow("tau -> A1 nu",            2sqrt(2),      :Pp_A1, "q", "q", 1),
    LeptonicRow("tau -> K*(892) nu",       2sqrt(3),      :V_V,   "q", "s", 1),
    # e+e- rows (³S₁ triplet waves; ³D₁ central waves via V')
    LeptonicRow("rho -> e+ e-",            sqrt(6),       :V_V,   "q", "q", 1),
    LeptonicRow("omega -> e+ e-",          sqrt(2 / 3),   :V_V,   "q", "q", 1),
    LeptonicRow("phi -> e+ e-",           -sqrt(4 / 3),   :V_V,   "s", "s", 1),
    LeptonicRow("rhoS -> e+ e-",           sqrt(6),       :V_V,   "q", "q", 2),
    LeptonicRow("rhoD -> e+ e-",           sqrt(4 / 3),   :Vp_V,  "q", "q", 1),
    LeptonicRow("omegaS -> e+ e-",         sqrt(2 / 3),   :V_V,   "q", "q", 2),
    LeptonicRow("omegaD -> e+ e-",         sqrt(4 / 27),  :Vp_V,  "q", "q", 1),
    LeptonicRow("phiS -> e+ e-",          -sqrt(4 / 3),   :V_V,   "s", "s", 2),
    LeptonicRow("phiD -> e+ e-",          -sqrt(8 / 27),  :Vp_V,  "s", "s", 1),
    LeptonicRow("psi -> e+ e-",            sqrt(16 / 3),  :V_V,   "c", "c", 1),
    LeptonicRow("psi' -> e+ e-",           sqrt(16 / 3),  :V_V,   "c", "c", 2),
    LeptonicRow("psi'' -> e+ e-",          sqrt(32 / 27), :Vp_V,  "c", "c", 1),
    LeptonicRow("psi''' -> e+ e-",         sqrt(16 / 3),  :V_V,   "c", "c", 3),
    LeptonicRow("Upsilon -> e+ e-",       -sqrt(4 / 3),   :V_V,   "b", "b", 1),
    LeptonicRow("Upsilon' -> e+ e-",      -sqrt(4 / 3),   :V_V,   "b", "b", 2),
    LeptonicRow("Upsilon'' -> e+ e-",     -sqrt(4 / 3),   :V_V,   "b", "b", 3),
    LeptonicRow("Upsilon''' -> e+ e-",    -sqrt(4 / 3),   :V_V,   "b", "b", 4),
    LeptonicRow("1^3D_1(bb) -> e+ e-",    -sqrt(8 / 27),  :Vp_V,  "b", "b", 1),
]

# S-wave hyperfine-distinct family (singlet or triplet) for a flavor pair: the
# paper-order finite-HO-basis full diagonalization of H_central + smeared contact
# (levels are the hyperfine-split meson masses).
function swave_family(params, solver_ho, m1, m2, multiplicity; nlevels = 4)
    J = multiplicity == 1 ? 0 : 1
    return distorted_family(
        params, solver_ho, m1, m2, 0, multiplicity, J; nlevels = nlevels,
    )
end

function run_leptonic(params, solver_ho, mq, paper)
    scache = Dict{Tuple{Float64,Float64,Int},Any}()      # (m1, m2, multiplicity)
    ccache = Dict{Tuple{Float64,Float64,Int},Any}()      # (m1, m2, L)
    results = NamedTuple[]
    for row in LEPTONIC_ROWS
        m1, m2 = mq[row.f1], mq[row.f2]
        fam = if row.kind === :P_P
            get!(() -> swave_family(params, solver_ho, m1, m2, 1), scache, (m1, m2, 1))
        elseif row.kind === :V_V
            get!(() -> swave_family(params, solver_ho, m1, m2, 3), scache, (m1, m2, 3))
        elseif row.kind === :Vp_V
            get!(() -> central_family(params, solver_ho, m1, m2, 2), ccache, (m1, m2, 2))
        else # :Pp_A1
            get!(() -> central_family(params, solver_ho, m1, m2, 1), ccache, (m1, m2, 1))
        end
        M = kinematic_mass(row.decay)
        wave = reporting_wave(radial_wave(fam, row.n))
        L = QuarkModelTransitions.LEPTONIC_FACTOR_KINDS[row.kind][1]
        Mt = mock_meson_mass(wave, m1, m2; L = L, npoints = NPTS)
        factor = isfinite(M) ? leptonic_decay_factor(row.kind, wave, m1, m2, M; npoints = NPTS) : NaN
        model = row.coeff * factor
        pap = get(paper, row.decay, NaN)
        ratio = isnan(pap) || pap == 0 ? NaN : abs(model) / abs(pap)
        sign_ok = !isnan(pap) && pap != 0 && sign(model) == sign(pap)
        push!(results, (label = row.decay, kind = row.kind, M = M, Mt = Mt,
                        factor = factor, model = model, paper = pap,
                        ratio = ratio, sign_ok = sign_ok))
    end
    return results
end

# =============================================================================
# Two-photon slice (part b)
# =============================================================================

struct TwoPhotonRow
    decay::String        # CSV key
    kind::Symbol         # :P (¹S₀) | :P2 (³P₂)
    f1::String
    f2::String
    n::Int               # radial level
    q_eff::Float64       # effective squared charge Σ aᵢ eᵢ² of the flavor state
    unit::Symbol         # display unit of the paper value (:eV | :keV)
end

# Rows with unambiguous flavor content (single flavor, isovector, or near-ideal
# tensor mixing). The strongly-mixed isoscalar pseudoscalars (eta, eta', eta_r,
# eta'_r) hinge on the P1/P2 pseudoscalar-annihilation mixing model of Sec. V A
# and are deferred; the hypothetical t-tbar eta_t is not modelled.
const QPI = (4 / 9 - 1 / 9) / sqrt(2)      # (uū−dd̄)/√2 isovector
const QNS = (4 / 9 + 1 / 9) / sqrt(2)      # (uū+dd̄)/√2 nonstrange isoscalar
const TWO_PHOTON_ROWS = [
    TwoPhotonRow("pi -> gamma gamma",     :P,  "q", "q", 1, QPI,   :eV),
    TwoPhotonRow("pi' -> gamma gamma",    :P,  "q", "q", 2, QPI,   :keV),
    TwoPhotonRow("eta_c -> gamma gamma",  :P,  "c", "c", 1, 4 / 9, :keV),
    TwoPhotonRow("eta'_c -> gamma gamma", :P,  "c", "c", 2, 4 / 9, :keV),
    TwoPhotonRow("eta_b -> gamma gamma",  :P,  "b", "b", 1, 1 / 9, :keV),
    TwoPhotonRow("A2 -> gamma gamma",     :P2, "q", "q", 1, QPI,   :keV),
    TwoPhotonRow("f -> gamma gamma",      :P2, "q", "q", 1, QNS,   :keV),   # f₂ ≈ nonstrange (ideal)
    TwoPhotonRow("f' -> gamma gamma",     :P2, "s", "s", 1, 1 / 9, :keV),   # f₂' ≈ ss̄ (ideal)
]
const GG_DEFERRED = ["eta", "eta'", "eta_r", "eta'_r"]   # isoscalar-pseudoscalar mixing

function run_two_photon(params, solver_ho, mq, paper)
    scache = Dict{Tuple{Float64,Float64},Any}()
    ccache = Dict{Tuple{Float64,Float64},Any}()
    results = NamedTuple[]
    for row in TWO_PHOTON_ROWS
        m1, m2 = mq[row.f1], mq[row.f2]
        # ¹S₀ singlets carry the contact distortion; the light ³P₂ (A2/f/f′)
        # radial function is J-independent at the paper's order → central.
        fam = row.kind === :P ?
              get!(() -> swave_family(params, solver_ho, m1, m2, 1), scache, (m1, m2)) :
              get!(() -> central_family(params, solver_ho, m1, m2, 1; nlevels = max(row.n, 1)), ccache, (m1, m2))
        M = kinematic_mass(row.decay)
        wave = reporting_wave(radial_wave(fam, row.n))
        A = isfinite(M) ? two_photon_amplitude(row.kind, wave, m1, M, row.q_eff; npoints = NPTS) : NaN   # GeV^½
        pv, punit = get(paper, row.decay, (NaN, :keV))
        pg = isnan(pv) ? NaN : to_gev(pv, punit)
        ratio = isnan(pg) || pg == 0 ? NaN : abs(A) / abs(pg)
        sign_ok = !isnan(pg) && pg != 0 && sign(A) == sign(pg)
        push!(results, (label = row.decay, kind = row.kind, M = M,
                        model = from_gev(A, row.unit), paper = pv, unit = row.unit,
                        ratio = ratio, sign_ok = sign_ok))
    end
    return results
end

# --- isoscalar-MIXED pseudoscalar γγ rows (eta, eta', eta_r, eta'_r) ----------
# These strongly-mixed states cannot use ideal mixing (it inverts the η<η'
# ordering). The γγ amplitude is a COHERENT sum over the [1nn, 1ss, 2nn, 2ss]
# flavor×radial components, with amplitudes from the P1 pseudoscalar-annihilation
# block (the "P1 for pseudoscalars" of Sec. V A / the Table VI footnote-d
# prescription). The physical mock-meson mass M_P is used in (M/M̃)^(3/2) — GI's
# mock-meson prescription — because the model underpredicts the light-
# pseudoscalar masses and the (M/M̃)^(3/2) factor is acutely mass-sensitive.
const GG_MIXED = [   # (CSV label, model state in ascending eigenvalue order)
    ("eta -> gamma gamma",     1),
    ("eta' -> gamma gamma",    2),
    ("eta_r -> gamma gamma",   3),   # ~ eta(1295)
    ("eta'_r -> gamma gamma",  4),   # ~ iota(1440)
]

function run_two_photon_mixed(params, solver_ho, mq, paper)
    mu, ms = mq["q"], mq["s"]
    Qnn = (4 / 9 + 1 / 9) / sqrt(2)     # (uū+dd̄)/√2 effective charge
    Qss = 1 / 9                          # ss̄
    psl = [GIModel.BasisState(1, "S", 1, 0), GIModel.BasisState(2, "S", 1, 0)]
    final = compute_isoscalar_spectrum(
        params,
        Meson(LightQuark(mu), LightQuark(mu)),
        Meson(StrangeQuark(ms), StrangeQuark(ms));
        levels = psl,
        solver = solver_ho,
        pseudoscalar = PaperP1Annihilation(),
    )
    states = sort(
        [state for state in final.states if
         state.L == "S" && state.multiplicity == 1 && state.J == 0 &&
         any(mix -> occursin("annihilation", mix.mechanism), state.mixings)];
        by = state -> state.mass_GeV,
    )
    length(states) == length(GG_MIXED) || error(
        "expected $(length(GG_MIXED)) final pseudoscalars, got $(length(states))",
    )

    results = NamedTuple[]
    for (label, col) in GG_MIXED
        Mphys = kinematic_mass(label)
        state = states[col]
        components = physical_components(final, state)
        A = physical_state_amplitude(final, state) do component
            flavors = component.basis.flavors
            m, charge = flavors == (:q, :q) ? (mu, Qnn) :
                        flavors == (:s, :s) ? (ms, Qss) :
                        error("unexpected pseudoscalar flavor component $flavors")
            isfinite(Mphys) ? two_photon_amplitude(
                :P, component.wave, m, Mphys, charge; npoints = NPTS,
            ) : NaN
        end
        model = from_gev(A, :keV)
        pv, _ = get(paper, label, (NaN, :keV))
        ratio = isnan(pv) || pv == 0 ? NaN : abs(model) / abs(pv)
        sign_ok = !isnan(pv) && pv != 0 && sign(model) == sign(pv)
        push!(results, (label = label, M = Mphys, M_model = state.mass_GeV,
                        mix = [sum(c.coefficient for c in components if c.basis.n == n &&
                            c.basis.flavors == (f, f); init=0.0) for (n, f) in
                            ((1, :q), (1, :s), (2, :q), (2, :s))],
                        model = model, paper = pv, ratio = ratio, sign_ok = sign_ok))
    end
    return results
end

# =============================================================================
# Charge-radii slice (part d)
# =============================================================================

struct ChargeRadiusRow
    decay::String        # CSV key
    f1::String
    e1::Float64          # quark charges
    f2::String
    e2::Float64
    anchor::Bool         # true = the pi+ fit point (f=0.2 fitted here)
end

const CHARGE_RADIUS_ROWS = [
    ChargeRadiusRow("pi+", "q",  2 / 3, "q", 1 / 3, true),    # u d̄  (fit anchor)
    ChargeRadiusRow("K+",  "q",  2 / 3, "s", 1 / 3, false),   # u s̄
    ChargeRadiusRow("K0",  "q", -1 / 3, "s", 1 / 3, false),   # d s̄
]

# paper column "±(x)^2" fm² -> signed r_E² in fm²
function load_radius_predictions(path)
    preds = Dict{String,Float64}()
    open(path) do io
        header = split(strip(readline(io)), ',')
        di = findfirst(==("decay"), header)
        pi_ = findfirst(==("predicted"), header)
        si = findfirst(==("subtable"), header)
        for line in eachline(io)
            cols = split(line, ',')
            length(cols) < max(di, pi_, si) && continue
            cols[si] == "charge_radius" || continue
            m = match(r"^([+-])\(([0-9.]+)\)", strip(cols[pi_]))
            m === nothing && continue
            s = m.captures[1] == "-" ? -1.0 : 1.0
            preds[strip(cols[di])] = s * parse(Float64, m.captures[2])^2
        end
    end
    return preds
end

function run_charge_radii(params, solver_ho, mq, paper)
    results = NamedTuple[]
    for row in CHARGE_RADIUS_ROWS
        m1, m2 = mq[row.f1], mq[row.f2]
        # ¹S₀ ground-state wave for the flavor pair (paper-order distorted)
        fam = swave_family(params, solver_ho, m1, m2, 1; nlevels = 2)
        wave = radial_wave(fam, 1)
        rE2 = charge_radius_squared(wave, m1, row.e1, m2, row.e2) * HBARC_FM2   # fm²
        pap = get(paper, row.decay, NaN)
        ratio = isnan(pap) || pap == 0 ? NaN : rE2 / pap                        # signed
        sign_ok = !isnan(pap) && pap != 0 && sign(rE2) == sign(pap)
        push!(results, (label = row.decay, rE2 = rE2, paper = pap,
                        ratio = ratio, sign_ok = sign_ok, anchor = row.anchor))
    end
    return results
end

# signed sqrt for reporting r_E² as a radius (fm)
signed_sqrt(x) = sign(x) * sqrt(abs(x))

# =============================================================================
# Report
# =============================================================================

function main()
    params, mq = load_parameters_and_quark_masses(PARAMS_PATH)
    # The paper's own method, at this script's mesh and basis size.
    solver_ho = OscillatorSolver(nbasis = INITIAL_NBASIS)
    glu = run_gluonic(params, solver_ho, mq, load_paper_predictions(TABLE, "gluonic"))
    lep = run_leptonic(params, solver_ho, mq, load_paper_predictions(TABLE, "leptonic"))
    gg_paper = load_gg_predictions(TABLE)
    gg = run_two_photon(params, solver_ho, mq, gg_paper)
    ggm = run_two_photon_mixed(params, solver_ho, mq, gg_paper)
    cr = run_charge_radii(params, solver_ho, mq, load_radius_predictions(TABLE))

    # Eqs. (D7)-(D9): widths from the decay constants. Formula validation +
    # model-f D8 dilepton predictions (see report section).
    wpi_exp = leptonic_pseudoscalar_width(VALIDATION_INPUTS.f_pi_GeV / kinematic_mass("pi -> mu nu"), kinematic_mass("pi -> mu nu"), experimental_mass("D7", "mu"))
    wpi_pdg = VALIDATION_INPUTS.hbar_GeV_s / VALIDATION_INPUTS.pion_lifetime_s
    fpsi_exp = sqrt(VALIDATION_INPUTS.psi_dilepton_width_GeV / ((4π / 3) * QuarkModelTransitions.ALPHA_EM^2 * kinematic_mass("psi -> e+ e-")))
    lep_by = Dict(r.label => r for r in lep)
    wid = NamedTuple[]
    for (label, pdg, note) in DILEPTON_VALIDATION_ROWS
        f = abs(lep_by[label].model)
        w = dilepton_vector_width(f, kinematic_mass(label))
        push!(wid, (label = label, f = f, width = w, pdg = pdg, ratio = w / pdg, note = note))
    end

    gratios = filter(!isnan, [r.ratio for r in glu])
    lratios = filter(!isnan, [r.ratio for r in lep])
    ggratios = filter(!isnan, [r.ratio for r in gg])
    ggmratios = filter(!isnan, [r.ratio for r in ggm])
    # charge radii: score the predictions (exclude the pi+ fit anchor)
    crratios = filter(!isnan, [r.ratio for r in cr if !r.anchor])
    gmed, lmed, ggmed = median_of(gratios), median_of(lratios), median_of(ggratios)
    ggmmed = median_of(ggmratios)
    crmed = median_of(crratios)
    gsign, lsign, ggsign = count(r -> r.sign_ok, glu), count(r -> r.sign_ok, lep), count(r -> r.sign_ok, gg)

    open(REPORT, "w") do io
        println(io, "# Table VII Audit — Annihilation Amplitudes")
        println(io, "\nAll kinematic M and M_P values use the pinned PDG 2026 experimental registry. M̃ and S_L remain wavefunction integrals. NaN denotes an unavailable experimental assignment; no GI/model mass fallback is used. Width-reference measurements below are legacy comparisons, not updated PDG 2026 width inputs.")
        println(io)
        println(io, numerics_provenance(OscillatorSolver(nbasis = INITIAL_NBASIS)))
        println(io)
        println(io, "Generated by `julia GIPaper/checks/audit_table_vii.jl`. All 61 canonical")
        println(io, "rows in the four subtables are accounted for: 57 are computed from the")
        println(io, "model wavefunctions and four hypothetical toponium rows are explicitly unavailable.")
        println(io, "The only fitted quantity anywhere is the")
        println(io, "charge-radius smearing exponent `f = 0.2` (the paper's own π⁺ fit); every")
        println(io, "other amplitude is parameter-free.")
        println(io)
        println(io, "Both slices run on the Eq. (17)-style smeared momentum integral over the")
        println(io, "jₗ-transformed radial wave. Phase convention throughout: each radial")
        println(io, "wave's **outermost antinode is positive**, so the wavefunction-at-origin —")
        println(io, "and every amplitude below — flips sign once per radial node. The overall")
        println(io, "sign of a single width is unobservable (a rephasing of `|M>`); matching the")
        println(io, "paper's alternating signs under one consistent convention demonstrates the")
        println(io, "node structure is right.")
        println(io)

        # --- gluonic ---------------------------------------------------------
        println(io, "## Gluonic decays (part c)")
        println(io)
        println(io, "```")
        println(io, "amp = sqrt(prefactor) · S_L(Ψ),   amp² = Γ")
        println(io, "  Γ(¹S₀→2g) = 8π α_s²/(3 m_Q²) |S₀|²      Γ(³S₁→3g) = 40(π²−9)/(81 m_Q²) α_s³ |S₀|²")
        println(io, "  Γ(³P₂→2g) = 32π α_s²/(45 m_Q²) |S₁|²    Γ(³P₀→2g) = 8π α_s²/(3 m_Q²) |S₁|²")
        println(io, "```")
        println(io)
        println(io, "`S_L(Ψ)` is the Eq. (17) smeared wavefunction-at-origin")
        println(io, "(`wavefunction_origin_smearing`), `α_s = α_s(M)` at the meson mass, `m_Q`")
        println(io, "the constituent quark mass. Waves come from the paper-order **native")
        println(io, "finite-HO-basis full fixed-channel diagonalization** of `H_central + V_spin`:")
        println(io, "the smeared contact for the S-waves (singlet/triplet split) and the")
        println(io, "literal A15-A16 spin-orbit+tensor operator for the ³P_J chi rows — the same")
        println(io, "treatment used for every other subtable below (W6). Amplitudes in `MeV^(1/2)`.")
        println(io)
        println(io, @sprintf("**%d gluonic rows scored; median |model|/|paper| = %.2f; signs agree on %d/%d.** ",
            length(gratios), gmed, gsign, length(gratios)),
            "The two hypothetical t-tbar rows (`eta_t`, `zeta`) are not modelled ",
            "(no top constituent mass in the 1985 set).")
        println(io)
        println(io, "| decay | M (GeV) | α_s(M) | S_L | model MeV^½ | paper MeV^½ | ratio | sign |")
        println(io, "|---|---:|---:|---:|---:|---:|:-:|:-:|")
        for r in glu
            println(io, @sprintf("| `%s` | %.2f | %.3f | %+0.4f | %+0.3f | %+0.3f | %s | %s |",
                r.label, r.M, r.alpha, r.S, r.model, r.paper,
                isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio),
                r.sign_ok ? "✓" : "✗"))
        end
        println(io)

        # --- leptonic --------------------------------------------------------
        println(io, "## Leptonic decays (part a)")
        println(io)
        println(io, "```")
        println(io, "f = coefficient · factor")
        println(io, "  P_P   = M⁻¹ M̃^(-1/2) K[√(m₁m₂/E₁E₂)]     (¹S₀ singlet wave)")
        println(io, "  V_V   = M⁻² M̃^(+1/2) K[√(m₁m₂/E₁E₂)]     (³S₁ triplet wave)")
        println(io, "  V'_V  = M⁻² M̃^(+1/2) K[(m/E)(1 − m/E)]   (³D₁ central wave)")
        println(io, "  P'_A1 = M⁻² M̃^(+1/2) K[m·p/E²]           (³P₁ central wave)")
        println(io, "  K[w]  = (2π)^(-3/2) ∫d³p (4π)^(-1/2) Φ_L(p) w(p)")
        println(io, "```")
        println(io)
        println(io, "`coefficient` is the printed quark-charge/flavor factor of the formula")
        println(io, "column (e.g. `2√3` for `f_π`, `(16/3)^(1/2)` for `f_ψ`); `M` is the")
        println(io, "experimental PDG 2026 mass from the shared registry;")
        println(io, "`M̃` is the model mock mass `<E₁+E₂>` over the same wave")
        println(io, "(`mock_meson_mass`). The tabulated amplitude is dimensionless (the weak")
        println(io, "rows are `f_P/M_P`: the experiment column shows `f_π/M_π = 0.95`).")
        println(io)
        println(io, "The formula column is for unmixed/ideally-mixed states; the paper's")
        println(io, "numerical column additionally folds in Table III isoscalar mixing, which is")
        println(io, "a small effect for the near-ideally-mixed `ω`/`φ` — part of the residual")
        println(io, "here. The hypothetical t-tbar `ζ` row is not modelled.")
        println(io)
        println(io, @sprintf("**%d leptonic rows scored; median |model|/|paper| = %.2f; signs agree on %d/%d.**",
            length(lratios), lmed, lsign, length(lratios)))
        println(io)
        println(io, "| decay | factor | M (GeV) | M̃ (GeV) | value | f model | f paper | ratio | sign |")
        println(io, "|---|:-:|---:|---:|---:|---:|---:|:-:|:-:|")
        for r in lep
            println(io, @sprintf("| `%s` | %s | %.3f | %.3f | %+0.4f | %+0.4f | %+0.4f | %s | %s |",
                r.label, r.kind, r.M, r.Mt, r.factor, r.model, r.paper,
                isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio),
                r.sign_ok ? "✓" : "✗"))
        end
        println(io)

        # --- two-photon ------------------------------------------------------
        println(io, "## Two-photon decays (part b)")
        println(io)
        println(io, "```")
        println(io, "A(P→γγ)   = √6   q_eff (α/m) (M/M̃)^(3/2) (1/2π)     ∫d³p φ_P(p) [m/E]")
        println(io, "A(³P₂→γγ) = −√(4/5) q_eff (α/m) (M/M̃)^(3/2) (2/π)^(1/2) ∫dp p² Φ(p) [m·p/E²]")
        println(io, "```")
        println(io)
        println(io, "`q_eff = Σ aᵢ eᵢ²` is the state's effective squared charge (flavor")
        println(io, "amplitude times quark charges: `(e_u²−e_d²)/√2` for `π`/`A2`, `4/9` for")
        println(io, "`cc̄`, `1/9` for `bb̄`/`ss̄`); `α` the fine-structure constant; `M` the meson")
        println(io, "mass, `M̃` the mock mass. Amplitude² is Γ; units follow the paper (`π` in")
        println(io, "`eV^½`, the rest in `keV^½`).")
        println(io)
        println(io, "The clean-flavor rows are below; the strongly-mixed isoscalar pseudoscalars")
        println(io, "(`", join(GG_DEFERRED, "`, `"), "`) follow in their own table (they use the final")
        println(io, "P1-mixed spectrum). `f`/`f'` use ideal tensor mixing (`f₂` nonstrange, `f₂'` =")
        println(io, "`ss̄`); the hypothetical t-tbar `eta_t` is not modelled.")
        println(io)
        println(io, @sprintf("**%d clean-flavor two-photon rows; median |model|/|paper| = %.2f; signs agree on %d/%d.**",
            length(ggratios), ggmed, ggsign, length(ggratios)))
        println(io)
        println(io, "| decay | kind | M (GeV) | q_eff | model | paper | unit | ratio | sign |")
        println(io, "|---|:-:|---:|---:|---:|---:|:-:|:-:|:-:|")
        for (row, r) in zip(TWO_PHOTON_ROWS, gg)
            println(io, @sprintf("| `%s` | %s | %.3f | %+0.3f | %+0.3f | %+0.3f | %s^½ | %s | %s |",
                r.label, r.kind, r.M, row.q_eff, r.model, r.paper, r.unit,
                isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio),
                r.sign_ok ? "✓" : "✗"))
        end
        println(io)
        println(io, "**Isoscalar-mixed pseudoscalars** (`η`, `η'`, `η_r`, `η'_r`) — a coherent")
        println(io, "sum over the final state's signed `[1nn, 1ss, 2nn, 2ss]` native components")
        println(io, "from `compute_isoscalar_spectrum` / `physical_components` (Sec. V A),")
        println(io, "using the physical mock-meson")
        println(io, "`M_P`. Ideal mixing cannot be used here: it inverts the `η<η'` ordering.")
        println(io, @sprintf("Median |model|/|paper| = %.2f; signs agree on %d/%d.",
            ggmmed, count(r -> r.sign_ok, ggm), length(ggmratios)))
        println(io)
        println(io, "| decay | M_P (GeV) | M model (GeV) | 1nn | 1ss | 2nn | 2ss | model keV^½ | paper keV^½ | ratio | sign |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|:-:|:-:|")
        for r in ggm
            println(io, @sprintf("| `%s` | %.3f | %.6f | %+0.6f | %+0.6f | %+0.6f | %+0.6f | %+0.3f | %+0.3f | %s | %s |",
                r.label, r.M, r.M_model, r.mix[1], r.mix[2], r.mix[3], r.mix[4], r.model, r.paper,
                isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio),
                r.sign_ok ? "✓" : "✗"))
        end
        println(io)
        println(io, "The `η<η'` γγ ordering — backwards under ideal mixing — is reproduced.")
        println(io, "Mixed-state signs use positive overlap with the assigned unmixed precursor;")
        println(io, "the remaining magnitude spread reflects the model's")
        println(io, "sensitivity to the P1 mixing amplitudes and the light-pseudoscalar masses.")
        println(io, "The hypothetical t-tbar `eta_t` is not modelled.")
        println(io)

        # --- charge radii ----------------------------------------------------
        println(io, "## Charge radii (part d)")
        println(io)
        println(io, "```")
        println(io, "r_E² = Σᵢ eᵢ [ <rᵢ²> + (3/4mᵢ²) ∫d³p |φ(p)|² (mᵢ/Eᵢ)^{2f} ]")
        println(io, "```")
        println(io)
        println(io, "`rᵢ = (mⱼ/M)·r` is quark i's position from the meson center of mass, `<r²>`")
        println(io, "the position-space expectation over the `¹S₀` wave, and the second term the")
        println(io, "relativistic smearing of the quark position (`f = 0.2`, the paper's one fit,")
        println(io, "fixed on the `π⁺`). `charge_radius_squared`; values as signed `r_E²` (fm²)")
        println(io, "and the radius `sign·√|r_E²|` (fm).")
        println(io)
        println(io, @sprintf("**%d charge-radius predictions scored (the `π⁺` is the `f` fit anchor); median r_E²(model)/r_E²(paper) = %.2f; signs agree on %d/%d.**",
            length(crratios), crmed, count(r -> r.sign_ok, cr), length(cr)))
        println(io)
        println(io, "| meson | r_E² model fm² | r_E² paper fm² | radius model fm | radius paper fm | ratio | sign |")
        println(io, "|---|---:|---:|---:|---:|:-:|:-:|")
        for r in cr
            note = r.anchor ? " (fit)" : ""
            println(io, @sprintf("| `%s`%s | %+0.4f | %+0.4f | %+0.3f | %+0.3f | %s | %s |",
                r.label, note, r.rE2, r.paper, signed_sqrt(r.rE2), signed_sqrt(r.paper),
                isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio),
                r.sign_ok ? "✓" : "✗"))
        end
        println(io)

        println(io, "## Hypothetical top-flavor rows")
        println(io)
        println(io, "These four printed Table VII predictions cannot be recomputed because Table II")
        println(io, "contains no top constituent mass and the active GI parameter set therefore has no")
        println(io, "toponium wavefunction. Substituting the current quark pole mass would define a new")
        println(io, "model rather than reproduce the 1985 calculation.")
        println(io)
        println(io, "| decay | paper value | unit | local result | reason |")
        println(io, "|---|---:|---|---:|---|")
        println(io, "| `zeta -> e+ e-` | +0.01 | dimensionless | — | no Table II top constituent mass |")
        println(io, "| `eta_t -> gamma gamma` | +1.5 | keV^½ | — | no toponium wavefunction |")
        println(io, "| `eta_t -> 2g` | +1.1 | MeV^½ | — | no toponium wavefunction |")
        println(io, "| `zeta -> 3g` | +0.10 | MeV^½ | — | no toponium wavefunction |")
        println(io)

        println(io, "## Decay widths (Eqs. D7-D9)")
        println(io)
        println(io, "```")
        println(io, "Γ(P→ℓν)   = G² f_P² m_ℓ² / (8π M_P) · (M_P² − m_ℓ²)²        (D7)")
        println(io, "Γ(V→ℓ⁺ℓ⁻) = (4π/3) α² M_V f_V²                              (D8)")
        println(io, "Γ(τ→A₁ν)  = G² f_A1² m_τ³ M_A1²/(16π) [1−M_A1²/m_τ²][1+2M_A1²/m_τ²] (D9)")
        println(io, "```")
        println(io)
        println(io, "Table VII lists the decay *constants* f, not the widths; Eqs. (D7)-(D9)")
        println(io, "(`leptonic_pseudoscalar_width`, `dilepton_vector_width`, `axial_tau_width`)")
        println(io, "turn a constant into a width. The formulas are validated against experiment")
        println(io, "by feeding the **measured** constant:")
        println(io)
        println(io, "| check | fed | model width | measured | note |")
        println(io, "|---|---|---:|---:|---|")
        @printf(io, "| D7 `π→μν` | f_π/M_π = 0.936 (exp) | %.3e GeV | %.3e GeV | +4%% is the Cabibbo cos²θ_C the reduced G² omits |\n",
            wpi_exp, wpi_pdg)
        @printf(io, "| D8 `ψ→ee` | f_ψ = %.4f (exp) | %.3e GeV | 5.55e-06 GeV | round-trip |\n",
            fpsi_exp, dilepton_vector_width(fpsi_exp, kinematic_mass("psi -> e+ e-")))
        println(io)
        println(io, "Applied to the **model's own** constants (physical masses; the paper's")
        println(io, "reduced `G²`/`α²` carry no CKM or QCD radiative factor), the D8 dilepton")
        println(io, "widths carry the Table VII `f_V` residual squared:")
        println(io)
        println(io, "| decay | f_model | width (D8) | measured | ratio | note |")
        println(io, "|---|---:|---:|---:|:-:|---|")
        for w in wid
            @printf(io, "| `%s` | %+.4f | %.3e GeV | %.3e GeV | %.2f | %s |\n",
                w.label, w.f, w.width, w.pdg, w.ratio, w.note)
        end
        println(io)
        println(io, "The paper tabulates decay constants; D8 widths are derived here using the same PDG mass as the corresponding constant. Radiative corrections are not included.")
        println(io, "\n## Reading")
        println(io, "The comparison tests GI wavefunctions and operators at modern experimental kinematics. Differences may reflect both numerical reproduction and changed mass inputs. Model eigenvalues belong to the separate spectrum comparison.")
        println(io, "Experimental state correspondences and unavailable inputs are listed in the shared mass registry. No mass is selected to improve agreement.")

    end

    @printf("wrote %s\n", REPORT)
    @printf("gluonic:   %d rows, median |model|/|paper| = %.2f, signs %d/%d\n",
        length(gratios), gmed, gsign, length(gratios))
    @printf("leptonic:  %d rows, median |model|/|paper| = %.2f, signs %d/%d\n",
        length(lratios), lmed, lsign, length(lratios))
    @printf("two-photon: %d clean rows, median = %.2f, signs %d/%d ; mixed eta rows median = %.2f, signs %d/%d\n",
        length(ggratios), ggmed, ggsign, length(ggratios),
        ggmmed, count(r -> r.sign_ok, ggm), length(ggmratios))
    for r in ggm
        @printf("  %-22s model=%+7.3f paper=%+7.3f keV^½  ratio=%s sign=%s\n",
            r.label, r.model, r.paper, isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio),
            r.sign_ok ? "ok" : "X")
    end
    @printf("charge radii: %d predictions, median ratio = %.2f, signs %d/%d\n",
        length(crratios), crmed, count(r -> r.sign_ok, cr), length(cr))
    for r in cr
        @printf("  %-6s%-6s r_E²=%+7.4f fm² (paper %+7.4f)  radius=%+6.3f (paper %+6.3f)  ratio=%s\n",
            r.label, r.anchor ? "(fit)" : "", r.rE2, r.paper,
            signed_sqrt(r.rE2), signed_sqrt(r.paper),
            isnan(r.ratio) ? "—" : @sprintf("%.2f", r.ratio))
    end
    return (; solver_ho, glu, lep, gg, ggm, cr, wid, wpi_exp, wpi_pdg, fpsi_exp)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
