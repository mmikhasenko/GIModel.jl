#!/usr/bin/env julia
# Prototype audit of Table VI photon-decay amplitudes computed from the
# model's own FD wavefunctions (no new fitted constants).
#
# Appendix D mock-meson matrix elements:
#
#   I_i(x,y)   = sqrt(4 M~x M~y)/(M~x + M~y)
#                * ∫ dp p² Φx(p) Φy(p) (1/m_i) (m_i/E_i)^0.7
#   E_n^i(x,y) = | m_i / sqrt(<E_i>_x <E_i>_y) |^0.5 * ∫ dr u_x(r) u_y(r) r^n
#
# with Φ(p) the normalized radial momentum wavefunction (j_L transform of the
# FD reduced radial wave u(r)), E_i = sqrt(m_i² + p²), and the mock mass
# M~ = <E_1> + <E_2>.  Table VI lists M1 moments in units of e/2, so
# μ/μ_N = coefficient * I * M_N.  E1 amplitudes are `formula * sqrt(α q)`
# in MeV^(1/2) with the q-dependence written explicitly in the formula
# column (q in GeV there, MeV inside the square root).
#
# The exponents 0.7 / 0.5 are the paper's, fitted there to ρ→πγ and A2→πγ;
# ρ→πγ therefore doubles as the normalization check of this whole pipeline.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf
using LinearAlgebra

root = dirname(@__DIR__)
using GIModel
const G = GIModel

params_path = joinpath(dirname(root), "data", "parameters.provisional.toml")
params, mq = load_parameters_and_quark_masses(params_path)

const NGRID = 450
const RMAX = 24.0
const M_N_GEV = 0.93827
const ALPHA_EM = 1 / 137.036
const PMAX = 30.0
const NP = 1501

trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i = 1:length(x)-1)

# --- radial solves ----------------------------------------------------------

"""S-wave solve with the nonperturbative smeared contact term for the given
spin multiplicity (1 = pseudoscalar, 3 = vector). Returns physically
normalized reduced radial waves (∫ u² dr = 1) with the GI phase Φ(0) > 0."""
function swave_waves(params, masses; multiplicity, nlevels = 3)
    H, r = G.relativistic_hamiltonian(params, masses, 0; ngrid = NGRID, rmax = RMAX)
    op = G.contact_hyperfine_operator(params, masses, "S", multiplicity, r)
    _, vecs = G.lowest_eigenpairs(Symmetric(Matrix(H) + Matrix(op)), nlevels)
    return physical_waves(vecs, r)
end

"""Central (spin-independent) solve for orbital L; used for the P-wave chi
states whose radial function is unperturbed at the paper's first order."""
function central_waves(params, masses, L; nlevels = 2)
    H, r = G.relativistic_hamiltonian(params, masses, L; ngrid = NGRID, rmax = RMAX)
    _, vecs = G.lowest_eigenpairs(Symmetric(Matrix(H)), nlevels)
    return physical_waves(vecs, r)
end

"""Model (GI predicted) meson masses in GeV: eigenvalues of the same
S-wave Hamiltonian (kinetic sqrt terms carry the constituent rest masses,
so eigenvalues ARE meson masses) plus the nonperturbative contact term.
Reuses `swave_waves`' operators; no new fitted constants."""
function swave_masses(params, masses; multiplicity, nlevels = 3)
    H, r = G.relativistic_hamiltonian(params, masses, 0; ngrid = NGRID, rmax = RMAX)
    op = G.contact_hyperfine_operator(params, masses, "S", multiplicity, r)
    vals, _ = G.lowest_eigenpairs(Symmetric(Matrix(H) + Matrix(op)), nlevels)
    return vals
end

"""Model (GI predicted) meson masses in GeV for the central orbital-L solve
(P-wave chi radial function is J-independent at the paper's first order)."""
function central_masses(params, masses, L; nlevels = 2)
    H, r = G.relativistic_hamiltonian(params, masses, L; ngrid = NGRID, rmax = RMAX)
    vals, _ = G.lowest_eigenpairs(Symmetric(Matrix(H)), nlevels)
    return vals
end

function physical_waves(vecs, r)
    h = r[2] - r[1]
    waves = RadialWaveOnUniformMesh[]
    for n = 1:size(vecs, 2)
        u = vecs[:, n] ./ sqrt(h)              # Euclidean → ∫u² dr = 1
        sum(r .* u) < 0 && (u = -u)            # Φ(0) > 0 phase convention
        push!(waves, RadialWaveOnUniformMesh(u, r, h))
    end
    return waves
end

# --- momentum-space overlaps ------------------------------------------------

struct MomentumWave
    p::Vector{Float64}
    phi::Vector{Float64}
end

function momentum_wave(wave::RadialWaveOnUniformMesh, L::Integer)
    p = collect(range(0.0, PMAX; length = NP))
    phi = [G._momentum_radial_wave(wave, pk, L) for pk in p]
    nrm = sqrt(trapz(p, p .^ 2 .* phi .^ 2))
    return MomentumWave(p, phi ./ nrm)
end

mean_energy(mw::MomentumWave, m) =
    trapz(mw.p, mw.p .^ 2 .* mw.phi .^ 2 .* sqrt.(m^2 .+ mw.p .^ 2))

mock_mass(mw::MomentumWave, m1, m2) = mean_energy(mw, m1) + mean_energy(mw, m2)

"""Appendix-D M1 overlap I_i(x,y) in GeV⁻¹; `m_i` is the emitting quark."""
function I_overlap(mwx::MomentumWave, mwy::MomentumWave, Mx, My, m_i)
    pref = sqrt(4 * Mx * My) / (Mx + My)
    kern = mwx.p .^ 2 .* mwx.phi .* mwy.phi .*
           (1 / m_i) .* (m_i ./ sqrt.(m_i^2 .+ mwx.p .^ 2)) .^ 0.7
    return pref * trapz(mwx.p, kern)
end

"""Appendix-D E1 moment E_n^i(x,y) in GeV⁻ⁿ (symmetric in x, y)."""
function E_moment(wx::RadialWaveOnUniformMesh, wy::RadialWaveOnUniformMesh,
                  Ex_mean, Ey_mean, m_i; n = 1)
    wx.r == wy.r || error("E_moment: meshes differ")
    radial = sum(wx.u .* wy.u .* wx.r .^ n) * wx.h
    return abs(m_i / sqrt(Ex_mean * Ey_mean))^0.5 * radial
end

photon_momentum(M_parent, M_child) = (M_parent^2 - M_child^2) / (2 * M_parent)

# --- assemble the flavor systems --------------------------------------------

m_ud = mq["u"]
m_s = mq["s"]
m_c = mq["c"]
m_b = mq["b"]

struct SWaveSystem
    label::String
    m1::Float64
    m2::Float64
    singlet::Vector{RadialWaveOnUniformMesh}
    triplet::Vector{RadialWaveOnUniformMesh}
    singlet_p::Vector{MomentumWave}
    triplet_p::Vector{MomentumWave}
end

function SWaveSystem(label, m1, m2; nlevels = 3)
    masses = ConstituentMasses(m1, m2)
    s = swave_waves(params, masses; multiplicity = 1, nlevels = nlevels)
    t = swave_waves(params, masses; multiplicity = 3, nlevels = nlevels)
    SWaveSystem(label, m1, m2, s, t,
        [momentum_wave(w, 0) for w in s], [momentum_wave(w, 0) for w in t])
end

println("solving S-wave systems (central + nonperturbative contact) ...")
sys = Dict(
    "nn" => SWaveSystem("nn", m_ud, m_ud),
    "ns" => SWaveSystem("ns", m_ud, m_s),
    "nc" => SWaveSystem("nc", m_ud, m_c),
    "sc" => SWaveSystem("sc", m_s, m_c),
    "nb" => SWaveSystem("nb", m_ud, m_b),
    "sb" => SWaveSystem("sb", m_s, m_b),
    "cc" => SWaveSystem("cc", m_c, m_c),
    "bb" => SWaveSystem("bb", m_b, m_b),
)

"""M1 moment μ/μ_N for V(nV) → P(nP) γ given `(coeff, m_i)` terms."""
function m1_moment(system::SWaveSystem, nP, nV, terms)
    mwx = system.singlet_p[nP]
    mwy = system.triplet_p[nV]
    Mx = mock_mass(mwx, system.m1, system.m2)
    My = mock_mass(mwy, system.m1, system.m2)
    return sum(c * I_overlap(mwx, mwy, Mx, My, m_i) for (c, m_i) in terms) * M_N_GEV
end

# --- M1 rows (mixing-free entries of Table VI) -------------------------------

m1_rows = Vector{Tuple{String,Float64,String}}()

push!(m1_rows, ("rho -> pi gamma (fit of exponent 0.7)",
    m1_moment(sys["nn"], 1, 1, [(1 / 3, m_ud)]), "+0.69"))
push!(m1_rows, ("omega -> pi gamma",
    m1_moment(sys["nn"], 1, 1, [(1.0, m_ud)]), "+2.07"))
# Open-flavor coefficients rebuilt from quark charges with the transition-
# moment rule mu = e_q I_q - e_qbar I_qbar (antiquark charge enters flipped).
# Paper values below are IMAGE-VERIFIED against the printed "Predicted mu"
# column on PDF page 24 (page_images/page-024.png, printed p. 212). The raw
# vision-OCR column (md lines 1124-1133) was displaced DOWN by one row vs the
# decay labels; every printed value read directly off the page confirms the
# earlier shift-correction exactly (K*+ +0.91, K*0 -1.20, D*+ -0.35,
# D*0 +1.78, F* -0.13, B* +1.37, B*0 -0.78, F_b* -0.55). No values changed.
# (Printed labels: the u-bbar row is B*- -> B- gamma, and F_b* carries note f.)
push!(m1_rows, ("K*+ -> K+ gamma   [+2/3 I_d - 1/3 I_s]",
    m1_moment(sys["ns"], 1, 1, [(2 / 3, m_ud), (-1 / 3, m_s)]), "+0.91 (image-verified page 24)"))
push!(m1_rows, ("K*0 -> K0 gamma   [-1/3 I_d - 1/3 I_s]",
    m1_moment(sys["ns"], 1, 1, [(-1 / 3, m_ud), (-1 / 3, m_s)]), "-1.20 (image-verified page 24)"))
push!(m1_rows, ("D*+ -> D+ gamma   [+2/3 I_c - 1/3 I_d]",
    m1_moment(sys["nc"], 1, 1, [(2 / 3, m_c), (-1 / 3, m_ud)]), "-0.35 (image-verified page 24)"))
push!(m1_rows, ("D*0 -> D0 gamma   [+2/3 I_c + 2/3 I_d]",
    m1_moment(sys["nc"], 1, 1, [(2 / 3, m_c), (2 / 3, m_ud)]), "+1.78 (image-verified page 24)"))
push!(m1_rows, ("F* -> F gamma     [+2/3 I_c - 1/3 I_s]",
    m1_moment(sys["sc"], 1, 1, [(2 / 3, m_c), (-1 / 3, m_s)]), "-0.13 (image-verified page 24)"))
push!(m1_rows, ("B*-(u bbar)       [+2/3 I_d - 1/3 I_b]",
    m1_moment(sys["nb"], 1, 1, [(2 / 3, m_ud), (-1 / 3, m_b)]), "+1.37 (image-verified page 24)"))
push!(m1_rows, ("B*0(d bbar)       [-1/3 I_d - 1/3 I_b]",
    m1_moment(sys["nb"], 1, 1, [(-1 / 3, m_ud), (-1 / 3, m_b)]), "-0.78 (image-verified page 24)"))
push!(m1_rows, ("F_b*(b sbar)      [-1/3 I_b - 1/3 I_s]",
    m1_moment(sys["sb"], 1, 1, [(-1 / 3, m_b), (-1 / 3, m_s)]), "-0.55 (image-verified page 24)"))
push!(m1_rows, ("psi -> eta_c gamma",
    m1_moment(sys["cc"], 1, 1, [(4 / 3, m_c)]), "+0.69"))
push!(m1_rows, ("psi' -> eta_c' gamma",
    m1_moment(sys["cc"], 2, 2, [(4 / 3, m_c)]), "+0.68"))
push!(m1_rows, ("Upsilon -> eta_b gamma",
    m1_moment(sys["bb"], 1, 1, [(-2 / 3, m_b)]), "-0.13"))
push!(m1_rows, ("Upsilon' -> eta_b' gamma",
    m1_moment(sys["bb"], 2, 2, [(-2 / 3, m_b)]), "-0.12"))
push!(m1_rows, ("Upsilon'' -> eta_b'' gamma",
    m1_moment(sys["bb"], 3, 3, [(-2 / 3, m_b)]), "-0.12"))

# hindered M1 with recoil term: psi' -> eta_c gamma
let s = sys["cc"]
    q = photon_momentum(3.686, 2.980)
    Mx = mock_mass(s.singlet_p[1], m_c, m_c)
    My = mock_mass(s.triplet_p[2], m_c, m_c)
    I_direct = I_overlap(s.singlet_p[1], s.triplet_p[2], Mx, My, m_c)
    E2 = E_moment(s.singlet[1], s.triplet[2],
        mean_energy(s.singlet_p[1], m_c), mean_energy(s.triplet_p[2], m_c), m_c; n = 2)
    mu = (4 / 3) * (I_direct - q^2 / (24 * m_c) * E2) * M_N_GEV
    push!(m1_rows, ("psi' -> eta_c gamma (hindered, with recoil)", mu, "-0.056"))
end

# --- E1 rows (quarkonium chi systems) ----------------------------------------

println("solving P-wave channels for E1 rows ...")
cc_P = central_waves(params, ConstituentMasses(m_c, m_c), 1; nlevels = 2)
bb_P = central_waves(params, ConstituentMasses(m_b, m_b), 1; nlevels = 2)
cc_P_p = [momentum_wave(w, 1) for w in cc_P]
bb_P_p = [momentum_wave(w, 1) for w in bb_P]

"""E1 amplitude `coeff(q_GeV) * E_1^i * sqrt(alpha * q_MeV)` in MeV^(1/2).
Pass an explicit `q` to override the measured-mass photon momentum."""
function e1_amplitude(wS, mwS, wP, mwP, m_i, coeff_of_q, M_parent, M_child; q = nothing)
    q === nothing && (q = photon_momentum(M_parent, M_child))
    E1 = E_moment(wS, wP, mean_energy(mwS, m_i), mean_energy(mwP, m_i), m_i; n = 1)
    return coeff_of_q(q) * E1 * sqrt(ALPHA_EM * 1000 * q)
end

e1_rows = Vector{Tuple{String,Float64,String}}()

let s = sys["cc"]
    ψ, ψp = s.triplet[1], s.triplet[2]
    ψ_p, ψp_p = s.triplet_p[1], s.triplet_p[2]
    for (name, Mχ, target) in
        [("chi_c2 -> psi gamma", 3.556, "+0.50"),
         ("chi_c1 -> psi gamma", 3.510, "+0.44"),
         ("chi_c0 -> psi gamma", 3.415, "+0.30")]
        amp = e1_amplitude(ψ, ψ_p, cc_P[1], cc_P_p[1], m_c, q -> 4q / 9, Mχ, 3.097)
        push!(e1_rows, (name, amp, target))
    end
    for (name, Mχ, cJ, target) in
        [("psi' -> chi_c2 gamma", 3.556, sqrt(5 / 3), "+0.14"),
         ("psi' -> chi_c1 gamma", 3.510, 1.0, "+0.15"),
         ("psi' -> chi_c0 gamma", 3.415, sqrt(1 / 3), "+0.14")]
        amp = e1_amplitude(ψp, ψp_p, cc_P[1], cc_P_p[1], m_c, q -> cJ * 4q / 9, 3.686, Mχ)
        push!(e1_rows, (name, amp, target))
    end
end

let s = sys["bb"]
    Υ, Υp = s.triplet[1], s.triplet[2]
    Υ_p, Υp_p = s.triplet_p[1], s.triplet_p[2]
    for (name, Mχ, target) in
        [("chi_b2 -> Upsilon gamma", 9.913, "-0.18"),
         ("chi_b1 -> Upsilon gamma", 9.892, "-0.17"),
         ("chi_b0 -> Upsilon gamma", 9.860, "-0.16")]
        amp = e1_amplitude(Υ, Υ_p, bb_P[1], bb_P_p[1], m_b, q -> -2q / 9, Mχ, 9.460)
        push!(e1_rows, (name, amp, target))
    end
    for (name, Mχ, cJ, target) in
        [("Upsilon' -> chi_b2 gamma", 9.913, sqrt(5 / 3), "-0.040"),
         ("Upsilon' -> chi_b1 gamma", 9.892, 1.0, "-0.038"),
         ("Upsilon' -> chi_b0 gamma", 9.860, sqrt(1 / 3), "-0.025")]
        amp = e1_amplitude(Υp, Υp_p, bb_P[1], bb_P_p[1], m_b, q -> -cJ * 2q / 9, 10.023, Mχ)
        push!(e1_rows, (name, amp, target))
    end
end

# --- 2S -> chi_0 momentum convention test ------------------------------------
# The `psi' -> chi_c0 gamma` and `Upsilon' -> chi_b0 gamma` rows sit 20-30%
# high with q from measured 1984 masses. These rows have the largest q of the
# 2S->chi block, so the amplitude (~ q^(3/2) via coeff*sqrt(q)) is the most
# q-sensitive. Test whether the paper used the MODEL predicted masses for the
# parent (2^3S_1) and daughter (1^3P_0) instead. Model masses are eigenvalues
# of the same solves (no new constants).
println("computing model masses for the 2S->chi_0 momentum test ...")
cc_S_masses = swave_masses(params, ConstituentMasses(m_c, m_c); multiplicity = 3, nlevels = 3)
bb_S_masses = swave_masses(params, ConstituentMasses(m_b, m_b); multiplicity = 3, nlevels = 3)
cc_P_masses = central_masses(params, ConstituentMasses(m_c, m_c), 1; nlevels = 2)
bb_P_masses = central_masses(params, ConstituentMasses(m_b, m_b), 1; nlevels = 2)

# rows: (name, target, system, m_i, cJ, coeff_sign, M_parent_meas, M_child_meas,
#        M_parent_model, M_child_model)
chi0_variant_rows = Vector{NTuple{7,Any}}()
let s = sys["cc"]
    q_meas = photon_momentum(3.686, 3.415)
    q_model = photon_momentum(cc_S_masses[2], cc_P_masses[1])
    amp_meas = e1_amplitude(s.triplet[2], s.triplet_p[2], cc_P[1], cc_P_p[1], m_c,
        q -> sqrt(1 / 3) * 4q / 9, 3.686, 3.415; q = q_meas)
    amp_model = e1_amplitude(s.triplet[2], s.triplet_p[2], cc_P[1], cc_P_p[1], m_c,
        q -> sqrt(1 / 3) * 4q / 9, 3.686, 3.415; q = q_model)
    push!(chi0_variant_rows,
        ("psi' -> chi_c0 gamma", "+0.14", q_meas, amp_meas, q_model, amp_model,
         (cc_S_masses[2], cc_P_masses[1])))
end
let s = sys["bb"]
    q_meas = photon_momentum(10.023, 9.860)
    q_model = photon_momentum(bb_S_masses[2], bb_P_masses[1])
    amp_meas = e1_amplitude(s.triplet[2], s.triplet_p[2], bb_P[1], bb_P_p[1], m_b,
        q -> -sqrt(1 / 3) * 2q / 9, 10.023, 9.860; q = q_meas)
    amp_model = e1_amplitude(s.triplet[2], s.triplet_p[2], bb_P[1], bb_P_p[1], m_b,
        q -> -sqrt(1 / 3) * 2q / 9, 10.023, 9.860; q = q_model)
    push!(chi0_variant_rows,
        ("Upsilon' -> chi_b0 gamma", "-0.025", q_meas, amp_meas, q_model, amp_model,
         (bb_S_masses[2], bb_P_masses[1])))
end

reldev(val, paper) = abs(val - paper) / abs(paper)

# --- report ------------------------------------------------------------------

outpath = joinpath(root, "docs", "residual_reports", "table_vi_photon_decays.md")
open(outpath, "w") do io
    println(io, "# Table VI Photon-Decay Prototype Audit")
    println(io)
    println(io, "Generated by `scripts/audit_table_vi_photon_decays.jl` on ",
        Dates.format(Dates.now(), "yyyy-mm-dd"), ".")
    println(io)
    println(io, "Appendix-D mock-meson overlaps `I_i(x,y)` / `E_n^i(x,y)` evaluated on")
    println(io, "the model's own FD wavefunctions (contact-distorted S waves, central P")
    println(io, "waves), with the paper's fitted exponents 0.7 / 0.5 and **no new fitted")
    println(io, "constants**. `rho -> pi gamma` is the paper's fit row for the 0.7")
    println(io, "exponent, so it doubles as the pipeline normalization check.")
    println(io)
    println(io, "Paper values in the open-flavor M1 block are marked **image-verified")
    println(io, "(page 24)**: each was read directly off the printed \"Predicted mu\"")
    println(io, "column of PDF page 24 (`paper/vision_ocr/page_images/page-024.png`,")
    println(io, "printed p. 212). The raw vision-OCR column (md lines 1124-1133) was")
    println(io, "displaced by one row against the decay labels; the printed page confirms")
    println(io, "the earlier one-row shift-correction exactly (K*+ +0.91, K*0 -1.20,")
    println(io, "D*+ -0.35, D*0 +1.78, F* -0.13, B* +1.37, B*0 -0.78, F_b* -0.55). No")
    println(io, "values changed by the crop audit. (Printed labels: the u-bbar row is")
    println(io, "B*- -> B- gamma; F_b* carries footnote f.) Open-flavor formula")
    println(io, "coefficients are rebuilt from quark charges via mu = e_q I_q - e_qbar I_qbar.")
    println(io)
    println(io, "## Magnetic-dipole moments (mu / mu_N)")
    println(io)
    println(io, "| Decay | Computed | Paper |")
    println(io, "|---|---|---|")
    for (name, val, target) in m1_rows
        @printf(io, "| %s | %+.3f | %s |\n", name, val, target)
    end
    println(io)
    println(io, "## E1 (and mixed) multipole amplitudes (MeV^(1/2))")
    println(io)
    println(io, "Photon momenta from measured 1984 masses; radial P waves are the")
    println(io, "central-channel solutions (paper's first-order treatment leaves the")
    println(io, "chi radial function J-independent, so rows differ only via q and the")
    println(io, "angular coefficient).")
    println(io)
    println(io, "| Decay | Computed | Paper |")
    println(io, "|---|---|---|")
    for (name, val, target) in e1_rows
        @printf(io, "| %s | %+.3f | %s |\n", name, val, target)
    end
    println(io)
    println(io, "## 2S -> chi_0 photon-momentum convention test")
    println(io)
    println(io, "The two `2S -> chi_0` E1 rows sit high with q from measured 1984")
    println(io, "masses. They carry the largest q of the 2S->chi block, so the amplitude")
    println(io, "(`~ q^{3/2}`, from `coeff(q) * sqrt(q)`) is the most q-sensitive row and")
    println(io, "the natural place to test the mass convention. Below, `q` is recomputed")
    println(io, "from the MODEL (GI predicted) parent/daughter masses -- eigenvalues of")
    println(io, "the same 2^3S_1 and central 1^3P_0 solves, no new constants -- and shown")
    println(io, "side by side with the measured-mass variant. Everything else (E_1 moment,")
    println(io, "coefficient) is held fixed.")
    println(io)
    println(io, "| Decay | q (meas) MeV | Computed (meas) | q (model) MeV | Computed (model) | Paper |")
    println(io, "|---|---|---|---|---|---|")
    for (name, target, q_meas, amp_meas, q_model, amp_model, mm) in chi0_variant_rows
        @printf(io, "| %s | %.1f | %+.3f | %.1f | %+.3f | %s |\n",
            name, q_meas * 1000, amp_meas, q_model * 1000, amp_model, target)
    end
    println(io)
    for (name, target, q_meas, amp_meas, q_model, amp_model, mm) in chi0_variant_rows
        p = parse(Float64, target)
        d_meas = 100 * reldev(amp_meas, p)
        d_model = 100 * reldev(amp_model, p)
        closer = d_model < d_meas ? "model masses" : "measured masses"
        @printf(io, "- `%s`: measured-mass q gives %.1f%% deviation, model-mass q (parent %.3f -> daughter %.3f GeV) gives %.1f%%. **%s closer.**\n",
            name, d_meas, mm[1], mm[2], d_model, closer)
    end
    println(io)
    println(io, "## Conventions used")
    println(io)
    println(io, "- `mu/mu_N = coefficient * I_i * M_N` (Table VI lists moments in units")
    println(io, "  of `e/2`; `M_N = ", M_N_GEV, " GeV`).")
    println(io, "- E1 amplitude = `coeff(q_GeV) * E_1^i * sqrt(alpha * q_MeV)`.")
    println(io, "- Mock masses `M~ = <E_1> + <E_2>` from the momentum-space waves.")
    println(io, "- Phase convention: `Phi(0) > 0` for every radial level (matches the")
    println(io, "  annihilation-block convention).")
end

println("wrote ", outpath)
println()
println("M1 rows:")
for (name, val, target) in m1_rows
    @printf("  %-46s %+8.3f   paper %s\n", name, val, target)
end
println()
println("E1 rows:")
for (name, val, target) in e1_rows
    @printf("  %-46s %+8.3f   paper %s\n", name, val, target)
end
println()
println("2S -> chi_0 q-convention test:")
for (name, target, q_meas, amp_meas, q_model, amp_model, mm) in chi0_variant_rows
    p = parse(Float64, target)
    @printf("  %-24s meas q=%.1f MeV -> %+.3f (%.1f%%) | model q=%.1f MeV -> %+.3f (%.1f%%) | paper %s\n",
        name, q_meas * 1000, amp_meas, 100 * reldev(amp_meas, p),
        q_model * 1000, amp_model, 100 * reldev(amp_model, p), target)
end
