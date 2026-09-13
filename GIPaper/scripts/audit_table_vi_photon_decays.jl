#!/usr/bin/env julia
# Canonical Table VI audit. Reference identities/kinematics belong to GIPaper;
# overlap kernels, multipole assembly and physical-state composition to GIModel.
using Pkg
Pkg.activate(dirname(@__DIR__))
using GIModel, GIPaper, Dates, Printf, Statistics

const ROOT = dirname(@__DIR__)
const G = GIModel
params, mq = load_parameters_and_quark_masses(default_parameters_path())
references = load_table_vi()
states = load_table_vi_states()
# Optional independent comparator uses the same state/operator pipeline.
solver_name = get(ENV, "GI_TABLE_VI_SOLVER", "ho")
solver_name in ("ho", "fd") || error("GI_TABLE_VI_SOLVER must be ho or fd")
solver = solver_name == "fd" ?
    FiniteDifferenceSolver(ngrid = 700, rmax = 24.0, nlevels_per_channel = 3) :
    OscillatorSolver(nlevels_per_channel = 3)
println("solving native Table VI sectors ...")
flush(stdout)

# Solve only the spectroscopic states in this table; no extra D/F mixing
# partners are introduced where Table III assumes an unlisted angle is zero.
flavor_pairs = sort!(unique([s.flavors for s in values(states)]); by = string)
spectra = Dict()
for flavors in flavor_pairs
    levels = unique([s.basis for s in values(states) if s.flavors == flavors])
    # Table III P1 includes 2cc even when only its allowed M1 moment is tabulated.
    spectra[flavors] = compute_spectrum(params, Meson(mq, flavors...);
        levels = levels, solver = solver)
    println("  solved ", flavors); flush(stdout)
end

# Table III: 1nn,1ss,1cc,1bb,2nn,2ss,2cc. The same solver-native
# components are then exposed through the final MixedSpectrum.
iso_flavors = [(:q, :q), (:s, :s), (:c, :c), (:b, :b)]
p1_basis = [
    BasisState(n, "S", 1, 0; flavors = flavors)
    for n in 1:2 for flavors in iso_flavors if !(n == 2 && flavors == (:b, :b))
]
println("composing four-flavor P1, vector and tensor isoscalars ..."); flush(stdout)
iso_spec = add_isoscalar_annihilation(params, [spectra[f] for f in iso_flavors];
    pseudoscalar = PaperP1Annihilation(), pseudoscalar_basis = p1_basis,
    amplitudes = Dict(("S", 3, 1) => params.annihilation.s1_A,
                      ("P", 3, 2) => params.annihilation.a_3p2),
    annihilation_radial_levels = (1,))

# Physical state lookup never reads or adjusts a mixing eigenvector.
function resolved(id)
    descriptor = states[id]
    spec = descriptor.mixed ? iso_spec : spectra[descriptor.flavors]
    return spec, spectrum_state(spec, descriptor.basis)
end
const ISOVECTORS = Set(["rho", "pi", "pi0", "A2", "A1", "B"])
mass(flavor) = mq[String(flavor)]
charge(flavor) = flavor in (:u, :c) ? 2 / 3 : -1 / 3

function kinematics(row)
    p, d = states[row.id[2]], states[row.id[3]]
    isnothing(p.mass_GeV) || isnothing(d.mass_GeV) ? nothing :
        photon_momentum(p.mass_GeV, d.mass_GeV)
end

# Cache representation-native momentum transforms/mean energies. No sampling
# of an HO state into a mesh is used for either overlap or annihilation.
momentum_cache = IdDict()
mom(w, L) = get!(momentum_cache, w) do
    G.momentum_wave(w, L)
end

function transition(kernel, parent_id, daughter_id)
    pspec, pstate = resolved(parent_id)
    dspec, dstate = resolved(daughter_id)
    # Compose independently when a physical isoscalar connects to an unmixed
    # isovector/open-flavor state. Shared coefficients enter exactly once.
    return physical_state_amplitude(pspec, pstate) do p
        physical_state_amplitude(dspec, dstate) do d
            kernel(p, d)
        end
    end
end

function m1_kernel(p, d, parent_id, daughter_id, q; recoil = false)
    p.basis.flavors == d.basis.flavors || return 0.0
    p.basis.L_label == d.basis.L_label == "S" || return 0.0
    p.basis.multiplicity != d.basis.multiplicity || return 0.0
    singlet, triplet = p.basis.multiplicity == 1 ? (p, d) : (d, p)
    f1, f2 = p.basis.flavors
    m1, m2 = mass(f1), mass(f2)
    if f1 == f2
        coefficient = neutral_m1_charge(f1;
            isovector_left = parent_id in ISOVECTORS,
            isovector_right = daughter_id in ISOVECTORS)
        return recoil ? m1_recoil_moment(
            singlet.wave, triplet.wave, m1, coefficient, q) :
            m1_transition_moment(mom(singlet.wave, 0), mom(triplet.wave, 0),
                m1, m2, [(coefficient, m1)])
    end
    recoil && error("no unequal-flavor recoil prescription in Table VI")
    return m1_transition_moment(mom(singlet.wave, 0), mom(triplet.wave, 0),
        m1, m2, [(charge(f1), m1), (charge(f2), m2)])
end

function multipole_kernel(p, d, parent_id, daughter_id, q, multipole)
    p.basis.flavors == d.basis.flavors || return 0.0
    Set([p.basis.L_label, d.basis.L_label]) == Set(["S", "P"]) || return 0.0
    sw, pw = p.basis.L_label == "S" ? (p, d) : (d, p)
    f1, f2 = sw.basis.flavors
    if (sw.basis.multiplicity, pw.basis.multiplicity) == (1, 3)
        expected = pw.basis.J == 2 ? :M2 : :E1
        multipole == expected || error("incorrect spin-flip multipole identity")
        # The charged A states carry e_u-e_d=1; K2 carries both unequal
        # emitting-quark terms. These are the printed spin-flip angular coefficients.
        terms = f1 == f2 ? [(1.0, mass(f1))] :
            [(charge(f1), mass(f1)), (-charge(f2), mass(f2))]
        return spin_flip_photon_amplitude(sw.wave, pw.wave, terms, pw.basis.J, q)
    end
    sw.basis.multiplicity == pw.basis.multiplicity || return 0.0
    f1 == f2 || error("unexpected open-flavor E1 row")
    c = neutral_m1_charge(f1; isovector_left = parent_id in ISOVECTORS,
        isovector_right = daughter_id in ISOVECTORS)
    angular = e1_angular_coefficient(pw.basis.J;
        singlet = sw.basis.multiplicity == 1, parent_is_S = p.basis.L_label == "S")
    return e1_transition_amplitude(sw.wave, mom(sw.wave, 0),
        pw.wave, mom(pw.wave, 1), mass(f1), q -> c * angular * q,
        1.0, 0.0; q = q)
end

function evaluate(row)
    multipole, parent, daughter = row.id
    q = kinematics(row)
    footnotes = split(row.footnotes, ",")
    if isnothing(q) && (multipole != :M1 || "c" in footnotes || "g" in footnotes)
        value, supplementary = NaN, 0.0
    elseif multipole == :M1
        value = transition(parent, daughter) do p, d
            m1_kernel(p, d, parent, daughter, q; recoil = "c" in footnotes)
        end
        # Explicit paper input, not fitted here: page 26 footnote a attributes
        # this additive moment to pi0-eta mixing.
        supplementary = "a" in footnotes ? 0.01 : 0.0
        value += supplementary
    else
        isnothing(q) && error("missing photon kinematics for $(row.decay)")
        q >= 0 || error("closed multipole channel $(row.decay)")
        value = transition(parent, daughter) do p, d
            multipole_kernel(p, d, parent, daughter, q, multipole)
        end
        supplementary = 0.0
    end
    if "g" in footnotes && !isnothing(q)
        isnothing(q) && error("missing footnote-g kinematics")
        value *= photon_recoil_form_factor(q)
    end
    p, d = states[parent], states[daughter]
    historical_q = isnothing(p.historical_mass_GeV) || isnothing(d.historical_mass_GeV) ? nothing :
        photon_momentum(p.historical_mass_GeV, d.historical_mass_GeV)
    return (reference = row, computed = value, q_GeV = q,
        q_display_GeV = something(q, NaN), q_historical_GeV = historical_q,
        q_source = isnothing(q) ? "experimental assignment unavailable" : "PDG 2026",
        supplementary = supplementary)
end

println("evaluating all canonical rows ..."); flush(stdout)
results = evaluate.(references)
@assert length(results) == 79
@assert Set(r.reference.id for r in results) == Set(r.id for r in references)
@assert all(isfinite(r.computed) || isnothing(r.q_GeV) for r in results)
@assert all(isnothing(r.q_GeV) || r.q_GeV >= 0 for r in results)

@assert all(!iszero(r.computed) || r.reference.approximate for r in results)

output = joinpath(ROOT, "docs", "residual_reports",
    solver isa OscillatorSolver ? "table_vi_photon_decays.md" : "table_vi_photon_decays_fd.md")
# Machine-readable regression/coverage artifact: one result per canonical row.
csvpath = replace(output, ".md" => ".csv")
open(csvpath, "w") do io
    println(io, "multipole,parent,daughter,computed,paper,parenthetical,approximate,q_GeV,supplementary_mu,q_display_GeV,q_source,q_historical_GeV")
    for r in results
        ref = r.reference
        println(io, join((ref.id..., r.computed, ref.predicted,
            ref.parenthetical, ref.approximate, something(r.q_GeV, ""), r.supplementary, r.q_display_GeV, r.q_source, something(r.q_historical_GeV, "")), ","))
    end
end
open(output, "w") do io
    println(io, "# Table VI photon-decay audit")
    println(io, "\nGenerated by `scripts/audit_table_vi_photon_decays.jl` on ", Dates.today(), ".")
    println(io, "\nCoverage: **79/79 canonical rows** (42 M1, 35 E1, 2 M2); no omitted rows.")
    println(io, "\n", numerics_provenance(solver))
    println(io, """

    Every radial state uses the shared fixed-channel solver, including the
    spin-dependent P waves. P1 pseudoscalar mixing contains 1nn/1ss/1cc/1bb
    and 2nn/2ss/2cc; the n=1 vector and tensor blocks contain all four flavors.
    The final MixedSpectrum owns each signed composition. Unlisted higher
    radial mixing and D/F partners are not introduced. Pure-flavor charges
    are composed exactly once: the old ground-state eta calculation counted
    the perfect-mixing 1/sqrt(2) factor a second time.

    The paper's exponents 0.7/0.5 and Table III annihilation parameters are
    unchanged. Footnote c retains the E2 recoil correction; footnote g applies
    exp(-q²/(16β²)), β=0.40 GeV. Footnote a's additive +0.01 μN contribution
    to phi -> pi gamma is a **paper-supplied pi0-eta mixing input**.
    Parenthesized paper predictions are order-of-magnitude estimates, not
    precision targets. Approximately zero is preserved as a qualifier.

    M1 values are μ/μN; E1/M2 values are MeV^(1/2), whose squares are widths
    in MeV. Signs are shown without row-dependent rephasing or adjustments.
    Mixed states retain a positive overlap with their assigned unmixed state,
    using the same ascending-mass assignment as Spectrum. This shared phase
    convention avoids reversing heavy states because of a tiny negative nn
    admixture. No per-row sign adjustment is made.
    """)
    for multipole in (:M1, :E1, :M2)
        block = filter(r -> first(r.reference.id) == multipole, results)
        exact = filter(r -> isfinite(r.computed) && !r.reference.parenthetical && !r.reference.approximate &&
            r.reference.predicted != 0, block)
        signs = count(r -> sign(r.computed) == sign(r.reference.predicted), exact)
        ratios = [abs(r.computed / r.reference.predicted) for r in exact]
        println(io, "\n## ", multipole)
        @printf(io, "\nNon-parenthesized targets: %d/%d signs; median |model/paper| %.3f.\n",
            signs, length(exact), median(ratios))
        println(io, "\n| Decay | Computed | Paper | q MeV | q reference MeV | src |")
        println(io, "|---|---:|---:|---:|---:|---|")
        for r in block
            qtext = @sprintf("%.1f", 1000r.q_display_GeV)
            qref = isnothing(r.q_historical_GeV) ? "—" : @sprintf("%.1f", 1000r.q_historical_GeV)
            @printf(io, "| %s | %+.6g | %s | %s | %s | %s |\n", r.reference.decay,
                r.computed, r.reference.predicted_text, qtext, qref, r.q_source)
        end
    end
    println(io, "\n## Largest magnitude residuals")
    ordinary = filter(r -> isfinite(r.computed) && !r.reference.parenthetical && r.reference.predicted != 0, results)
    ranked = sort(ordinary; by = r -> abs(abs(r.computed / r.reference.predicted) - 1),
        rev = true)
    for r in ranked[1:min(8, length(ranked))]
        @printf(io, "\n- %s: |model/paper| %.3f.", r.reference.decay,
            abs(r.computed / r.reference.predicted))
    end
    println(io, "\n\n## Kinematics provenance")
    println(io, "\nAll operator and displayed momenta use the same PDG 2026 registry. No model eigenvalue is used for kinematics.")
    println(io, "Historical reference q values are reconstructed from legacy audit inputs, some of which were GI model estimates; they are not printed GI momenta.")
    println(io, "Missing experimental assignments leave q and q-dependent amplitudes unavailable. Momentum-independent M1 moments remain calculable.")
    println(io, "\n## Encoding corrections")
    println(io, "\nThe 3S -> 1S bottomonium target is -0.004; +0.007 belongs to 3S -> 2S.")
    println(io, "The A2 -> pi M2 denominator is sqrt(60) m_u, not sqrt(60 m_u).")
    println(io, "A1 -> pi is spin-flip E1, not M2: J=1 -> J=0 allows a dipole photon.")
    println(io, "Canonical identities, rather than CSV row positions, select every target.")
end
println("wrote ", output)
println("wrote ", csvpath)
