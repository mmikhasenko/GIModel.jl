#!/usr/bin/env julia
# Canonical Table VI audit. Reference identities/kinematics belong to GIPaper;
# operators belong to QuarkModelTransitions; states and waves belong to GIModel.
mkpath(joinpath(dirname(@__DIR__), "reports"))

using Pkg
Pkg.activate(dirname(@__DIR__))
using GIModel, GIModel.QuarkModelTransitions, GIPaper, Dates, Printf, Statistics
using GIModel.QuarkModelTransitions: photon_momentum

const ROOT = dirname(@__DIR__)
const G = GIModel
const QMT = QuarkModelTransitions
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
const TABLE_POLICY = load_table_policy()["table_vi"]
const ISOVECTORS = Set(TABLE_POLICY["isovectors"])

function resolved_physical(id; charged = false)
    spec, state = resolved(id)
    native = physical_state(spec, state)
    components = []
    for c in native.components
        b = c.basis
        if b.flavors == (:q, :q)
            flavors = charged ? ((:u, :d),) : ((:u, :u), (:d, :d))
            for (i, pair) in enumerate(flavors)
                phase = !charged && id in ISOVECTORS && i == 2 ? -1 : 1
                coefficient = c.coefficient * phase / (charged ? 1 : sqrt(2))
                push!(components, (basis = BasisState(b.n, b.L_label, b.multiplicity, b.J;
                    label = b.label, flavors = pair), coefficient, wave = c.wave))
            end
        else
            push!(components, c)
        end
    end
    return PhysicalState(native.label, native.mass_GeV, components;
                         provenance = (source = :TableVI_explicit_isospin,))
end

function kinematics(row)
    p, d = states[row.id[2]], states[row.id[3]]
    isnothing(p.mass_GeV) || isnothing(d.mass_GeV) ? nothing :
        photon_momentum(p.mass_GeV, d.mass_GeV)
end

# Table VI footnote c requests the paper's relative q^2 recoil correction.
# This lookup is audit policy, not transition-package behavior.
table_vi_recoil_order(footnotes) = "c" in footnotes ? 2 : 0

function evaluate(row)
    multipole, parent, daughter = row.id
    q = kinematics(row)
    footnotes = split(row.footnotes, ",")
    if isnothing(q) && (multipole != :M1 || "c" in footnotes || "g" in footnotes)
        value, supplementary = NaN, 0.0
    else
        # The A1/A2 -> pi spin-flip rows are charged transitions. Their old
        # single-emitter coefficient was 1; explicit u dbar gives e_u-e_d.
        pdesc, ddesc = states[parent], states[daughter]
        charged = pdesc.flavors == ddesc.flavors == (:q, :q) &&
                  pdesc.basis.multiplicity != ddesc.basis.multiplicity && multipole != :M1
        parent_state = resolved_physical(parent; charged)
        daughter_state = resolved_physical(daughter; charged)
        operator = PhotonEmission(
            mq;
            recoil_order = table_vi_recoil_order(footnotes),
            recoil_form_factor = "g" in footnotes,
        )
        amplitude = matrix_element(daughter_state, operator, parent_state)
        # A plain M1 moment has no q dependence; missing target masses leave it
        # at the model value. Other rows require the explicit comparison factor.
        correction = isnothing(q) ? 1.0 : mass_correction_factor(
            daughter_state, operator, parent_state; target_momentum=q,
        )
        amplitude.multipole == multipole || error(
            "Table VI labels $(parent) -> $(daughter) as $multipole, " *
            "but spectroscopy selects $(amplitude.multipole)",
        )
        value = real(correction * amplitude.value)
        # Explicit paper input, not fitted here: page 26 footnote a attributes
        # this additive moment to pi0-eta mixing.
        supplementary = multipole == :M1 &&
            TABLE_POLICY["supplementary_footnote"] in footnotes ?
            TABLE_POLICY["supplementary_mu_N"] : 0.0
        value += supplementary
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

output = joinpath(ROOT, "reports",
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
    println(io, "\nGenerated by `checks/audit_table_vi_photon_decays.jl` on ", Dates.today(), ".")
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
