# Shared, dependency-light definitions for the formal Table-IV A-class survey.
# This file is included by the runner and its tests; it does not change Table V.
# Archived study-branch API: coupling_coefficients/topology/integral_cache are
# not part of the current native matrix_element API. See the study README.
module AClassKernels
using GIModel, GIModel.QuarkModelTransitions, LinearAlgebra
using GIModel.QuarkModelTransitions: CMKinematics
export TRANSITIONS, flavor_routes, coefficients, compare_couplings, threshold_limit

const TRANSITIONS = (
    (family="A", id="A1", L="S", multiplicity=3, J=1, daughter_L="S", daughter_spin=1, daughter_J=0, sho_hg=0.25, relative_L=1),
    (family="A", id="A2", L="P", multiplicity=3, J=2, daughter_L="S", daughter_spin=1, daughter_J=0, sho_hg=0.25, relative_L=2),
    (family="A", id="A3", L="P", multiplicity=3, J=2, daughter_L="S", daughter_spin=3, daughter_J=1, sho_hg=0.25, relative_L=2),
    (family="A", id="A4", L="P", multiplicity=3, J=1, daughter_L="S", daughter_spin=3, daughter_J=1, sho_hg=0.25, relative_L=2),
    (family="A", id="A5", L="P", multiplicity=1, J=1, daughter_L="S", daughter_spin=3, daughter_J=1, sho_hg=0.25, relative_L=2),
    (family="A", id="A6", L="D", multiplicity=3, J=3, daughter_L="S", daughter_spin=1, daughter_J=0, sho_hg=0.25, relative_L=3),
    (family="A", id="A7", L="D", multiplicity=3, J=3, daughter_L="S", daughter_spin=3, daughter_J=1, sho_hg=0.25, relative_L=3),
    (family="Aprime", id="A8", L="P", multiplicity=1, J=1, daughter_L="P", daughter_spin=3, daughter_J=0, sho_hg=-0.25, relative_L=1),
    (family="Adoubleprime", id="A9", L="D", multiplicity=3, J=3, daughter_L="P", daughter_spin=1, daughter_J=1, sho_hg=0.125, relative_L=2),
)

"""Finite census: n=1 basis states, all u/d/s emissions, both charge conjugates.
Isoscalar rows are resolved uu/dd/ss basis components, not mixed physical poles.
Heavy quarks are spectators; quarkonia and Bc are absent by construction.
"""
function flavor_routes()
    pairs = [(sector="isovector", pair=p) for p in ((:u,:d), (:d,:u))]
    append!(pairs, [(sector="isoscalar", pair=p) for p in ((:u,:u),(:d,:d),(:s,:s))])
    append!(pairs, [(sector="strange", pair=p) for p in ((:u,:s),(:d,:s),(:s,:u),(:s,:d))])
    for (sector, heavy) in (("charmed", :c), ("bottom-flavored", :b)), light in (:u,:d,:s)
        push!(pairs, (sector=sector, pair=(heavy,light)))
        push!(pairs, (sector=sector, pair=(light,heavy)))
    end
    rows = NamedTuple[]
    for entry in pairs, topology in (:quark,:antiquark)
        e = topology == :quark ? 1 : 2
        entry.pair[e] in (:u,:d,:s) || continue
        for produced in (:u,:d,:s)
            daughter = e == 1 ? (produced,entry.pair[2]) : (entry.pair[1],produced)
            emitted = e == 1 ? (entry.pair[1],produced) : (produced,entry.pair[2])
            push!(rows, (;entry.sector, parent=entry.pair, daughter, emitted, topology,
                emitter=entry.pair[e], spectator=entry.pair[3-e],
                light_constituent=entry.sector in ("charmed","bottom-flavored") ?
                    (entry.pair[e] == :s ? "strange" : "nonstrange") : "light"))
        end
    end
    return rows
end

function coefficients(transition, route, parent_wave, daughter_wave, masses, q;
                      parent_mass=2.0, daughter_mass=1.0, beta=0.4, integral_cache=nothing)
    b = BasisState(1, transition.L, transition.multiplicity, transition.J; flavors=route.parent)
    d = BasisState(1, transition.daughter_L, transition.daughter_spin, transition.daughter_J; flavors=route.daughter)
    p = BasisState(1, "S", 1, 0; flavors=route.emitted)
    state(b, w, m, tag) = PhysicalState(tag * b.label, m, [(basis=b, coefficient=1.0, wave=w)];
                                      provenance=(source=:A_class_unmixed_basis,))
    initial = state(b, parent_wave, parent_mass, "parent:")
    final = TwoMesonChannel(state(d, daughter_wave, daughter_mass, "daughter:"),
                           state(p, OscillatorWave(0,beta,[1.0]), 0.14, "field:"))
    c = coupling_coefficients(final, PseudoscalarEmission(0,0,masses), initial;
                              kinematics=CMKinematics(q), topology=route.topology, integral_cache)
    pw = PartialWave(transition.relative_L, transition.daughter_J)
    return (g=Dict(c.g.partial_wave_amplitudes)[pw], h=Dict(c.h.partial_wave_amplitudes)[pw])
end

"""Complex ratios plus a phase-invariant fallback. A tiny vector has no direction.
The scalar gate is relative to the coupling-vector norm (after q^L reduction).
"""
function compare_couplings(native, sho; tolerance=1e-8)
    v, s = ComplexF64[native.g,native.h], ComplexF64[sho.g,sho.h]
    nv, ns = norm(v), norm(s)
    vector_valid = min(nv,ns) > 1e-12
    stable = vector_valid && minimum(abs.(v)) > tolerance*nv && minimum(abs.(s)) > tolerance*ns
    overlap = vector_valid ? dot(s/ns,v/nv) : 0im
    phase = abs(overlap) > tolerance ? conj(overlap)/abs(overlap) : 1.0+0im
    distance = vector_valid ? norm(phase*v/nv-s/ns) : NaN
    return (Rg=stable ? v[1]/s[1] : ComplexF64(NaN),
            Rh=stable ? v[2]/s[2] : ComplexF64(NaN),
            U_hg=stable ? (v[2]/v[1])/(s[2]/s[1]) : ComplexF64(NaN),
            scalar_stable=stable, vector_valid, vector_distance=distance)
end

"""Richardson q² extrapolation with two successive independent estimates.
Input samples remain in the raw trace. Certification requires stability of BOTH
reduced coefficients, not merely a cancellation in their ratio.
"""
function threshold_limit(q, values, L; tolerance=2e-3)
    length(q) >= 3 || throw(ArgumentError("at least three nonzero q samples required"))
    all(>(0),q) && all(<(0),diff(q)) || throw(ArgumentError("q must decrease and stay positive"))
    reduced = [(g=v.g/x^L,h=v.h/x^L) for (x,v) in zip(q,values)]
    extrap(i, piece) = begin
        a,b = q[i-1]^2,q[i]^2
        (a*getproperty(reduced[i],piece)-b*getproperty(reduced[i-1],piece))/(a-b)
    end
    n = length(q)
    last_estimate = (g=extrap(n,:g),h=extrap(n,:h))
    previous = (g=extrap(n-1,:g),h=extrap(n-1,:h))
    error = maximum(abs(getproperty(last_estimate,p)-getproperty(previous,p))/
                    max(abs(getproperty(last_estimate,p)),1e-12) for p in (:g,:h))
    return (value=last_estimate, error=error, certified=error<tolerance)
end
end
