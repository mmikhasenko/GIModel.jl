include(joinpath(@__DIR__, "common.jl"))

"""Tensor-mixed caption states, with signed amplitudes and squared probabilities.

Use the native fixed-channel HO waves and the full 3S + 3D block. Select
by overlap with the named precursor, never by proximity to the paper mass.
"""
function compute_caption_mixing(; solver=OscillatorSolver())
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    levels = [BasisState(n, L, 3, 1) for L in ("S", "D") for n in 1:3]
    lines = ["figure\tsector\tstate\tmass_GeV\tcomponent\tamplitude\tpaper_amplitude\tprobability"]
    cases = [("iii", :q, :q, "2^3S_1", Dict("2^3S_1"=>1.0, "1^3D_1"=>0.04)),
             ("iv", :q, :s, "2^3S_1", Dict("2^3S_1"=>1.0, "1^3D_1"=>0.04)),
             ("vi", :c, :c, "1^3D_1", Dict("1^3D_1"=>1.0, "1^3S_1"=>0.01,
                "2^3S_1"=>-0.03, "3^3S_1"=>-0.01))]
    for (figure, f1, f2, precursor, paper) in cases
        spec = compute_spectrum(params, Meson(mq, f1, f2); levels, solver)
        # The spectrum labels are mass-ranked; identify the physical eigenstate
        # by its dominant precursor instead of assuming that ranking survives.
        weights = [sum(abs2(c.coefficient) for c in physical_components(spec, state)
            if c.basis.label == precursor; init=0.0) for state in spec.states]
        state = spec.states[argmax(weights)]
        components = physical_components(spec, state)
        anchor = only(c.coefficient for c in components if c.basis.label == precursor)
        phase = anchor < 0 ? -1 : 1
        @assert isapprox(sum(abs2(c.coefficient) for c in components), 1; atol=1e-8)
        for c in components
            amplitude = phase * c.coefficient
            push!(lines, join((figure, "$f1 $f2", precursor, state.mass_GeV, c.basis.label,
                amplitude, get(paper, c.basis.label, ""), abs2(amplitude)), '\t'))
        end
    end
    return write_paper_table("caption_mixing.tsv", join(lines, '\n'))
end

abspath(PROGRAM_FILE) == (@__FILE__) && compute_caption_mixing()
