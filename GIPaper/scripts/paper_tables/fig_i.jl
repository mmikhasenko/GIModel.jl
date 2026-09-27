include(joinpath(@__DIR__, "common.jl"))

"""Evaluate the coupling entering the Fig. 1 annihilation graph at explicit scales.

M is an evaluation-scale input in GeV, not a predicted meson mass. The paper
writes alpha_s(M²); the GIModel API takes M, and squares it internally.
"""
function compute_fig_i()
    lines = ["M_GeV\talpha_s_M"]
    for M in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0)
        push!(lines, join((M, GIModel.alpha_s_q(M)), '\t'))
    end
    return write_paper_table("fig_i.tsv", join(lines, '\n'))
end
abspath(PROGRAM_FILE) == (@__FILE__) && compute_fig_i()
