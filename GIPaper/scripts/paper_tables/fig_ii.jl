include(joinpath(@__DIR__, "common.jl"))

function matched_lambdas()
    masses = (0.419, 1.628, 4.977)
    lambdas = [0.2]
    nf = 2
    for mass in masses
        q = 2mass
        logratio = ((33 - 2nf) / (33 - 2(nf + 1))) * log(q^2 / lambdas[end]^2)
        push!(lambdas, q * exp(-logratio / 2))
        nf += 1
    end
    return masses, lambdas
end

function alpha_lo_qcd(q, masses, lambdas)
    nf = 2 + count(mass -> q > 2mass, masses)
    denominator = (33 - 2nf) * log(q^2 / lambdas[nf - 1]^2)
    return denominator <= 0 ? NaN : 12pi / denominator
end

function compute_fig_ii()
    masses, lambdas = matched_lambdas()
    lines = ["Q_GeV\talpha_GIModel\talpha_LO_QCD"]
    for q in exp.(range(log(0.45), log(20.0); length=900))
        push!(lines, join((q, GIModel.alpha_s_q(q), alpha_lo_qcd(q, masses, lambdas)), '\t'))
    end
    return write_paper_table("fig_ii.tsv", join(lines, '\n'))
end

abspath(PROGRAM_FILE) == (@__FILE__) && compute_fig_ii()
