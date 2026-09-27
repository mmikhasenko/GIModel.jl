include(joinpath(@__DIR__, "common.jl"))
using LinearAlgebra

const HBARC_GEV_FM = 0.1973269804

"""
Solve the S-wave nonrelativistic Coulomb-plus-linear problem used to make the
motivational estimate in GI Table I,

    H = p²/m_Q - (4/3)alpha_s/r + b*r.

The returned confinement share is the change of the linear-potential
expectation value between 2S and 1S divided by the total level splitting.  This
is the Hellmann--Feynman contribution at the quoted value of `b`; GI do not
state an operational definition for their rounded percentage column.
"""
function table_i_confinement_share(mass, alpha, b; ngrid=6000)
    reduced_mass = mass / 2
    coulomb_strength = (4 / 3) * alpha
    a0 = inv(reduced_mass * coulomb_strength)
    linear_length = cbrt(inv(mass * b))
    rmax = 18 * max(a0, linear_length)
    h = rmax / (ngrid + 1)
    r = collect(1:ngrid) .* h

    diagonal = @. inv(reduced_mass * h^2) - coulomb_strength / r + b * r
    offdiagonal = fill(-inv(2 * reduced_mass * h^2), ngrid - 1)
    levels = eigen(SymTridiagonal(diagonal, offdiagonal), 1:2)
    energies = levels.values
    linear_expectations = [sum(abs2.(levels.vectors[:, n]) .* (b .* r)) for n in 1:2]
    splitting = energies[2] - energies[1]
    linear_splitting = linear_expectations[2] - linear_expectations[1]
    return (; splitting, linear_splitting, share_percent=100 * linear_splitting / splitting)
end

function compute_table_i()
    source = joinpath(GIPAPER_DIR, "data", "raw", "digitized_tables",
        "table_i_confinement", "table_i_confinement.csv")
    rows = filter(!isempty, readlines(source)[2:end])
    b = GIModel.load_parameters(GIModel.default_parameters_path()).potential.b
    lines = ["mass_GeV\talpha_s\tb_GeV2\tcomputed_a0_fm\tgi_a0_fm\t" *
             "computed_splitting_GeV\tcomputed_linear_splitting_GeV\t" *
             "computed_confinement_percent\tgi_confinement_percent"]
    for row in rows
        fields = split(row, ','; limit=9)
        mass = parse(Float64, fields[5])
        alpha = parse(Float64, fields[6])
        gi_radius = parse(Float64, fields[7])
        share = parse(Int, fields[8])
        computed = HBARC_GEV_FM / ((2 / 3) * alpha * mass)
        result = table_i_confinement_share(mass, alpha, b)
        push!(lines, join((mass, alpha, b, computed, gi_radius, result.splitting,
            result.linear_splitting, result.share_percent, share), '\t'))
    end
    return write_paper_table("table_i.tsv", join(lines, '\n'))
end

abspath(PROGRAM_FILE) == (@__FILE__) && compute_table_i()
