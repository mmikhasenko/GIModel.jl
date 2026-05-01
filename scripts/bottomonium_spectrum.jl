#!/usr/bin/env julia
# Prints GI-model bottomonium masses (GeV) for states in Fig. 8 reference list.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DataFrames

root = dirname(@__DIR__)
using GIModel

params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
ref = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_bottomonium.csv"))
computed = compute_sector(params, ref, "b"; kinetic = :relativistic)
rows = compare(
    computed,
    ref;
    contact_hyperfine = true,
    use_fine_structure = params.fine_structure,
)

df = DataFrame(rows)
df = transform(
    df,
    [:n, :multiplicity, :L, :J] =>
        ByRow((n, mult, L, J) -> string(n, '^', mult, L, '_', J)) => :state,
    :predicted_GeV => ByRow(x -> round(x; digits = 3)) => :model_GeV,
    :reference_GeV => ByRow(x -> round(x; digits = 3)) => :ref_GeV,
    :residual_MeV => ByRow(x -> round(x; digits = 1)) => :delta_MeV,
)
df = select(df, :state, :model_GeV, :ref_GeV, :delta_MeV)

println("Bottomonium (predicted masses, GeV)\n")
show(stdout, df; allrows = true, show_row_number = false)
println()
