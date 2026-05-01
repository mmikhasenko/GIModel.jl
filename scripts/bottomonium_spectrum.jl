#!/usr/bin/env julia
# Prints GI-model bottomonium masses (GeV) for states in Fig. 8 reference list.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf

root = dirname(@__DIR__)
using GIModel

params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
ref = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_bottomonium.csv"))
rows = compare_sector(
    params,
    ref,
    "b";
    kinetic = :relativistic,
    contact_hyperfine = true,
    use_fine_structure = params.fine_structure,
)

println("Bottomonium (predicted masses, GeV)\n")
@printf("%-12s  %10s  %10s  %11s\n", "state", "model", "ref", "Δ MeV")
@printf("%s\n", repeat("-", 48))
for row in rows
    lab = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
    @printf(
        "%-12s  %10.3f  %10.3f  %+10.1f\n",
        lab,
        row.predicted_GeV,
        row.reference_GeV,
        row.residual_MeV,
    )
end
