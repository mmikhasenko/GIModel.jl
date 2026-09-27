# Invoked by `Pkg.test()` with the GIPaper environment active.
# To run manually: `julia --project=GIPaper GIPaper/test/runtests.jl` from the
# repository root.

using Test
using CSV
using TOML
using GIModel
using QuarkModelTransitions
# Internal kernels for reference calculations.
using QuarkModelTransitions: DecayChannel
using GIPaper

const GIPAPER_ROOT = dirname(@__DIR__)
const REPOSITORY_ROOT = dirname(GIPAPER_ROOT)

const load_ref = GIPaper.load_reference_spectrum

const params, mq = load_parameters_and_quark_masses(GIPaper.model_parameters_path())

@testset "GIPaper" begin
    include("package_boundary.jl")
    include("data_validation.jl")
    include("reference_loading.jl")
    include("table_v.jl")
    include("table_vi.jl")
    include("mass_inputs.jl")
    include("comparison.jl")
    include("annihilation.jl")
    include("residual_reports.jl")
end
