# Invoked by `Pkg.test()` with the package environment already active.
# To run manually: `julia --project=. test/runtests.jl` from the repository root.
#
# The default suite is the quick one CI runs. Solver-convergence sweeps live in
# test/heavy/ and run only with GI_HEAVY_TESTS=true:
#   GI_HEAVY_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'

using Test
using FiniteDifferences
using LinearAlgebra
using QuadGK
using SpecialFunctions: erf
using GIModel

const TEST_ROOT = dirname(@__DIR__)
const root = TEST_ROOT

include("testutils.jl")

@testset "GIModel" begin
    include("parameters_and_core.jl")
    include("central_potentials.jl")
    include("fine_structure.jl")
    include("spin_kernels.jl")
    include("contact_operator.jl")
    include("quarks.jl")
    include("oscillator_basis.jl")
    include("radial_waves.jl")
    include("fixed_channel_solvers.jl")
    include("spectrum.jl")
    include("aqua.jl")
end

if lowercase(get(ENV, "GI_HEAVY_TESTS", "false")) in ("1", "true", "yes")
    @testset "GIModel heavy" begin
        include(joinpath("heavy", "radial_waves.jl"))
        include(joinpath("heavy", "fixed_channel_solvers.jl"))
    end
end
