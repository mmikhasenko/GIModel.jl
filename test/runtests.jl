# Invoked by `Pkg.test()` with the package environment already active.
# To run manually: `julia --project=. test/runtests.jl` from the repository root.

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
    include("strong_decays.jl")
    include("transition_amplitudes.jl")
    include("quarks.jl")
    include("oscillator_basis.jl")
    include("radial_waves.jl")
    include("fixed_channel_solvers.jl")
    include("spectrum.jl")
    include("observables.jl")
    include("radiative_decays.jl")
end
