using Test
using GIModel
using QuarkModelTransitions

const REPOSITORY_ROOT = dirname(dirname(@__DIR__))
const root = REPOSITORY_ROOT

@testset "QuarkModelTransitions" begin
    include("package_boundary.jl")
    include("strong_decays.jl")
    include("transition_amplitudes.jl")
end
