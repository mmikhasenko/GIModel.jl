using Test
using LinearAlgebra
using GIModel
using QuarkModelTransitions

const REPOSITORY_ROOT = dirname(dirname(@__DIR__))
const root = REPOSITORY_ROOT

@testset "QuarkModelTransitions" begin
    include("package_boundary.jl")
    include("flavor_algebra.jl")
    include("spin_algebra.jl")
    include("algebraic_decomposition.jl")
    include("pseudoscalar_emission.jl")
    include("observables.jl")
    include("radiative_decays.jl")
    include("strong_decays.jl")
    include("transition_amplitudes.jl")
end
