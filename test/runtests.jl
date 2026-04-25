using Test

root = dirname(@__DIR__)
include(joinpath(root, "src", "GIModel", "GIModel.jl"))
using .GIModel

@testset "reference loading" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    @test params.masses["c"] ≈ 1.628
    @test params.masses["b"] ≈ 4.977
    @test params.b ≈ 0.18

    ccbar = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    bbbar = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_bottomonium.csv"))
    @test length(ccbar) == 28
    @test length(bbbar) == 30
end

@testset "baseline solver shape" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    cc = solve_sector(params, "c"; maxn = 4, ngrid = 250, rmax = 20.0)
    bb = solve_sector(params, "b"; maxn = 4, ngrid = 250, rmax = 16.0)

    @test cc[(1, "S")] < cc[(2, "S")] < cc[(3, "S")]
    @test bb[(1, "S")] < bb[(2, "S")] < bb[(3, "S")]
    @test cc[(1, "S")] < cc[(1, "P")] < cc[(1, "D")]
    @test bb[(1, "S")] < bb[(1, "P")] < bb[(1, "D")]
end
