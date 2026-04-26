using Test

root = dirname(@__DIR__)
include(joinpath(root, "src", "GIModel", "GIModel.jl"))
using .GIModel

@testset "reference loading" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    @test params.masses["c"] ≈ 1.628
    @test params.masses["b"] ≈ 4.977
    @test params.b ≈ 0.18
    @test params.appendix_a_central == false
    @test params.epsilon_c ≈ -0.168
    @test params.fine_structure == true
    @test params.k_spin_orbit > 0.0

    ccbar = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_charmonium.csv"))
    bbbar = load_reference_spectrum(joinpath(root, "data", "reference_spectrum_bottomonium.csv"))
    @test length(ccbar) == 28
    @test length(bbbar) == 30
    @test ccbar[1].quark_content == "c cbar"
    @test ccbar[1].composition == "1^1S_0"
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

@testset "quark mass resolution" begin
    params = load_parameters(joinpath(root, "data", "parameters.provisional.toml"))
    m1, m2 = parse_quark_masses(params, "ccbar", "ignore")
    @test m1 ≈ m2 ≈ params.masses["c"]
    m1, m2 = parse_quark_masses(params, "charmonium", "c cbar")
    @test m1 ≈ m2 ≈ params.masses["c"]
    m1, m2 = parse_quark_masses(params, "charmed", "-c dbar; c ubar")
    @test m1 ≈ params.masses["c"] && m2 ≈ params.masses["d"]
    m1, m2 = parse_quark_masses(params, "charmed_strange", "c sbar")
    @test m1 ≈ params.masses["c"] && m2 ≈ params.masses["s"]
end

@testset "Appendix A 3D smearing (constant preserves norm)" begin
    h = 0.02
    n = 2000
    r = collect(h:h:(h * n))
    σ = 1.0
    n_tail = max(0, Int(ceil(8 * σ / h)))
    r_ext = n_tail > 0 ? vcat(r, collect((r[end] + h):h:(r[end] + n_tail * h))) : r
    w = GIModel.smear_3d_radial(ones(length(r_ext)), r_ext, σ)
    w = w[1:n]
    @test maximum(abs.(w .- 1.0)) < 0.01
    @test w[1] ≈ 1.0 atol = 0.01
    @test w[div(n, 2)] ≈ 1.0 atol = 0.01
end

@testset "triplet fine-structure angular factors" begin
    for L in 1:4
        js = collect((L - 1):(L + 1))
        weights = [2J + 1 for J in js]
        ldot = [GIModel.LdotS(L, 1, J) for J in js]
        tensor = [GIModel.tensor_triplet_LJ(L, J, 1) for J in js]
        @test sum(weights .* ldot) ≈ 0.0 atol = 1e-12
        @test sum(weights .* tensor) ≈ 0.0 atol = 1e-12
    end
    @test GIModel.tensor_triplet_LJ(1, 0, 1) ≈ -4.0
    @test GIModel.tensor_triplet_LJ(1, 1, 1) ≈ 2.0
    @test GIModel.tensor_triplet_LJ(1, 2, 1) ≈ -0.4
end
