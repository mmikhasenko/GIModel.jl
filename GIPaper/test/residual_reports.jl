@testset "write_residual_report renders rows" begin
    reference = load_ref(reference_spectrum_path("charmonium"))
    rows = compare_reference(params, mq, reference[1:4]; ngrid = 100, rmax = 10.0)
    mktempdir() do d
        path = joinpath(d, "report.md")
        GIPaper.write_residual_report(path, "Residuals: test", rows)
        text = read(path, String)
        @test occursin("# Residuals: test", text)
        @test occursin("Contribution Breakdown", text)
        @test occursin("Spin-Averaged Diagnostics", text)
    end
end
