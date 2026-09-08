@testset "Table V paper adapter" begin
    amplitude_row = (
        parent = "rho",
        daughter1 = "omega",
        daughter2 = "pi",
        coefficient = 1.25,
        amp_class = "P",
        qbar_power = 1,
        decay = "rho -> omega pi",
        section = "light",
    )
    ignored_row = merge(amplitude_row, (amp_class = "mixing_only",))

    channels = GIPaper.load_table_v((amplitude_row, ignored_row))
    @test length(channels) == 1
    @test only(channels) isa DecayChannel
    @test only(channels).class == :P
    @test only(channels).label == "rho -> omega pi"

    heavy_row = merge(amplitude_row, (amp_class = "Ac", section = "charmed"))
    @test_throws ArgumentError GIPaper.load_table_v((heavy_row,))
    heavy = only(GIPaper.load_table_v((heavy_row,); heavy_fraction_for = _ -> 0.8))
    @test heavy.class == :A_c
    @test heavy.heavy_fraction == 0.8
end
