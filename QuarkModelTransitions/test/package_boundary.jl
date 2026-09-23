@testset "Package boundary" begin
    transition_names = (
        :TransitionOperator,
        :StrongDecayOperator,
        :PseudoscalarEmission,
        :PhysicalState,
        :TwoMesonChannel,
        :TransitionAmplitude,
        :matrix_element,
        :decay_width,
    )
    @test all(name -> Base.isexported(QuarkModelTransitions, name), transition_names)
    @test all(name -> !Base.isexported(GIModel, name), transition_names)

    gi_project = read(joinpath(REPOSITORY_ROOT, "Project.toml"), String)
    @test !occursin("QuarkModelTransitions", gi_project)
end
