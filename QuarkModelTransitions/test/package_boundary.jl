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
        :m1_transition_moment, :e1_transition_amplitude, :m1_radiative_width,
        :gluonic_annihilation_width, :leptonic_decay_factor, :two_photon_amplitude,
        :charge_radius_squared, :ALPHA_EM, :G_FERMI_GEV,
        :ELECTROMAGNETIC_DEFAULTS, :STRONG_DECAY_DEFAULTS,
    )
    @test all(name -> Base.isexported(QuarkModelTransitions, name), transition_names)
    @test all(name -> !Base.isexported(GIModel, name), transition_names)

    @test all(name -> !isdefined(GIModel, name), transition_names)
    for name in (:radial_overlap, :momentum_overlap, :momentum_functional,
                 :compute_spectrum, :add_isoscalar_annihilation)
        @test Base.isexported(GIModel, name)
    end

    gi_project = read(joinpath(REPOSITORY_ROOT, "Project.toml"), String)
    @test !occursin("QuarkModelTransitions", gi_project)
end
