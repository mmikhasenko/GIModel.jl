@testset "Linear couplings and emitter topology" begin
    masses = QuarkMassTable("u"=>0.22,"d"=>0.22,"s"=>0.419,"c"=>1.628)
    state(L,S,J,pair,tag) = PhysicalState(tag,2.0,[(basis=BasisState(1,L,S,J;flavors=pair),
        coefficient=1.0,wave=OscillatorWave(L=="S" ? 0 : 1,0.4,[1.0]))])
    cache=Dict()
    for pair in ((:u,:u),(:c,:u),(:u,:c)), q in (0.23,0.1+0.03im)
        initial=state("P",3,1,pair,"parent")
        final=TwoMesonChannel(state("S",3,1,pair,"surviving"),state("S",1,0,(:u,:u),"field"))
        op=PseudoscalarEmission(-0.73,1.29,masses)
        cached=coupling_coefficients(final,op,initial;kinematics=CMKinematics(q),integral_cache=cache)
        allc=coupling_coefficients(final,op,initial;kinematics=CMKinematics(q))
        @test cached.g.partial_wave_amplitudes == allc.g.partial_wave_amplitudes
        @test cached.h.partial_wave_amplitudes == allc.h.partial_wave_amplitudes
        old_length=length(cache)
        again=coupling_coefficients(final,op,initial;kinematics=CMKinematics(q),integral_cache=cache)
        @test length(cache)==old_length
        @test again.g.partial_wave_amplitudes==allc.g.partial_wave_amplitudes
        shifted=matrix_element(final,PseudoscalarEmission(2,-3,masses),initial;kinematics=CMKinematics(q),integral_cache=cache)
        uncached=matrix_element(final,PseudoscalarEmission(2,-3,masses),initial;kinematics=CMKinematics(q))
        @test shifted.partial_wave_amplitudes==uncached.partial_wave_amplitudes
        for topology in (:both,:quark,:antiquark)
            c=coupling_coefficients(final,op,initial;kinematics=CMKinematics(q),topology)
            a=matrix_element(final,op,initial;kinematics=CMKinematics(q),topology)
            for field in (:helicity,:partial_wave_amplitudes), i in eachindex(getproperty(a,field))
                @test last(getproperty(a,field)[i]) ≈ op.g*last(getproperty(c.g,field)[i])+op.h*last(getproperty(c.h,field)[i]) atol=1e-12
            end
        end
        cq=coupling_coefficients(final,op,initial;kinematics=CMKinematics(q),topology=:quark)
        ca=coupling_coefficients(final,op,initial;kinematics=CMKinematics(q),topology=:antiquark)
        for piece in (:g,:h), i in eachindex(allc.g.partial_wave_amplitudes)
            value(c)=last(getproperty(c,piece).partial_wave_amplitudes[i])
            @test value(allc) ≈ value(cq)+value(ca) atol=1e-12
        end
        @test_throws ArgumentError coupling_coefficients(final,op,initial;topology=:invalid)
    end
end
