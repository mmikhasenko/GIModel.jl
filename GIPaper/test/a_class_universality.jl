include(joinpath(@__DIR__,"../scripts/a_class_kernels.jl"))
using .AClassKernels
@testset "A-class survey contracts" begin
    routes=flavor_routes()
    @test length(routes)==90
    @test length(TRANSITIONS)==9
    @test Set(r.sector for r in routes)==Set(["isovector","isoscalar","strange","charmed","bottom-flavored"])
    @test all(r.emitter in (:u,:d,:s) for r in routes)
    @test all(!all(x in (:c,:b) for x in r.parent) for r in routes)
    @test length(unique((r.sector,r.parent,r.daughter,r.emitted,r.topology) for r in routes))==90
    masses=QuarkMassTable("u"=>0.22,"d"=>0.22,"s"=>0.419,"c"=>1.628,"b"=>4.977)
    # Independent Gaussian integral prediction under the EXISTING Eq19 routing.
    # This explicitly exposes, rather than silently fixing, the emitter/spectator
    # disagreement in the requested unequal-mass extrapolation.
    for t in filter(t->t.family=="A", TRANSITIONS), topology in (:quark,:antiquark), pair in ((:u,:d),(:s,:u),(:u,:c))
        matches=filter(r->r.parent==pair && r.topology==topology && r.parent==r.daughter,routes)
        isempty(matches) && continue
        r=first(matches)
        wp=OscillatorWave(orbital_angular_momentum(t.L),0.4,[1.0])
        ws=OscillatorWave(0,0.4,[1.0])
        c=coefficients(t,r,wp,ws,masses,0.31)
        expected=masses[string(r.spectator)]/(2*(masses[string(r.emitter)]+masses[string(r.spectator)]))
        @test c.h/c.g ≈ expected rtol=1e-10
        @test compare_couplings(c,c).U_hg ≈ 1
    end
    # All three SHO families, both emitter lines, and finite-q/threshold limits.
    for t in TRANSITIONS, topology in (:quark,:antiquark), q in (0.0075,0.15,0.6)
        r=first(filter(r->r.parent==(:u,:d) && r.daughter==r.parent && r.topology==topology,routes))
        wp=OscillatorWave(orbital_angular_momentum(t.L),0.4,[1.0])
        wd=OscillatorWave(orbital_angular_momentum(t.daughter_L),0.4,[1.0])
        c=coefficients(t,r,wp,wd,masses,q)
        if t.family=="A"
            @test c.h/c.g ≈ t.sho_hg rtol=1e-8
        else
            @test isfinite(c.h/c.g)
        end
        @test compare_couplings(c,c).U_hg ≈ 1
    end
    for t in TRANSITIONS[8:9], topology in (:quark,:antiquark)
        r=first(filter(r->r.parent==(:u,:d) && r.daughter==r.parent && r.topology==topology,routes))
        wp=OscillatorWave(orbital_angular_momentum(t.L),0.4,[1.0])
        wd=OscillatorWave(orbital_angular_momentum(t.daughter_L),0.4,[1.0])
        qs=[0.12,0.06,0.03,0.015,0.0075]
        values=[coefficients(t,r,wp,wd,masses,q) for q in qs]
        limit=threshold_limit(qs,values,t.relative_L)
        @test limit.certified
        @test limit.value.h/limit.value.g ≈ t.sho_hg atol=1e-8
        # Independent sampled-wave derivative/overlap representation.
        mesh=collect(range(0.01,30.0;step=0.01))
        sampled=coefficients(t,r,sample_wave(wp,mesh),sample_wave(wd,mesh),masses,0.3)
        analytic=coefficients(t,r,wp,wd,masses,0.3)
        @test sampled.g ≈ analytic.g rtol=3e-5
        @test sampled.h ≈ analytic.h rtol=3e-5
    end
    c=(g=2.0+1im,h=0.4-0.2im)
    phased=(g=im*c.g,h=im*c.h)
    @test compare_couplings(phased,c).vector_distance < 1e-14
    @test compare_couplings(phased,c).U_hg ≈ 1
    tiny=compare_couplings((g=1e-14,h=1.0),(g=2e-14,h=1.0))
    @test !tiny.scalar_stable
    @test tiny.vector_valid
    @test tiny.vector_distance < 1e-12
    @test isnan(real(tiny.U_hg))
    @test !compare_couplings((g=0,h=0),(g=0,h=0)).vector_valid
    q=[0.12,0.06,0.03,0.015]
    for L in 1:3
        values=[(g=x^L*(2+3x^2),h=x^L*(0.5-x^2)) for x in q]
        limit=threshold_limit(q,values,L)
        @test limit.certified
        @test limit.value.g ≈ 2
        @test limit.value.h ≈ 0.5
    end
    @test_throws ArgumentError threshold_limit([0.1,0.0,0.01],[(g=1,h=1) for _ in 1:3],1)
end
