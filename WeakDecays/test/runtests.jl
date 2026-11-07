using Test, GIModel, LinearAlgebra
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
@testset "Exploratory weak current: independent checks" begin
    x,w=gauss(12,-1.,1.)
    @test sum(w.*x.^8) ≈ 2/9 atol=1e-14
    wave=prepare_wave(OscillatorWave(0,0.5,[1.]);n=12001)
    c=CurrentModel(wave,wave,1.5,1.5,.3;nr=180,nz=48)
    # Flavor diagonal charge normalization at zero spatial transfer (M != m
    # only because decay API requires distinct external masses).
    h=current_components(c,2.1,2.0,.1^2-1e-15)
    @test h.V0/(2sqrt(2.1*2.0)) ≈ 1 rtol=3e-6
    @test abs(h.Vz)<1e-6
    f=form_factors(c,3.,2.,0.)
    @test f.fp==f.f0
    @test differential_rate(c,3.,2.,0.,0.,.9).total ≈ GF^2*.9^2*recoil(3.,2.,0.)^3*f.fp^2/(24pi^3) rtol=1e-12
    for q2 in (.1,.5,.9), ml in (0.,.105)
        ff=form_factors(c,3.,2.,q2); p=recoil(3.,2.,q2)
        # Independently expressed invariant-form-factor rate.
        lam=(2*3*p)^2
        expected=GF^2*.9^2*p/(96pi^3*3^2)*(1-ml^2/q2)^2 *
            ((1+ml^2/(2q2))*lam*ff.fp^2+3ml^2/(2q2)*(3^2-2^2)^2*ff.f0^2)
        @test differential_rate(c,3.,2.,q2,ml,.9).total ≈ expected rtol=1e-12
        vf=form_factors(c,3.,2.,q2;vector=true)
        hp=(3+2)*vf.A1-sqrt(lam)/(3+2)*vf.V
        hm=(3+2)*vf.A1+sqrt(lam)/(3+2)*vf.V
        h0=((3^2-2^2-q2)*(3+2)*vf.A1-lam/(3+2)*vf.A2)/(2*2*sqrt(q2))
        ht=sqrt(lam/q2)*vf.A0
        vexpected=GF^2*.9^2*p*q2/(96pi^3*3^2)*(1-ml^2/q2)^2 *
            ((1+ml^2/(2q2))*(hp^2+hm^2+h0^2)+3ml^2/(2q2)*ht^2)
        @test differential_rate(c,3.,2.,q2,ml,.9;vector=true).total ≈ vexpected rtol=1e-12
    end
    v=form_factors(c,3.,2.,0.;vector=true)
    @test v.A0 ≈ ((3+2)*v.A1-(3-2)*v.A2)/(2*2) atol=1e-12
    @test integrated_rate(c,3.,2.,1.1,.9).width==0
    @test differential_rate(c,3.,2.,.5,0.,0.).total==0
    @test_throws ArgumentError form_factors(c,3.,2.,1.)
    @test_throws ArgumentError CurrentModel(wave,wave,1.,1.,.3;scheme=:unknown)
    # Radial orthogonality tests node signs and Fourier transform phase.
    excited=prepare_wave(OscillatorWave(0,.5,[0.,1.]);n=12001)
    node=CurrentModel(wave,excited,1.5,1.5,.3;nr=180,nz=48)
    @test abs(current_components(node,2.1,2.,.1^2-1e-15).V0)<2e-5
    # Independent spinor-matrix calculation at generic three-dimensional momenta.
    I2=Matrix{ComplexF64}(I,2,2); Z=zeros(ComplexF64,2,2)
    sig=[[0 1;1 0], [0 -im;im 0], [1 0;0 -1]]
    g0=[I2 Z;Z -I2]; gi=[[Z s;-s Z] for s in sig]
    g5=im*g0*gi[1]*gi[2]*gi[3]
    pi3=[.2,.3,-.1]; pf3=pi3+[0,0,.7]; mi=1.5; mf=.4
    ei=sqrt(mi^2+sum(abs2,pi3)); ef=sqrt(mf^2+sum(abs2,pf3))
    ui=sqrt((ei+mi)/(2ei))*[I2;sum(sig[j]*pi3[j] for j in 1:3)/(ei+mi)]
    uf=sqrt((ef+mf)/(2ef))*[I2;sum(sig[j]*pf3[j] for j in 1:3)/(ef+mf)]
    norm=sqrt((ei+mi)*(ef+mf)/(4ei*ef)); den=(ei+mi)*(ef+mf)
    @test real(tr(uf'*ui)/2) ≈ norm*(1+dot(pi3,pf3)/den)
    @test real(tr(sig[3]*uf'*g0*gi[3]*g5*ui)/2) ≈ norm*(1+(2pi3[3]*pf3[3]-dot(pi3,pf3))/den)
    @test imag(tr(sig[1]*uf'*g0*gi[2]*ui)/2) ≈ norm*(pi3[3]/(ei+mi)-pf3[3]/(ef+mf))
    # Folded width converges to stable width, without arbitrary renormalization.
    stable=integrated_rate(c,3.,2.,0.,.9;vector=true,nq=24).width
    folded=folded_rate(c,3.,2.,1e-6,.5,.2,0.,.9;ns=40,nq=24,running=false)
    @test folded/stable ≈ 1 rtol=1e-4
    @test spectral_density(.1,1.,.1,.5,.2)==0
end
@testset "Coherent S-wave state adapter" begin
    using GIModel.QuarkModelTransitions
    w=OscillatorWave(0,.5,[1.])
    initial=PhysicalState("D",1.87,[(basis=BasisState(1,"S",1,0;flavors=(:c,:u)),coefficient=1.,wave=w)])
    final=PhysicalState("K",.494,[(basis=BasisState(1,"S",1,0;flavors=(:s,:u)),coefficient=1.,wave=w)])
    mq=QuarkMassTable("c"=>1.628,"s"=>.419,"u"=>.22)
    c=current_from_states(final,initial,mq;active=:c=>:s)
    width=integrated_rate(c,1.87,.494,0.,.97).width
    phase=CoherentCurrent([(im,c.terms[1][2])])
    @test integrated_rate(phase,1.87,.494,0.,.97).width ≈ width
    cancel=CoherentCurrent([(1.0 + 0im,c.terms[1][2]),(-1.0 + 0im,c.terms[1][2])])
    @test integrated_rate(cancel,1.87,.494,0.,.97).width == 0
    @test_throws ArgumentError current_from_states(final,initial,mq;active=:b=>:s)
    p=PhysicalState("axial",1.2,[(basis=BasisState(1,"P",3,1;flavors=(:s,:u)),coefficient=1.,wave=OscillatorWave(1,.5,[1.]))])
    @test_throws ArgumentError current_from_states(p,initial,mq;active=:c=>:s)
end
