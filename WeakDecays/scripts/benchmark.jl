using GIModel, GIModel.QuarkModelTransitions, Serialization, TOML, LinearAlgebra, Printf
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); method=isempty(ARGS) ? "ho48" : ARGS[1]
data=deserialize(joinpath(root,"results","waves_"*method*".jls"))
channels=TOML.parsefile(joinpath(root,"data","channels.toml"))["channel"]
function model(r,scheme;nr=100,nz=40,pmax=10.)
    CurrentModel(data.waves[r["wi"]],data.waves[r["wf"]],data.masses[r["mi"]],data.masses[r["mf"]],data.masses[r["ms"]];scheme,nr,nz,pmax)
end
open(joinpath(root,"results","rates_"*method*".csv"),"w") do io
    println(io,"channel,scheme,width_GeV,bf,bf_reference,ratio,FL,Rtau_mu,fp0,V0,A00,A10,A20")
    for r in channels, scheme in (:galilean,:boosted_spectator)
        c=model(r,scheme); M,m=r["M"],r["m"]
        result=integrated_rate(c,M,m,r["ml"],r["ckm"];vector=r["vector"])
        ff=form_factors(c,M,m,0.;vector=r["vector"])
        bf=result.width*r["tau"]/HBAR
        ratio=r["bf_ref"]>0 ? bf/r["bf_ref"] : NaN
        tau=integrated_rate(c,M,m,1.77686,r["ckm"];vector=r["vector"]).width
        mu=integrated_rate(c,M,m,.1056583755,r["ckm"];vector=r["vector"]).width
        nums=r["vector"] ? [NaN,ff.V,ff.A0,ff.A1,ff.A2] : [ff.fp,NaN,NaN,NaN,NaN]
        println(io,join([r["name"],scheme,result.width,bf,r["bf_ref"],ratio,result.FL,tau/mu,nums...],","))
        @info "Rate" method channel=r["name"] scheme bf ratio
    end
end
open(joinpath(root,"results","leptonic_"*method*".csv"),"w") do io
    println(io,"meson,mass_GI,mass_external,f_MeV,reference_MeV,ratio,width_mu_GeV,bf_mu")
    for r in TOML.parsefile(joinpath(root,"data","leptonic.toml"))["meson"]
        original=data.states[r["key"]]
        state=PhysicalState(r["name"],r["M"],original.components)
        f=leptonic_constant(state,data.masses)
        ml=.1056583755
        width=r["ckm"]>0 ? GF^2*r["ckm"]^2*f^2*ml^2*r["M"]/(8pi)*(1-ml^2/r["M"]^2)^2 : NaN
        println(io,join([r["name"],original.mass_GeV,r["M"],1000f,r["fref_MeV"],r["fref_MeV"]>0 ? 1000f/r["fref_MeV"] : NaN,width,width*r["tau"]/HBAR],","))
    end
end
open(joinpath(root,"results","curves_"*method*".csv"),"w") do io
    println(io,"channel,scheme,q2,fp,f0,V,A0,A1,A2,rate,FL")
    for r in channels, scheme in (:galilean,:boosted_spectator)
        c=model(r,scheme); M,m=r["M"],r["m"]
        for q2 in range(0.,.999*(M-m)^2,length=61)
            ff=form_factors(c,M,m,q2;vector=r["vector"])
            rr=differential_rate(c,M,m,q2,0.,r["ckm"];vector=r["vector"])
            nums=r["vector"] ? [NaN,NaN,ff.V,ff.A0,ff.A1,ff.A2] : [ff.fp,ff.f0,NaN,NaN,NaN,NaN]
            println(io,join([r["name"],scheme,q2,nums...,rr.total,rr.total>0 ? rr.longitudinal/rr.total : 0.],","))
        end
    end
end
if method=="ho48"
    open(joinpath(root,"results","convergence.csv"),"w") do io
        println(io,"channel,scheme,nr,nz,nq,pmax,width_GeV")
        for name in ("D0_K","D0_Kstar","B0_pi","B0_Dstar"), scheme in (:galilean,:boosted_spectator)
            r=only(filter(r->r["name"]==name,channels))
            for (nr,nz,nq,pmax) in [(64,24,32,8.),(100,40,48,10.),(160,64,80,12.)]
                c=model(r,scheme;nr,nz,pmax)
                rate=integrated_rate(c,r["M"],r["m"],r["ml"],r["ckm"];vector=r["vector"],nq)
                println(io,join([name,scheme,nr,nz,nq,pmax,rate.width],","))
            end
        end
    end
    open(joinpath(root,"results","resonances.csv"),"w") do io
        println(io,"channel,scheme,running,R,width_scale,folded_GeV,stable_GeV,ratio")
        for (name,gamma,ma,mb) in [("D0_Kstar",.0485,.493677,.13957039),("B0_rho",.1491,.13957039,.13957039)]
            r=only(filter(r->r["name"]==name,channels))
            for scheme in (:galilean,:boosted_spectator)
                c=model(r,scheme;nr=80,nz=32)
                stable=integrated_rate(c,r["M"],r["m"],0.,r["ckm"];vector=true,nq=32).width
                for (running,R,scale) in [(false,3.,.0001),(false,3.,1.),(true,1.,1.),(true,3.,1.),(true,5.,1.)]
                    folded=folded_rate(c,r["M"],r["m"],gamma*scale,ma,mb,0.,r["ckm"];ns=48,nq=32,R,running)
                    println(io,join([name,scheme,running,R,scale,folded,stable,folded/stable],","))
                end
            end
        end
    end
    open(joinpath(root,"results","nonleptonic.csv"),"w") do io
        println(io,"channel,scheme,a1,fpion_source,bf")
        fgi=leptonic_constant(PhysicalState("pi",.13957039,data.states["ud_1_0"].components),data.masses)
        for name in ("B0_D","D0_K"), scheme in (:galilean,:boosted_spectator), a1 in (.9,1.,1.1), (label,f) in (("external",.1302),("GI",fgi))
            r=only(filter(r->r["name"]==name,channels)); c=model(r,scheme)
            width=factorized_pp_width(c,r["M"],r["m"],.13957039,f,r["ckm"]*.97435;a1)
            println(io,join([name*"_pi",scheme,a1,label,width*r["tau"]/HBAR],","))
        end
    end
    # Radial excitations use MODEL external daughter masses: sensitivity studies,
    # no identification with a specific observed broad resonance is asserted.
    open(joinpath(root,"results","radial_excitations.csv"),"w") do io
        println(io,"channel,scheme,mass_GI,bf,overlap_zero_recoil")
        for name in ("B0_D","B0_Dstar","B0_pi","B0_rho","D0_K","D0_Kstar"), scheme in (:galilean,:boosted_spectator)
            r=copy(only(filter(r->r["name"]==name,channels)))
            r["wf"]=replace(r["wf"],"_1_"=>"_2_"); m=data.states[r["wf"]].mass_GeV
            c=model(r,scheme); rate=integrated_rate(c,r["M"],m,0.,r["ckm"];vector=r["vector"])
            overlap=sum(gauss(160,0.,12.) |> x-> [w*p^2*c.wi(p)*c.wf(p) for (p,w) in zip(x...)])
            println(io,join([name*"_2S",scheme,m,rate.width*r["tau"]/HBAR,overlap],","))
        end
    end
end
