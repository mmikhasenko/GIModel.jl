using GIModel,GIModel.QuarkModelTransitions,Serialization,TOML,LinearAlgebra
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,".."));d=deserialize(joinpath(root,"results","waves_ho48.jls"))
channels=TOML.parsefile(joinpath(root,"data","channels.toml"))["channel"]
model(r,scheme)=CurrentModel(d.waves[r["wi"]],d.waves[r["wf"]],d.masses[r["mi"]],d.masses[r["mf"]],d.masses[r["ms"]];scheme,nr=100,nz=40)
open(joinpath(root,"results","extended_convergence.csv"),"w") do io
    println(io,"experiment,channel,scheme,coarse,fine,relative_difference")
    for (name,gamma,ma,mb) in [("D0_Kstar",.0485,.493677,.13957039),("B0_rho",.1491,.13957039,.13957039)],scheme in (:galilean,:boosted_spectator)
        r=only(filter(r->r["name"]==name,channels));c=model(r,scheme)
        a=folded_rate(c,r["M"],r["m"],gamma,ma,mb,0.,r["ckm"];ns=48,nq=32)
        b=folded_rate(c,r["M"],r["m"],gamma,ma,mb,0.,r["ckm"];ns=96,nq=64)
        delta=abs(a/b-1);@assert delta<.001
        println(io,join(["spectral",name,scheme,a,b,delta],","))
    end
    for name in ("D0_pi","B0_pi")
        r=only(filter(r->r["name"]==name,channels));c=model(r,:galilean)
        a=lightfront_fp(c,4.;nx=100,nk=100,nphi=32)
        b=lightfront_fp(c,4.;nx=160,nk=160,nphi=64)
        delta=abs(a/b-1);@assert delta<.001
        println(io,join(["lightfront_Q2_4",name,"plus_current",a,b,delta],","))
    end
end
