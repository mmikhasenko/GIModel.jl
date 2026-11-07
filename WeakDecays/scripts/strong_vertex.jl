using GIModel,GIModel.QuarkModelTransitions,Serialization,LinearAlgebra
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); d=deserialize(joinpath(root,"results","waves_ho48.jls"))
function state(key,label,m,flavors)
    c=only(d.states[key].components); b=c.basis
    PhysicalState(label,m,[(basis=BasisState(b.n,b.L_label,b.multiplicity,b.J;flavors),coefficient=1.,wave=c.wave)])
end
piw=only(d.states["ud_1_0"].components).wave
pi0=PhysicalState("pi0",.1349768,[(basis=BasisState(1,"S",1,0;flavors=(f,f)),coefficient=sgn/sqrt(2),wave=piw) for (f,sgn) in ((:u,1),(:d,-1))])
piminus=state("ud_1_0","pi-",.13957039,(:d,:u))
km=state("su_1_0","K-",.493677,(:s,:u)); k0=state("su_1_0","Kbar0",.497611,(:s,:d))
function total(m,h_over_g)
    initial=state("su_1_1","K*-",m,(:s,:u))
    op=PseudoscalarEmission(1.,h_over_g,d.masses)
    a=matrix_element(TwoMesonChannel(km,pi0),op,initial)
    b=matrix_element(TwoMesonChannel(k0,piminus),op,initial)
    (decay_width(a)+decay_width(b))/1000 # QMT MeV -> GeV
end
open(joinpath(root,"results","strong_vertex.csv"),"w") do io
    println(io,"mass,h_over_g,gfit_GeVminus1,width_raw_GeV,width_calibrated_GeV")
    for h in (0.,.5,1.)
        reference=total(.89188,h); gfit=sqrt(.0485/reference)
        for m in range(.65,1.4,length=31)
            width=total(m,h)
            println(io,join([m,h,gfit,width,width*gfit^2],","))
        end
    end
end
