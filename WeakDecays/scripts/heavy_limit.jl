using GIModel,GIModel.QuarkModelTransitions,Serialization,LinearAlgebra
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); data=deserialize(joinpath(root,"results","waves_ho48.jls"))
# Fixed common light-cloud wave: isolates current symmetry from finite-mass
# wavefunction differences. This is NOT a newly solved static-limit spectrum.
wave=data.waves["bu_1_0"]
function hfunctions(c,M,m,w)
    q2=M*M+m*m-2M*m*w
    h=current_components(c,M,m,q2)
    scale=sqrt(M*m); rootw=sqrt(w*w-1)
    hp=(h.V0+(1-w)*h.Vz/rootw)/(2scale)
    hm=hp-h.Vz/(scale*rootw)
    hA1=h.Ax/(scale*(w+1)); hV=h.Vy/(scale*rootw)
    (;hp,hm,hA1,hV)
end
open(joinpath(root,"results","heavy_limit.csv"),"w") do io
    println(io,"scheme,heavy_scale,w,hplus,hminus,hA1,hV,slope")
    for scheme in (:galilean,:boosted_spectator), scale in (1.,4.,16.,64.)
        mi=4.977scale; mf=1.628scale; M=mi+.4; m=mf+.4
        c=CurrentModel(wave,wave,mi,mf,.22;scheme,nr=160,nz=64,pmax=12.)
        a=hfunctions(c,M,m,1.001); b=hfunctions(c,M,m,1.002)
        slope=-(b.hp-a.hp)/.001
        for w in (1.001,1.1,1.2,1.4)
            h=hfunctions(c,M,m,w)
            println(io,join([scheme,scale,w,h.hp,h.hm,h.hA1,h.hV,slope],","))
        end
    end
end
