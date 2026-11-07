using GIModel,GIModel.QuarkModelTransitions,Serialization,LinearAlgebra,Test
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); d=deserialize(joinpath(root,"results","waves_ho48.jls"))
M=1.86484; ma=.493677; mb=.13957039; poles=[.89188,1.5791029147667093]; widths=[.0485,.2]; q2=.02
rho(s)=s>(ma+mb)^2 ? 2recoil(sqrt(s),ma,mb^2)/sqrt(s) : 0.
g=sqrt.([poles[i]*widths[i]/rho(poles[i]^2) for i in 1:2])
# Production strengths are diagnostic transverse-current overlaps; these are
# not normalized 1->2 weak amplitudes, and strong relative signs are unknown.
beta=[]
for (key,m) in zip(["su_1_1","su_2_1"],poles)
    c=CurrentModel(d.waves["cu_1_0"],d.waves[key],d.masses["c"],d.masses["s"],d.masses["u"])
    push!(beta,current_components(c,M,m,q2).Ax)
end
function amplitudes(s,sign)
    D=poles.^2 .-s; r=rho(s)
    den=D[1]*D[2]-im*r*(g[1]^2*D[2]+g[2]^2*D[1])
    t=(g[1]^2*D[2]+g[2]^2*D[1])/den
    F=(beta[1]*g[1]*D[2]+sign*beta[2]*g[2]*D[1])/den
    terms=[beta[1]*g[1]/(D[1]-im*r*g[1]^2),sign*beta[2]*g[2]/(D[2]-im*r*g[2]^2)]
    (;t,F,bw=sum(terms),incoherent=sum(abs2,terms))
end
@testset "Elastic production unitarity (two-pole toy)" begin
    for s in range((ma+mb+.01)^2,(M-sqrt(q2))^2,length=19), sign in (-1,1)
        a=amplitudes(s,sign)
        @test imag(a.t) ≈ rho(s)*abs2(a.t) atol=1e-12
        @test imag(a.F) ≈ conj(a.t)*rho(s)*a.F atol=1e-12
    end
end
open(joinpath(root,"results","rescattering.csv"),"w") do io
    println(io,"mass,sign,incoherent,coherent_bw,kmatrix,production_beta1,production_beta2")
    for m in range(ma+mb+.001,M-sqrt(q2)-.001,length=500), sign in (-1,1)
        a=amplitudes(m*m,sign); phase=rho(m*m)*recoil(M,m,q2)
        println(io,join([m,sign,phase*a.incoherent,phase*abs2(a.bw),phase*abs2(a.F),beta...],","))
    end
end
