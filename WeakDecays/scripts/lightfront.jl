using GIModel,GIModel.QuarkModelTransitions,Serialization,TOML,LinearAlgebra,Test
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"..")); data=deserialize(joinpath(root,"results","waves_ho48.jls"))
channels=TOML.parsefile(joinpath(root,"data","channels.toml"))["channel"]
@testset "Light-front normalization and symmetry" begin
    w=prepare_wave(OscillatorWave(0,.5,[1.]);n=12001)
    c=CurrentModel(w,w,1.5,1.5,.22)
    @test lightfront_fp(c;nx=160,nk=160) ≈ 1 rtol=1e-5
    c1=CurrentModel(w,w,1.5,.4,.22); c2=CurrentModel(w,w,.4,1.5,.22)
    @test lightfront_fp(c1) ≈ lightfront_fp(c2) rtol=1e-14
    @test 0<lightfront_fp(c1)<1
    @test_throws ArgumentError lightfront_fp(c,-1.)
    x=.7; kt2=.12; a=lf_kinematics(x,kt2,1.5,.22); epsx=1e-6
    kz(x)=((x-.5)*sqrt((1.5^2+kt2)/x+(.22^2+kt2)/(1-x))+(.22^2-1.5^2)/(2sqrt((1.5^2+kt2)/x+(.22^2+kt2)/(1-x))))
    @test (kz(x+epsx)-kz(x-epsx))/(2epsx) ≈ a.jac rtol=1e-8
end
open(joinpath(root,"results","lightfront.csv"),"w") do io
    println(io,"channel,Q2,n,spin,fp")
    for r in filter(r->!r["vector"],channels)
        c=CurrentModel(data.waves[r["wi"]],data.waves[r["wf"]],data.masses[r["mi"]],data.masses[r["mf"]],data.masses[r["ms"]])
        for (Q2,n,spin) in [(0.,64,true),(0.,100,true),(0.,160,true),(0.,100,false),(1.,100,true),(2.,100,true),(4.,100,true)]
            fp=lightfront_fp(c,Q2;nx=n,nk=n,spin)
            println(io,join([r["name"],Q2,n,spin,fp],","))
        end
    end
end
