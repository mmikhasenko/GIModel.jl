using GIModel, GIModel.QuarkModelTransitions, Serialization, Printf
include("../src/WeakDecayExploration.jl")
using .WeakDecayExploration
LinearAlgebra = GIModel.LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,".."))
params,masses=load_parameters_and_quark_masses(default_parameters_path())
method=length(ARGS)>0 ? ARGS[1] : "ho32"
solver = startswith(method,"fd") ? FiniteDifferenceSolver(ngrid=parse(Int,method[3:end]),rmax=24.,nlevels_per_channel=2) :
    OscillatorSolver(nbasis=parse(Int,method[3:end]),converge=false,nlevels_per_channel=2)
levels=[BasisState(n,"S",mult,J) for n in 1:2 for (mult,J) in ((1,0),(3,1))]
states=Dict(); waves=Dict(); rows=[]
for (sector,a,b) in [("ud",:u,:d),("su",:s,:u),("cu",:c,:u),("cs",:c,:s),("bu",:b,:u),("bs",:b,:s),("bc",:b,:c),("ss",:s,:s)]
    @info "Solving" method sector
    spec=compute_spectrum(params,Meson(masses,a,b);levels,solver)
    for l in levels
        state=physical_state(spec,l); comp=only(state.components)
        key=sector*"_"*string(l.n)*"_"*string(l.J)
        states[key]=state; waves[key]=prepare_wave(comp.wave)
        push!(rows,(key,state.mass_GeV,wave_mean_squares(comp.wave,0).r2))
    end
end
serialize(joinpath(root,"results","waves_"*method*".jls"),(;states,waves,masses,method))
open(joinpath(root,"results","spectrum_"*method*".csv"),"w") do io
    println(io,"state,mass_GeV,r2_GeVminus2")
    for r in rows; println(io,join(r,",")); end
end
