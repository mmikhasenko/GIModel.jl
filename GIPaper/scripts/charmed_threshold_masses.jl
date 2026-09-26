#!/usr/bin/env julia
# julia --project=GIPaper/scripts GIPaper/scripts/charmed_threshold_masses.jl
using GIModel, TOML, SHA, LinearAlgebra
BLAS.set_num_threads(1)
root=normpath(joinpath(@__DIR__,"../.."))
params,mq=load_parameters_and_quark_masses(default_parameters_path())
levels=BasisState[]
for (L,nmax) in (("S",3),("P",2),("D",2),("F",1),("G",1))
    append!(levels,spectrum_levels(nmax;L_labels=(L,)))
end
refined="--refined" in ARGS
solver=FiniteDifferenceSolver(ngrid=refined ? 900 : 450,rmax=24.0,kinetic=:relativistic)
rows=Dict[]
for (sector,f) in (("D",:q),("Ds",:s))
    println("Computing ",sector);flush(stdout)
    spectrum=compute_spectrum(params,Meson(mq,:c,f);levels,solver)
    for s in spectrum.states
        ((s.L=="S" && s.n<=2) || (s.n==1 && s.L in ("P","D"))) || continue
        push!(rows,Dict("sector"=>sector,"n"=>s.n,"L"=>String(s.L),
            "multiplicity"=>s.multiplicity,"J"=>s.J,"label"=>String(s.label),
            "unmixed_GeV"=>s.corrected.mass_GeV,"mixed_GeV"=>s.mass_GeV,
            "mixing_count"=>length(s.mixings)))
    end
end
out=joinpath(root,"GIPaper/docs/pseudoscalar_census/charmed_threshold_masses"*(refined ? "_refined" : "")*".toml")
open(out,"w") do io
    TOML.print(io,Dict("solver"=>repr(solver),"parameter_sha256"=>bytes2hex(sha256(read(default_parameters_path()))),
        "prescription"=>"fixed-sector corrected masses and final mixed masses from compute_spectrum",
        "level_inventory"=>"3S,2P,2D,1F,1G (same as full spectrum figure)",
        "julia_version"=>string(VERSION),"states"=>rows))
end
println(out)
