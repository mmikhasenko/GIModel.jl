#!/usr/bin/env julia
# julia --project=GIPaper/scripts GIPaper/scripts/five_sector_threshold_masses.jl [--refined]
using GIModel, TOML, SHA, LinearAlgebra
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const OUT=joinpath(ROOT,"GIPaper/docs/pseudoscalar_census")
params,mq=load_parameters_and_quark_masses(default_parameters_path())
levels=BasisState[]
for (L,nmax) in (("S",3),("P",2),("D",2),("F",1),("G",1))
    append!(levels,spectrum_levels(nmax;L_labels=(L,)))
end
refined="--refined" in ARGS
suffix=refined ? "_refined" : ""
solver=FiniteDifferenceSolver(ngrid=refined ? 900 : 450,rmax=24.0,kinetic=:relativistic)
# Reuse exactly the independently recorded charm calculation.
charm=TOML.parsefile(joinpath(OUT,"charmed_threshold_masses"*suffix*".toml"))
@assert charm["solver"]==repr(solver)
@assert charm["parameter_sha256"]==bytes2hex(sha256(read(default_parameters_path())))
rows=copy(charm["states"])
for (sector,f1,f2) in (("nn",:q,:q),("ss",:s,:s),("K",:q,:s),("B",:b,:q),("Bs",:b,:s))
    println("Computing ",sector," ngrid=",refined ? 900 : 450);flush(stdout)
    spectrum=compute_spectrum(params,Meson(mq,f1,f2);levels,solver)
    for s in spectrum.states
        ((s.L=="S" && s.n<=2) || (s.n==1 && s.L in ("P","D"))) || continue
        push!(rows,Dict("sector"=>sector,"n"=>s.n,"L"=>String(s.L),"multiplicity"=>s.multiplicity,
            "J"=>s.J,"label"=>String(s.label),"unmixed_GeV"=>s.corrected.mass_GeV,
            "mixed_GeV"=>s.mass_GeV,"mixing_count"=>length(s.mixings)))
    end
end
open(joinpath(OUT,"five_sector_masses"*suffix*".toml"),"w") do io
    TOML.print(io,Dict("solver"=>repr(solver),"parameter_sha256"=>charm["parameter_sha256"],
        "prescription"=>"fixed-sector basis masses, no flavor annihilation; charm reused from prior audit",
        "level_inventory"=>charm["level_inventory"],"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "states"=>rows))
end
