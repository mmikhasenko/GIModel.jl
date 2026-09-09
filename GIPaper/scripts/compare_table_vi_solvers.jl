#!/usr/bin/env julia
# Compare already-generated native-HO and independent-FD Table VI results.
using Pkg
Pkg.activate(dirname(@__DIR__))
using CSV, Printf, Statistics

root = dirname(@__DIR__)
folder = joinpath(root, "docs", "residual_reports")
ho = collect(CSV.File(joinpath(folder, "table_vi_photon_decays.csv")))
fd = collect(CSV.File(joinpath(folder, "table_vi_photon_decays_fd.csv")))
key(r) = (String(r.multipole), String(r.parent), String(r.daughter))
@assert length(ho) == length(fd) == 79
fmap = Dict(key(r) => r for r in fd)
@assert Set(key.(ho)) == Set(keys(fmap))
path = joinpath(folder, "table_vi_solver_comparison.md")
open(path, "w") do io
    println(io, "# Table VI independent solver comparison")
    print(io, """
    
    Compares the complete native-HO run with the independent FD run at
    ngrid=700, rmax=24 GeV^-1. This is an observable cross-check at the stated
    resolution, not a new FD precision certification. Both runs use identical
    canonical states, charges, kinematic inputs, and shared phase conventions.
    
    | Multipole | Rows | Median relative difference | Largest absolute difference |
    |---|---:|---:|---:|
    """)
    for kind in ("M1", "E1", "M2")
        block = filter(r -> r.multipole == kind, ho)
        differences = [abs(r.computed - fmap[key(r)].computed) for r in block]
        relative = [d / max(abs(r.computed), abs(fmap[key(r)].computed))
                    for (r, d) in zip(block, differences)]
        @printf(io, "| %s | %d | %.3f%% | %.6g |\n",
            kind, length(block), 100median(relative), maximum(differences))
    end
    print(io, """
    
    Absolute differences are in μN for M1 and MeV^(1/2) for E1/M2.
    The following rows expose radial mixing and cancellation sensitivity.
    
    | Parent -> daughter | HO | FD | Paper |
    |---|---:|---:|---:|
    """)
    for r in ho
        if r.parent in ("eta_r", "etaprime_r") ||
           (r.parent == "Upsilondoubleprime" && r.daughter in ("eta_b", "etaprime_b", "chi_0b")) ||
           (r.parent == "Upsilon" && r.daughter == "etaprime")
            @printf(io, "| %s -> %s | %+.6g | %+.6g | %+.6g |\n",
                r.parent, r.daughter, r.computed, fmap[key(r)].computed, r.paper)
        end
    end
    print(io, """
    
    The excited eta magnitude residuals survive the change of numerical
    representation. They are not resolved by replacing the old central P
    waves or fixing the duplicated perfect-mixing coefficient. Near-zero
    amplitudes are shown in absolute units so a large relative error does not
    conceal the cancellation scale. Canonical paper values are never replaced
    by either solver's answer.
    """)
end
println("wrote ", path)
