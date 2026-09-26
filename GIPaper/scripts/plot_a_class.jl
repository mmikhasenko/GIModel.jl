#!/usr/bin/env julia
# julia --project=GIPaper/scripts GIPaper/scripts/plot_a_class.jl
using CSV, CairoMakie, Printf, TOML, SHA
const OUT=normpath(joinpath(@__DIR__,"../docs/a_class_universality"))
const SECTORS=["isovector","isoscalar","strange","charmed","bottom-flavored"]
const COLORS=["#0072B2","#E69F00","#009E73","#CC79A7","#D55E00"]
const LABELS=["³S₁ → ¹S₀ · P wave","³P₂ → ¹S₀ · D wave","³P₂ → ³S₁ · D wave",
              "³P₁ → ³S₁ · D wave","¹P₁ → ³S₁ · D wave","³D₃ → ¹S₀ · F wave","³D₃ → ³S₁ · F wave"]
readtable(name) = begin
    path=joinpath(OUT,name*".csv")
    CSV.File(isfile(path) ? path : read(`gzip -dc $(path*".gz")`))
end
rows=collect(readtable("comparisons"))
primary=filter(r->r.backend=="HO",rows)
set_theme!(Theme(font="TeX Gyre Heros",fontsize=17,Axis=(xgridvisible=false,ygridcolor=(:gray,0.15),)))
valid=filter(r->r.scalar_stable,primary)
lo=min(1.0,minimum(r.U_hg_re for r in valid)); hi=max(1.0,maximum(r.U_hg_re for r in valid))
pad=max(0.03,0.03*(hi-lo)); limits=(lo-pad,hi+pad)
bins=collect(range(limits...;length=35))
membership=NamedTuple[]
fig=Figure(size=(1250,870))
for (i,q) in enumerate([0.0,0.15,0.30,0.60])
    selected=filter(r->r.q_GeV==q,primary)
    title=q==0 ? "Threshold: qᴸ-reduced extrapolation" : "Common q = $q GeV"
    ax=Axis(fig[(i-1)÷2+1,(i-1)%2+1],title=title,xlabel="Uₕɡ = (Cₕ/Cɡ)native / (Cₕ/Cɡ)SHO",ylabel="Basis-route count")
    counts=zeros(Int,length(bins)-1,length(SECTORS))
    for (s,sector) in enumerate(SECTORS), r in filter(r->r.sector==sector,selected)
        bin=r.scalar_stable ? clamp(searchsortedlast(bins,r.U_hg_re),1,length(bins)-1) : 0
        bin>0 && (counts[bin,s]+=1)
        push!(membership,(figure="histograms",panel=i,point_id=String(r.point_id),
            bin=bin,x=r.scalar_stable ? r.U_hg_re : r.vector_distance,y=q,
            scalar_stable=r.scalar_stable,certified=r.certified,
            certificate_id=String(r.certificate_id),native_g_trace_ids=String(r.native_g_trace_ids),
            native_h_trace_ids=String(r.native_h_trace_ids),sho_g_trace_ids=String(r.sho_g_trace_ids),
            sho_h_trace_ids=String(r.sho_h_trace_ids)))
    end
    centers=(bins[1:end-1]+bins[2:end])/2
    barplot!(ax,repeat(centers,length(SECTORS)),vec(counts);stack=repeat(1:5,inner=length(centers)),
             color=repeat(COLORS,inner=length(centers)),width=0.95*(bins[2]-bins[1]))
    vlines!(ax,[1.0];color=:black,linestyle=:dash,linewidth=2)
    xlims!(ax,limits...)
    failed=count(r->!r.certified,selected); unstable=count(r->!r.scalar_stable,selected)
    text!(ax,0.98,0.97;text="$(length(selected)) routes · $failed uncertified · $unstable unstable",
          space=:relative,align=(:right,:top),fontsize=12)
end
Legend(fig[3,1:2],[PolyElement(color=c) for c in COLORS],SECTORS;orientation=:horizontal,framevisible=false)
Label(fig[0,1:2],"A-class kernels · common β = 0.40 GeV SHO comparator",fontsize=25)
Label(fig[4,1:2],"Formal unmixed basis survey; counts are correlated flavor routes, not independent measurements.\nExisting Eq. (19) routing; unequal-mass emitter-formula mismatch remains unresolved.",fontsize=14)
save(joinpath(OUT,"histograms.pdf"),fig); save(joinpath(OUT,"histograms.png"),fig;px_per_unit=2)

panels=Figure(size=(1550,1250))
for (i,label) in enumerate(LABELS)
    ax=Axis(panels[(i-1)÷3+1,(i-1)%3+1],title="A$i: $label",xlabel="q [GeV] (0 = extrapolation)",ylabel="Uₕɡ")
    selected=filter(r->r.transition=="A$i",primary)
    hlines!(ax,[1.0];color=:black,linestyle=:dash,linewidth=1.5)
    for (s,sector) in enumerate(SECTORS)
        sectorrows=filter(r->r.sector==sector,selected)
        for route in unique(r.route_id for r in sectorrows)
            route_rows=sort(filter(r->r.route_id==route,sectorrows);by=r->r.q_GeV)
            stable=filter(r->r.scalar_stable,route_rows)
            lines!(ax,[r.q_GeV for r in stable],[r.U_hg_re for r in stable];color=(COLORS[s],0.3),linewidth=1.0)
            for r in stable
                scatter!(ax,[r.q_GeV],[r.U_hg_re];color=COLORS[s],markersize=7,
                    marker=r.certified ? :circle : :xcross)
            end
        end
    end
    ylims!(ax,limits...)
    for r in selected
        push!(membership,(figure="transition_panels",panel=i,point_id=String(r.point_id),bin=0,
            x=r.q_GeV,y=r.scalar_stable ? r.U_hg_re : r.vector_distance,
            scalar_stable=r.scalar_stable,certified=r.certified,certificate_id=String(r.certificate_id),
            native_g_trace_ids=String(r.native_g_trace_ids),native_h_trace_ids=String(r.native_h_trace_ids),
            sho_g_trace_ids=String(r.sho_g_trace_ids),sho_h_trace_ids=String(r.sho_h_trace_ids)))
    end
end
Legend(panels[3,2:3],[PolyElement(color=c) for c in COLORS],SECTORS;framevisible=false)
Label(panels[0,1:3],"Seven A-class transitions · all resolved flavor routes",fontsize=26)
Label(panels[4,1:3],"Lines connect the same flavor/topology kernel. Crosses fail numerical certification.\nCoincident charge/isospin routes remain separate entries in plot_membership.csv.",fontsize=16)
save(joinpath(OUT,"transition_panels.pdf"),panels); save(joinpath(OUT,"transition_panels.png"),panels;px_per_unit=2)

# Unstable ratios are never silently dropped: their normalized-vector distance
# has its own figure and membership. Empty in a well-conditioned run.
unstable=filter(r->!r.scalar_stable,primary)
if !isempty(unstable)
    f=Figure(size=(900,500)); ax=Axis(f[1,1],xlabel="Route index",ylabel="Phase-aligned vector distance")
    for (i,r) in enumerate(unstable)
        r.vector_valid && scatter!(ax,[i],[r.vector_distance];color=COLORS[findfirst(==(r.sector),SECTORS)])
    end
    save(joinpath(OUT,"unstable_vectors.pdf"),f)
end
CSV.write(joinpath(OUT,"plot_membership.csv"),membership)
cert=collect(readtable("certificates"))
open(joinpath(OUT,"numerical_summary.md"),"w") do io
    println(io,"# Generated numerical findings\n")
    println(io,"Primary HO estimates: $(length(primary)); scalar-unstable: $(length(unstable)).")
    println(io,"Numerically certified: $(count(r->r.certified,primary))/$(length(primary)).")
    println(io,"Maximum HO/FD coefficient discrepancy: ", maximum(r.ho_fd_error for r in cert),".")
    println(io,"Maximum FD refinement discrepancy: ", maximum(r.fd_refinement_error for r in cert),".")
    println(io,"Maximum threshold extrapolation discrepancy: ", maximum(r.threshold_error for r in cert),".\n")
    println(io,"| Sector | Threshold routes | U_hg range | Requested emitter formula failures |\n|---|---:|---:|---:|")
    for sector in SECTORS
        rs=filter(r->r.sector==sector && r.q_GeV==0,primary)
        us=[r.U_hg_re for r in rs if r.scalar_stable]
        println(io,"| $sector | $(length(rs)) | ",@sprintf("%.6f–%.6f",minimum(us),maximum(us))," | $(count(r->!r.requested_formula_pass,rs)) |")
    end
    println(io,"\nThese are formal kernels under the existing spectator-coordinate convention. Numerical certification does not certify the disputed emitter-mass extension or physical decay kinematics.")
end
println("Saved figures, plot membership and numerical summary to $OUT")

# Preserve exact numeric tables compactly in git. gzip -n removes timestamp and
# filename headers; the uncompressed local tables are convenient for inspection.
for name in ("traces","comparisons","certificates","waves","wave_samples","plot_membership")
    path=joinpath(OUT,name*".csv")
    isfile(path) && run(pipeline(`gzip -n -c $path`;stdout=path*".gz"))
end

metadata=TOML.parsefile(joinpath(OUT,"run.toml"))
if haskey(metadata,"source_revision")
    metadata["base_revision"]=pop!(metadata,"source_revision")
end
metadata["source_identity"]="working source identified by SHA-256; base_revision precedes study edits"
root=normpath(joinpath(@__DIR__,"../.."))
metadata["source_sha256"]=Dict(p=>bytes2hex(sha256(read(joinpath(root,p)))) for p in
    ("QuarkModelTransitions/src/pseudoscalar_emission.jl","GIPaper/scripts/a_class_kernels.jl",
     "GIPaper/scripts/study_a_class.jl","GIPaper/scripts/plot_a_class.jl"))
open(joinpath(OUT,"run.toml"),"w") do io
    TOML.print(io,metadata)
end
