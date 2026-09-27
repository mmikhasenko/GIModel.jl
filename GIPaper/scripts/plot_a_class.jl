#!/usr/bin/env julia
# julia --project=GIPaper/scripts GIPaper/scripts/plot_a_class.jl
using CSV, CairoMakie, Printf, TOML, SHA, LaTeXStrings
const OUT=normpath(joinpath(@__DIR__,"../docs/a_class_universality"))
const SECTORS=["isovector","isoscalar","strange","charmed","bottom-flavored"]
const COLORS=["#0072B2","#E69F00","#009E73","#CC79A7","#D55E00"]
const LABELS=["³S₁ → ¹S₀ · P wave","³P₂ → ¹S₀ · D wave","³P₂ → ³S₁ · D wave",
              "³P₁ → ³S₁ · D wave","¹P₁ → ³S₁ · D wave","³D₃ → ¹S₀ · F wave","³D₃ → ³S₁ · F wave",
              "¹P₁ → ³P₀ · P wave","³D₃ → ¹P₁ · D wave"]
const XLABEL=L"U_{hg} = \frac{(C_h/C_g)_{\mathrm{native}}}{(C_h/C_g)_{\mathrm{SHO}}}"
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

# Publication target: threshold only, with matching LaTeX subscript typography.
threshold=filter(r->r.q_GeV==0,primary)
threshold_valid=filter(r->r.scalar_stable,threshold)
threshold_lo=min(1.0,minimum(r.U_hg_re for r in threshold_valid))
threshold_hi=max(1.0,maximum(r.U_hg_re for r in threshold_valid))
threshold_pad=max(0.025,0.03*(threshold_hi-threshold_lo))
threshold_bins=collect(range(threshold_lo-threshold_pad,threshold_hi+threshold_pad;length=36))
publication=Figure(size=(800,600),fontsize=23)
publication_axis=Axis(publication[1,1],xlabel=XLABEL,ylabel="Number of combinations",
        xlabelsize=28,ylabelsize=23,xticklabelsize=23,yticklabelsize=23)
publication_counts=zeros(Int,length(threshold_bins)-1,length(SECTORS))
for (s,sector) in enumerate(SECTORS), r in filter(r->r.sector==sector,threshold)
    bin=r.scalar_stable ? clamp(searchsortedlast(threshold_bins,r.U_hg_re),1,length(threshold_bins)-1) : 0
    bin>0 && (publication_counts[bin,s]+=1)
    push!(membership,(figure="threshold_histogram",panel=1,point_id=String(r.point_id),
        bin=bin,x=r.scalar_stable ? r.U_hg_re : r.vector_distance,y=0.0,
        scalar_stable=r.scalar_stable,certified=r.certified,certificate_id=String(r.certificate_id),
        native_g_trace_ids=String(r.native_g_trace_ids),native_h_trace_ids=String(r.native_h_trace_ids),
        sho_g_trace_ids=String(r.sho_g_trace_ids),sho_h_trace_ids=String(r.sho_h_trace_ids)))
end
publication_centers=(threshold_bins[1:end-1]+threshold_bins[2:end])/2
barplot!(publication_axis,repeat(publication_centers,length(SECTORS)),vec(publication_counts);stack=repeat(1:5,inner=length(publication_centers)),
         color=repeat(COLORS,inner=length(publication_centers)),width=0.94*(threshold_bins[2]-threshold_bins[1]))
vlines!(publication_axis,[1.0];color=:black,linestyle=:dash,linewidth=2)
xlims!(publication_axis,first(threshold_bins),last(threshold_bins))
Legend(publication[1,1],[PolyElement(color=c) for c in COLORS],SECTORS;
       halign=:right,valign=:top,tellwidth=false,tellheight=false,
       margin=(12,16,12,16),framevisible=false,labelsize=23)
publication_failed=count(r->!r.certified,threshold)
save(joinpath(OUT,"threshold_histogram.pdf"),publication)
save(joinpath(OUT,"threshold_histogram.png"),publication;px_per_unit=2)
CSV.write(joinpath(OUT,"threshold_bins.csv"),[(bin=i,left=threshold_bins[i],right=threshold_bins[i+1]) for i in 1:length(publication_centers)])
open(joinpath(OUT,"threshold_caption.md"),"w") do io
    println(io,"# Threshold publication figure\n")
    println(io,"**Candidate caption.** Distribution of the threshold ratio \\(U_{hg}= (a_h/a_g)_{\\mathrm{native}}/(a_h/a_g)_{\\mathrm{SHO}}\\), where \\(a_x=\\lim_{q\\to0}C_x(q)/q^L\\). Native HO-basis GI waves are compared with common-beta SHO waves (beta=0.40 GeV), using the same constituent masses, emitter topology and operator conventions. The nine transitions comprise seven A, one A′ and one A″ transition; their equal-mass SHO threshold ratios are +1/4, −1/4 and +1/8, respectively. The dashed line denotes U_hg=1. Colors identify the parent sector. All $(length(threshold)) basis/flavor/emitter combinations enter with unit weight, including charge/isospin-related entries; counts are not independent measurements or decay probabilities. The $publication_failed numerically flagged combinations are retained and individually identified in the trace ledger. This is a formal unmixed-basis kernel comparison, not an on-shell decay sample.")
    println(io,"\n**Numerical status.** Flags require further investigation: HO/FD coefficient disagreement ≥3%, FD refinement disagreement ≥0.3%, or threshold extrapolation disagreement ≥0.2%. See certificates.csv.gz for each diagnostic separately. The existing unequal-mass operator convention remains conditional as discussed in README.md.")
end

fig=Figure(size=(1250,870))
for (i,q) in enumerate([0.0,0.15,0.30,0.60])
    selected=filter(r->r.q_GeV==q,primary)
    title=q==0 ? "Threshold: qᴸ-reduced extrapolation" : "Common q = $q GeV"
    ax=Axis(fig[(i-1)÷2+1,(i-1)%2+1],title=title,xlabel=XLABEL,ylabel="Basis-route count")
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
    text!(ax,0.98,0.97;text="$(length(selected)) combinations · $failed numerical flags · $unstable unstable",
          space=:relative,align=(:right,:top),fontsize=12)
end
Legend(fig[3,1:2],[PolyElement(color=c) for c in COLORS],SECTORS;orientation=:horizontal,framevisible=false)
Label(fig[0,1:2],"A, A′, A″ kernels · finite-q validation",fontsize=25)
Label(fig[4,1:2],"Formal unmixed basis survey; counts are correlated flavor routes, not independent measurements.\nExisting Eq. (19) routing; unequal-mass emitter-formula mismatch remains unresolved.",fontsize=14)
save(joinpath(OUT,"histograms.pdf"),fig); save(joinpath(OUT,"histograms.png"),fig;px_per_unit=2)

panels=Figure(size=(1550,1250))
for (i,label) in enumerate(LABELS)
    ax=Axis(panels[(i-1)÷3+1,(i-1)%3+1],title="$(i==8 ? "A′" : i==9 ? "A″" : "A$i"): $label",xlabel="q [GeV] (0 = extrapolation)",ylabel=L"U_{hg}")
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
Legend(panels[4,1:3],[PolyElement(color=c) for c in COLORS],SECTORS;orientation=:horizontal,framevisible=false)
Label(panels[0,1:3],"Nine transitions · momentum-dependence validation",fontsize=26)
Label(panels[5,1:3],"Lines connect the same flavor/topology kernel. Crosses fail numerical certification.\nCoincident charge/isospin routes remain separate entries in plot_membership.csv.",fontsize=16)
save(joinpath(OUT,"transition_panels.pdf"),panels); save(joinpath(OUT,"transition_panels.png"),panels;px_per_unit=2)

# Unstable ratios are never silently dropped: their normalized-vector distance
# has its own figure and membership. Empty in a well-conditioned run.
unstable=filter(r->!r.scalar_stable,primary)
if !isempty(unstable)
    f=Figure(size=(900,500)); fallback_axis=Axis(f[1,1],xlabel="Route index",ylabel="Phase-aligned vector distance")
    for (i,r) in enumerate(unstable)
        r.vector_valid && scatter!(fallback_axis,[i],[r.vector_distance];color=COLORS[findfirst(==(r.sector),SECTORS)])
    end
    save(joinpath(OUT,"unstable_vectors.pdf"),f)
end
CSV.write(joinpath(OUT,"plot_membership.csv"),membership)
cert=collect(readtable("certificates"))
open(joinpath(OUT,"numerical_summary.md"),"w") do io
    println(io,"# Generated numerical findings\n")
    println(io,"Threshold publication estimates: $(length(threshold)); numerically certified: $(count(r->r.certified,threshold))/$(length(threshold)).\n")
    println(io,"Primary HO estimates: $(length(primary)); scalar-unstable: $(length(unstable)).")
    println(io,"Numerically certified: $(count(r->r.certified,primary))/$(length(primary)).")
    println(io,"Maximum HO/FD coefficient discrepancy: ", maximum(r.ho_fd_error for r in cert),".")
    println(io,"Maximum FD refinement discrepancy: ", maximum(r.fd_refinement_error for r in cert),".")
    println(io,"Maximum threshold extrapolation discrepancy: ", maximum(r.threshold_error for r in cert),".\n")
    println(io,"| Sector | Threshold routes | U_hg range | A-only emitter formula failures |\n|---|---:|---:|---:|")
    for sector in SECTORS
        rs=filter(r->r.sector==sector && r.q_GeV==0,primary)
        us=[r.U_hg_re for r in rs if r.scalar_stable]
        println(io,"| $sector | $(length(rs)) | ",@sprintf("%.6f–%.6f",minimum(us),maximum(us))," | $(count(r->r.mass_formula_applicable && r.requested_formula_pass===false,rs)) |")
    end
    println(io,"\n| Family | Threshold combinations | U_hg range | Numerical flags |\n|---|---:|---:|---:|")
    for family in ("A","Aprime","Adoubleprime")
        rs=filter(r->r.family==family,threshold)
        us=[r.U_hg_re for r in rs if r.scalar_stable]
        println(io,"| $family | $(length(rs)) | ",@sprintf("%.6f–%.6f",minimum(us),maximum(us))," | $(count(r->!r.certified,rs)) |")
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
# A render of archived traces must preserve the calculation's source identity.
# Only the renderer is new; it has not rerun the numerical kernels.
metadata["render_source_sha256"]=bytes2hex(sha256(read(@__FILE__)))
open(joinpath(OUT,"run.toml"),"w") do io
    TOML.print(io,metadata)
end
