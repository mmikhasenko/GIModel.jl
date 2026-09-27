#!/usr/bin/env julia
using GIModel, CSV, DataFrames, LinearAlgebra, Statistics, Printf, CairoMakie

root = dirname(@__DIR__)
out = joinpath(root, "docs/investigations")
mkpath(joinpath(out, "figures"))
# The plot's recorded identities define scope; no second spectrum envelope.
scope = DataFrame(CSV.File(joinpath(@__DIR__, "spectrum_plots/ten_meson_calculated_spectrum.csv")))
params, mq = load_parameters_and_quark_masses(default_parameters_path())
solver = FiniteDifferenceSolver(ngrid=450, rmax=24.0)
blocks, shifts, angles = NamedTuple[], NamedTuple[], NamedTuple[]
specs = Dict()
if "--reuse" in ARGS
    append!(blocks, NamedTuple.(eachrow(DataFrame(CSV.File(joinpath(out,"mixing_blocks.csv"))))))
    append!(shifts, NamedTuple.(eachrow(DataFrame(CSV.File(joinpath(out,"mixing_state_shifts.csv"))))))
    append!(angles, NamedTuple.(eachrow(DataFrame(CSV.File(joinpath(out,"mixing_projected_angles.csv"))))))
else
for group in groupby(scope, [:flavor1, :flavor2])
    f1, f2 = Symbol(group.flavor1[1]), Symbol(group.flavor2[1])
    levels = unique(b -> (b.n,b.L_label,b.multiplicity,b.J),
        [BasisState(r.n, r.L, r.multiplicity, r.J) for r in eachrow(group)])
    println("Solving ", f1, " ", f2, " (", length(levels), " levels)"); flush(stdout)
    spec = compute_spectrum(params, Meson(mq,f1,f2); levels, solver)
    specs[(f1,f2)] = spec
    seen = IdDict()
    sector = string(f1,"/",f2)
    for s in spec.states, m in s.mixings
        push!(shifts, (sector, state=s.label, mechanism=m.mechanism,
            before_MeV=1000m.unmixed_GeV, after_MeV=1000m.result.masses[m.eigenstate],
            shift_MeV=1000(m.result.masses[m.eigenstate]-m.unmixed_GeV),
            impurity=1-maximum(abs2,m.components)))
        haskey(seen,m.result) && continue
        seen[m.result] = true
        r = m.result; H = r.block.matrix; U = r.vectors; e = r.masses
        residual = opnorm(H*U-U*Diagonal(e),Inf)
        orth = opnorm(U'U-I,Inf)
        trace_error = abs(sum(e)-tr(H))
        spectral_error = abs(sum(abs2,e)-sum(abs2,H))
        @assert residual < 1e-10 && orth < 1e-10 && trace_error < 1e-10 && spectral_error < 1e-9
        # A fixed basis rephasing changes no mass or component probability.
        D = Diagonal([isodd(i) ? -1.0 : 1.0 for i in axes(H,1)])
        @assert maximum(abs.(eigvals(Symmetric(D*H*D))-e)) < 1e-10
        push!(blocks,(sector,mechanism=m.mechanism,dimension=length(e),
            basis=join([b.label for b in r.block.basis],";"),residual,orth,trace_error,spectral_error,
            max_shift_MeV=1000maximum(abs.(e-sort(diag(H)))),
            lower_shift_MeV=1000(first(e)-minimum(diag(H))),
            upper_shift_MeV=1000(last(e)-maximum(diag(H)))))
        @assert first(e) <= minimum(diag(H))+1e-12 && last(e) >= maximum(diag(H))-1e-12
    end
    for L in ("P","D","F","G"), n in 1:2
        ss = filter(s -> s.n == n && s.L == L && s.multiplicity == 1, spec.states)
        isempty(ss) && continue
        s = only(ss)
        ms = filter(m -> m.mechanism == "antisymmetric_spin_orbit",s.mixings)
        isempty(ms) && continue
        m = only(ms); r = m.result
        i = findfirst(b -> b.n == n && b.multiplicity == 1,r.block.basis)
        j = findfirst(b -> b.n == n && b.multiplicity == 3,r.block.basis)
        t = spectrum_state(spec,n,L,3,s.J)
        mt = only(filter(m -> m.mechanism == "antisymmetric_spin_orbit",t.mixings))
        col = s.mass_GeV < t.mass_GeV ? m.eigenstate : mt.eigenstate
        cs,ct = r.vectors[i,col],r.vectors[j,col]
        phase = cs < 0 ? -1 : 1
        theta = atand(phase*ct,phase*cs)
        H = r.block.matrix[[i,j],[i,j]]
        gap = abs(H[1,1]-H[2,2]); c=H[1,2]
        split = hypot(gap,2c)
        repulsion = (split-gap)/2
        exact = eigvals(Symmetric(H))
        @assert maximum(abs.(exact-[(tr(H)-split)/2,(tr(H)+split)/2])) < 1e-12
        push!(angles,(sector,nL=string(n,L),theta,
            singlet_low=cs^2/(cs^2+ct^2),outside=max(0.0,1-cs^2-ct^2),
            gap_MeV=1000gap,c_MeV=1000c,pair_repulsion_MeV=1000repulsion,
            lower_shift_MeV=1000(min(s.mass_GeV,t.mass_GeV)-minimum(diag(H))),
            upper_shift_MeV=1000(max(s.mass_GeV,t.mass_GeV)-maximum(diag(H))),
            low_MeV=1000min(s.mass_GeV,t.mass_GeV),high_MeV=1000max(s.mass_GeV,t.mass_GeV)))
    end
end
CSV.write(joinpath(out,"mixing_blocks.csv"),blocks)
CSV.write(joinpath(out,"mixing_state_shifts.csv"),shifts)
CSV.write(joinpath(out,"mixing_projected_angles.csv"),angles)
end

# Reuse the existing paper-angle ledger rather than transcribing targets again.
paper = NamedTuple[]
sector_map = Dict("u sbar"=>"q/s","c ubar"=>"c/q","c sbar"=>"c/s",
    "b ubar"=>"b/q","b sbar"=>"b/s","b cbar"=>"b/c")
refmap = Dict("q/s"=>("strange","strange"),"c/q"=>("charmed","charmed"),
    "c/s"=>("charmed","charmed_strange"),"b/q"=>("b_flavored","bottom_light"),
    "b/s"=>("b_flavored","bottom_strange"),"b/c"=>("b_flavored","bottom_charm"))
for line in eachline(joinpath(root,"docs/residual_reports/mixing_angles.md"))
    cells = strip.(split(line,'|'))
    length(cells)==7 && haskey(sector_map,cells[2]) || continue
    sec = sector_map[cells[2]]; nL=cells[3]; target=parse(Float64,cells[5])
    a = only(filter(a -> a.sector==sec && a.nL==nL,angles))
    file,refsector=refmap[sec]
    refs=DataFrame(CSV.File(joinpath(root,"data/reference_spectrum_"*file*".csv")))
    candidates=filter(r -> r.sector==refsector && r.n==parse(Int,nL[1:1]) &&
        r.L==nL[2:2] && r.J==findfirst(==(nL[2:2]),["S","P","D","F","G"])-1,refs)
    @assert nrow(candidates) in (0,2) "Incomplete reference pair for $sec $nL"
    lo,hi=nrow(candidates)==2 ? sort(1000 .* candidates.mass_GeV) : (NaN,NaN)
    # Conditional inverse 2x2 reconstruction, not published unmixed GI masses.
    gi_gap=(hi-lo)*abs(cosd(2target))
    gi_repulsion=((hi-lo)-gi_gap)/2
    push!(paper,merge(a,(paper_theta=target,angle_distance=abs(mod(a.theta-target+90,180)-90),
        paper_low_MeV=lo,paper_high_MeV=hi,
        paper_inferred_repulsion_MeV=gi_repulsion,
        paper_inferred_c_MeV=-(hi-lo)*sind(2target)/2)))
end
@assert length(paper)==13
CSV.write(joinpath(out,"mixing_paper_comparison.csv"),paper)

fig=Figure(size=(1250,950),fontsize=15)
ax=Axis(fig[1,1],xlabel="Angle (degrees)",
    yticks=(1:13,[r.sector*" "*r.nL for r in paper]),yreversed=true,
    title="13 quoted angles: fixed paper convention")
for (i,r) in enumerate(paper)
    lines!(ax,[r.paper_theta,r.theta],[i,i],color=:gray)
end
scatter!(ax,[r.paper_theta for r in paper],1:13,color=:darkorange,markersize=10,label="GI")
scatter!(ax,[r.theta for r in paper],1:13,color=:steelblue,markersize=10,label="Computed")
Legend(fig[0,1],ax,orientation=:horizontal,labelsize=12,tellwidth=false)
ax2=Axis(fig[1,2],xlabel="Lower-state singlet probability",ylabel="Number of projected pairs",
    title="All computed same-n spin–orbit pairs")
hist!(ax2,[r.singlet_low for r in angles],bins=0:0.1:1,color=:steelblue)
ax3=Axis(fig[2,1],xlabel="GI-inferred two-level repulsion (MeV)",ylabel="Computed two-level repulsion (MeV)",
    title="Shift magnitude is not generally reproduced")
comparable=filter(r->isfinite(r.paper_inferred_repulsion_MeV),paper)
limit=maximum(r.paper_inferred_repulsion_MeV for r in comparable)*1.1
lines!(ax3,[0,limit],[0,limit],color=:gray,linestyle=:dash)
scatter!(ax3,[r.paper_inferred_repulsion_MeV for r in comparable],[r.pair_repulsion_MeV for r in comparable],color=:darkorange,markersize=12)
for r in comparable
    offset = r.sector=="q/s" && r.nL=="1F" ? (-45,10) :
        r.sector=="q/s" && r.nL=="1G" ? (5,-15) : (5,3)
    text!(ax3,r.paper_inferred_repulsion_MeV,r.pair_repulsion_MeV,text=r.sector*" "*r.nL,fontsize=10,offset=offset)
end
ax4=Axis(fig[2,2],xlabel="Mixing-only mass shift (MeV)",ylabel="Number of states",
    title="Full blocks; each physical state counted once")
for (mech,color) in (("tensor_mixing",(:steelblue,0.65)),("antisymmetric_spin_orbit",(:darkorange,0.55)))
    hist!(ax4,[r.shift_MeV for r in shifts if r.mechanism==mech],bins=-40:2:40,color=color,label=mech)
end
axislegend(ax4,position=:rt,labelsize=11)
Label(fig[3,1:2],"FD 450, rmax=24 GeV⁻¹ • ten-panel envelope • isoscalar annihilation excluded by plot pipeline\nGI shifts inferred from rounded physical masses + angles; not separately published unmixed masses",fontsize=13)
save(joinpath(out,"figures/mixing_layer_review.png"),fig)
save(joinpath(out,"figures/mixing_layer_review.svg"),fig)

open(joinpath(out,"mixing_layer_results.md"),"w") do io
    println(io,"# Mixing-layer numerical census\n\nGenerated by `GIPaper/scripts/investigate_mixing_layer.jl`.\n")
    println(io,"Unique flavor sectors: $(length(unique(r.sector for r in blocks))); blocks: $(length(blocks)); participating states: $(length(shifts)); projected same-n spin–orbit angles: $(length(angles)); paper angles: $(length(paper)).\n")
    println(io,"| Sector | SO blocks | Tensor blocks | Participating states | Projected SO angles |\n|---|---:|---:|---:|---:|")
    for sec in unique(r.sector for r in blocks)
        nso=count(r->r.sector==sec && r.mechanism=="antisymmetric_spin_orbit",blocks)
        nt=count(r->r.sector==sec && r.mechanism=="tensor_mixing",blocks)
        println(io,"| $sec | $nso | $nt | $(count(r->r.sector==sec,shifts)) | $(count(r->r.sector==sec,angles)) |")
    end
    println(io,"\nMaximum eigen-equation residual: $(maximum(r.residual for r in blocks)) GeV; orthogonality defect: $(maximum(r.orth for r in blocks)); trace defect: $(maximum(r.trace_error for r in blocks)) GeV. All blocks pass basis-rephasing invariance and extreme-eigenvalue repulsion checks.\n")
    println(io,"| Mechanism | States | Median absolute shift MeV | Maximum absolute shift MeV | Maximum impurity |\n|---|---:|---:|---:|---:|")
    for mechanism in ("tensor_mixing","antisymmetric_spin_orbit")
        rows=filter(r->r.mechanism==mechanism,shifts)
        @printf(io,"| %s | %d | %.4f | %.4f | %.4f |\n",mechanism,length(rows),median(abs(r.shift_MeV) for r in rows),maximum(abs(r.shift_MeV) for r in rows),maximum(r.impurity for r in rows))
    end
    println(io,"| Sector | nL | angle | GI | lower shift MeV | upper shift MeV | 2×2 repulsion MeV | GI-inferred repulsion MeV | outside pair % |\n|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in paper
        @printf(io,"| %s | %s | %.2f | %.0f | %+.3f | %+.3f | %.3f | %.3f | %.4f |\n",r.sector,r.nL,r.theta,r.paper_theta,r.lower_shift_MeV,r.upper_shift_MeV,r.pair_repulsion_MeV,r.paper_inferred_repulsion_MeV,100max(0,r.outside))
    end
    println(io,"\nAngular distance is modulo 180°, with no low/high exchange. Median: $(median(r.angle_distance for r in paper)) degrees. Full matrices are multi-level: projected angles are not independent block rotation parameters. The CSV files retain every block, every participating state's signed shift, every projected angle and all 13 paper comparisons.\n")
    println(io,"Only $(length(comparable))/13 angle rows have both physical masses in the canonical numeric spectrum catalog. NaN denotes an unavailable inferred shift, not zero; no different-J mass is substituted.")
end
println("Wrote mixing census and figure to ",out)
