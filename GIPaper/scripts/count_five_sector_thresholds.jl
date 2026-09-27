#!/usr/bin/env julia
# julia --project=GIPaper/scripts GIPaper/scripts/count_five_sector_thresholds.jl
using TOML, CSV, SHA, Printf
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const OUT=joinpath(ROOT,"GIPaper/docs/pseudoscalar_census")
const SECTORS=["isovector","isoscalar","strange","charmed","bottom-flavored"]
const S0=(1,0,0,0); const S1=(1,0,1,1); const P0=(1,1,1,0)
const P1=(1,1,1,1); const B1=(1,1,0,1); const P2=(1,1,1,2)
const D1=(1,2,1,1); const D3=(1,2,1,3)
const TEMPLATES=[
 ("T01",S1,S0,[1]),("T02",P2,S0,[2]),("T03",P2,S1,[2]),
 ("T04",P1,S1,[0,2]),("T05",B1,S1,[0,2]),("T06",D3,S0,[3]),
 ("T07",D3,S1,[3]),("T08",B1,P0,[1]),("T09",D3,B1,[2]),
 ("T10",P1,P0,[1]),("T11",P0,S0,[0]),("T12",D1,S0,[1]),
 ("T13",D1,S1,[1]),("T14",(2,0,0,0),S1,[1]),
 ("T15",(2,0,1,1),S0,[1]),("T16",(2,0,1,1),S1,[1])]
spec(x)="$(x[1])^$(2x[3]+1)$("SPDFG"[x[2]+1])_$(x[4])"
state(name,sector,multiplet,conjugate,masssector,coeff)=
    (;name,sector,multiplet,conjugate,masssector,coeff=Dict(coeff))
function flavors(;ground=false)
    a=[state("pi+","isovector","pi","pi-","nn",[("u","d")=>1]),
       state("pi0","isovector","pi","pi0","nn",[("u","u")=>1,("d","d")=>-1]),
       state("pi-","isovector","pi","pi+","nn",[("d","u")=>1]),
       state("eta","isoscalar","eta","eta","nn",ground ? [("u","u")=>1,("d","d")=>1,("s","s")=>1] : [("u","u")=>1,("d","d")=>1]),
       state("eta-prime","isoscalar","eta-prime","eta-prime","ss",ground ? [("u","u")=>1,("d","d")=>1,("s","s")=>-2] : [("s","s")=>1])]
    for (name,mp,cc,pair) in [("K+","K","K-",("u","s")),("K0","K","K0bar",("d","s")),
        ("K-","Kbar","K+",("s","u")),("K0bar","Kbar","K0",("s","d"))]
        push!(a,state(name,"strange",mp,cc,"K",[pair=>1]))
    end
    for (heavy,sector,base) in [("c","charmed","D"),("b","bottom-flavored","B")], light in ["u","d","s"], anti in [false,true]
        pair=anti ? (light,heavy) : (heavy,light)
        conjugate=join(reverse(pair))*"bar"
        m=base*(light=="s" ? "s" : "")
        push!(a,state(join(pair)*"bar",sector,m*(anti ? "_anti" : ""),conjugate,m,[pair=>1]))
    end
    a
end
const BASIS=flavors(); const GROUND=flavors(ground=true); const FIELDS=GROUND[1:9]
const BYNAME=Dict(s.name=>s for s in BASIS)
const Q3=Dict("u"=>2,"d"=>-1,"s"=>-1,"c"=>2,"b"=>-1)
charge(s)=only(unique([Q3[a]-Q3[b] for (a,b) in keys(s.coeff)]))
mc(f)=f in ("u","d") ? "n" : f
function blocks(i,f,p,sign)
    b=Dict{Tuple{String,String,String},Int}()
    for ((a,c),ci) in i.coeff, z in ("u","d","s")
        if a in ("u","d","s")
            k=(mc(a),mc(c),mc(z));b[k]=get(b,k,0)+ci*get(f.coeff,(z,c),0)*get(p.coeff,(a,z),0)
        end
        if c in ("u","d","s")
            k=(mc(c),mc(a),mc(z));b[k]=get(b,k,0)+sign*ci*get(f.coeff,(a,z),0)*get(p.coeff,(z,c),0)
        end
    end
    filter!(kv->last(kv)!=0,b)
end
const PDGPATH=joinpath(ROOT,"GIPaper/data/mass_inputs/pdg-2026.csv")
const PDG=Dict(String(r.key)=>r.mass_GeV for r in CSV.File(PDGPATH))
const EMASS=Dict(e=>PDG[k] for (e,k) in ["pi0"=>"pi0","pi+"=>"pi+","pi-"=>"pi+","eta"=>"eta","eta-prime"=>"eta_prime", "K+"=>"K+","K-"=>"K+","K0"=>"K0","K0bar"=>"K0"])
function masses(path)
    doc=TOML.parsefile(path)
    Dict((r["sector"],(r["n"],findfirst(==(only(r["L"])),"SPDFG")-1,(r["multiplicity"]-1)÷2,r["J"]))=>r["unmixed_GeV"] for r in doc["states"])
end
const FINEPATH=joinpath(OUT,"five_sector_masses_refined.toml")
const COARSEPATH=joinpath(OUT,"five_sector_masses.toml")
const FINE=masses(FINEPATH); const COARSE=masses(COARSEPATH)
mass(t,s,m)=t==S0 && haskey(EMASS,s.name) ? EMASS[s.name] : m[(s.masssector,t)]
function key(ti,tf,p,f,e;ell=nothing,iso=false,cc=false)
    names=String[]
    for name in (p,f,e)
        s=BYNAME[cc ? BYNAME[name].conjugate : name]
        push!(names,iso ? s.multiplet : s.name)
    end
    tf==S0 && sort!(view(names,2:3))
    join(vcat([spec(ti),spec(tf)],isnothing(ell) ? String[] : [string(ell)],names),"|")
end
function display(t,name)
    if t==S0 && haskey(EMASS,name);return name;end
    names=Dict("cubar"=>"D0","cdbar"=>"D+","csbar"=>"Ds+","ucbar"=>"D0bar","dcbar"=>"D-","scbar"=>"Ds-",
        "bubar"=>"B-","bdbar"=>"B0bar","bsbar"=>"Bs0bar","ubbar"=>"B+","dbbar"=>"B0","sbbar"=>"Bs0")
    base=name=="eta" ? "nn(I=0)" : name=="eta-prime" ? "ss(I=0)" : get(names,name,name)
    "$(base)[$(spec(t))]"
end
function main()
    channels=Dict{String,Dict}();ordered=Dict(s=>0 for s in SECTORS)
    for (tid,ti,tf,ells) in TEMPLATES
        @assert all(abs(ti[4]-tf[4])<=l<=ti[4]+tf[4] && isodd(ti[2]-tf[2]-l) for l in ells)
        sign=(-1)^(ti[2]+ti[3]+tf[2]+tf[3])
        for i in BASIS, f in (tf==S0 ? GROUND : BASIS), e in FIELDS
            b=blocks(i,f,e,sign);isempty(b) && continue
            @assert charge(i)==charge(f)+charge(e)
            ordered[i.sector]+=length(ells)
            mi=mass(ti,i,FINE);mf=mass(tf,f,FINE);me=EMASS[e.name]
            delta=mi-mf-me
            cd=mass(ti,i,COARSE)-mass(tf,f,COARSE)-me
            k=key(ti,tf,i.name,f.name,e.name)
            isokey=key(ti,tf,i.name,f.name,e.name;iso=true)
            cckey=min(isokey,key(ti,tf,i.name,f.name,e.name;iso=true,cc=true))
            conjugate=key(ti,tf,i.name,f.name,e.name;cc=true)
            if haskey(channels,k)
                r=channels[k]
                @assert isapprox(r["excess_GeV"],delta;atol=1e-12)
                r["waves"]=sort!(unique(vcat(r["waves"],ells)))
                push!(r["roles"],f.name*" + "*e.name)
            else
                channels[k]=Dict("key"=>k,"sector"=>i.sector,"template"=>tid,
                    "parent"=>display(ti,i.name),"daughter"=>display(tf,f.name),"emitted"=>e.name,
                    "parent_GeV"=>mi,"daughter_GeV"=>mf,"emitted_GeV"=>me,
                    "excess_GeV"=>delta,"coarse_excess_GeV"=>cd,"open"=>delta>0,
                    "waves"=>copy(ells),"isospin_key"=>isokey,"isospin_C_key"=>cckey,
                    "conjugate_key"=>conjugate,"roles"=>[f.name*" + "*e.name])
            end
        end
    end
    for r in values(channels)
        c=channels[r["conjugate_key"]]
        @assert isapprox(r["excess_GeV"],c["excess_GeV"];atol=1e-12)
        @assert r["open"]==c["open"]
    end
    rows=sort!(collect(values(channels));by=r->(findfirst(==(r["sector"]),SECTORS),r["key"]))
    summary=Dict[]
    for sector in SECTORS
        allrows=filter(r->r["sector"]==sector,rows);good=filter(r->r["open"],allrows)
        iso=unique(r["isospin_key"] for r in good);isoC=unique(r["isospin_C_key"] for r in good)
        partial=length(unique((r["isospin_key"],l) for r in good for l in r["waves"]))
        partialC=length(unique((r["isospin_C_key"],l) for r in good for l in r["waves"]))
        split=count(k->any(!r["open"] for r in allrows if r["isospin_key"]==k),iso)
        push!(summary,Dict("sector"=>sector,"before_channels"=>length(allrows),
            "before_partial_waves"=>sum(length(r["waves"]) for r in allrows),"ordered_partial_waves"=>ordered[sector],
            "open_channels"=>length(good),"open_partial_waves"=>sum(length(r["waves"]) for r in good),
            "isospin_channels"=>length(iso),"isospin_partial_waves"=>partial,
            "isospin_C_channels"=>length(isoC),"isospin_C_partial_waves"=>partialC,
            "partly_open_isospin_groups"=>split,
            "grid_flips"=>count(r->r["open"]!=(r["coarse_excess_GeV"]>0),allrows)))
    end
    charm=only(filter(r->r["sector"]=="charmed",summary))
    # The algebraic inventory is fixed; kinematically open counts must follow
    # the current computed masses rather than a pre-correction numerical pin.
    @assert (charm["before_channels"],charm["before_partial_waves"])==(448,504)
    ledger=joinpath(OUT,"five_sector_thresholds.toml")
    open(ledger,"w") do io
        TOML.print(io,Dict("scope"=>"18 Table-IV partial waves; unmixed spin and ideal excited-isoscalar flavor basis",
            "mass_rule"=>"physical ground light pseudoscalars in either final role, all other masses calculated",
            "group_rule"=>"isospin group open if at least one charge member is open; charge conjugation separately reduced",
            "source_hashes"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in (FINEPATH,COARSEPATH,PDGPATH,@__FILE__)),
            "summary"=>summary,"channels"=>rows))
    end
    Base.run(pipeline(`gzip -n -c $ledger`;stdout=ledger*".gz"));rm(ledger)
    write_report(rows,summary)
    for r in summary;println(r);end
end
function write_report(rows,summary)
    open(joinpath(OUT,"five_sector_scratchpad.md"),"w") do io
        println(io,"# Five-sector threshold and isospin census\n")
        println(io,"Generated by `count_five_sector_thresholds.jl`. This is the same 18 partial waves / 16 basis transition templates audited for charm, not an exhaustive spectrum-wide decay list.\n")
        println(io,"## Conventions\n")
        println(io,"- Count a channel only if M(parent) > M(daughter) + m(P), using unrounded masses. No widths or finite-width tails are calculated.")
        println(io,"- All ground light pseudoscalars (π, K, η, η′) use the stored physical masses in either final role. Other states use 900-point calculated fixed-sector masses; 450 points provide a refinement check.")
        println(io,"- Excited isoscalars are separate ideal nn and ss basis states with their own model masses. The ground η/η′ flavor states contain both nn and ss. No arbitrary mixed flavor state is assigned a bare mass. This changes the scope from the earlier generic-isoscalar census.")
        println(io,"- Spin states remain unmixed, including the two axial basis states. Thus these counts are model/basis candidates, not a census of measured resonances. The charm subset is checked independently by the Python threshold audit.")
        println(io,"- Quark/antiquark contributions are combined before counting. A final pair of ground pseudoscalars is unordered; S and D waves of the same final channel count once in the channel column and twice in the partial-wave column.")
        println(io,"- Isospin reduction replaces each charged member by its multiplet: π triplet, K doublet, anti-K doublet, D doublet, anti-D doublet, etc. Ds and Bs are separate singlets; nn and ss are distinct isoscalars. It groups only charge variants of the SAME spectroscopic transition and final multiplets.")
        println(io,"- Charge conjugates stay separate in the isospin-only columns. The additional I+C columns merge them explicitly. A reduced group survives if ANY charge member is open; partially open groups are listed separately. These group counts are not counts of independent coupling constants.\n")
        println(io,"## Final-channel counts\n\n| Parent sector | Before cut | Open charge-resolved | Isospin reduced | Isospin + C reduced | Partly open I groups |\n|---|---:|---:|---:|---:|---:|")
        for r in summary;println(io,"| $(r["sector"]) | $(r["before_channels"]) | $(r["open_channels"]) | $(r["isospin_channels"]) | $(r["isospin_C_channels"]) | $(r["partly_open_isospin_groups"]) |");end
        println(io,"| Total | ",join([sum(r[k] for r in summary) for k in ("before_channels","open_channels","isospin_channels","isospin_C_channels","partly_open_isospin_groups")]," | ")," |")
        println(io,"\n## Partial-wave counts\n\n| Parent sector | Before cut | Open charge-resolved | Isospin reduced | Isospin + C reduced | Grid flips (final channels) |\n|---|---:|---:|---:|---:|---:|")
        for r in summary;println(io,"| $(r["sector"]) | $(r["before_partial_waves"]) | $(r["open_partial_waves"]) | $(r["isospin_partial_waves"]) | $(r["isospin_C_partial_waves"]) | $(r["grid_flips"]) |");end
        println(io,"| Total | ",join([sum(r[k] for r in summary) for k in ("before_partial_waves","open_partial_waves","isospin_partial_waves","isospin_C_partial_waves","grid_flips")]," | ")," |")
        println(io,"\n## Counts by transition template\n\nEach cell is open charge channels / isospin groups. The two axial -> vector rows each carry S and D partial waves.\n\n| Template | Initial → final | isovector | isoscalar | strange | charmed | bottom-flavored |\n|---|---|---:|---:|---:|---:|---:|")
        for (tid,ti,tf,ells) in TEMPLATES
            cells=String[]
            for sector in SECTORS
                good=filter(r->r["sector"]==sector && r["template"]==tid && r["open"],rows)
                push!(cells,"$(length(good)) / $(length(unique(r["isospin_key"] for r in good)))")
            end
            println(io,"| $tid | $(spec(ti)) → $(spec(tf)) | ",join(cells," | ")," |")
        end
        println(io,"\n## Mass and grouping audit\n\nFull masses and all rejected as well as retained channels are recorded in `five_sector_masses*.toml` and `five_sector_thresholds.toml.gz`. The latter includes the exact grouping key, the emission-role assignments, both-grid threshold excesses and source hashes.\n")
        near=filter(r->abs(r["excess_GeV"])<0.010 || r["open"]!=(r["coarse_excess_GeV"]>0),rows)
        println(io,"### Within 10 MeV of threshold, or changed under refinement\n\n| Sector | Channel | Excess (900), MeV | Excess (450), MeV | Open |\n|---|---|---:|---:|---|")
        for r in near;println(io,"| $(r["sector"]) | $(r["parent"]) → $(r["daughter"]) + $(r["emitted"]) | ",@sprintf("%+.6f | %+.6f",1000r["excess_GeV"],1000r["coarse_excess_GeV"])," | $(r["open"]) |");end
        println(io,"\n## Isospin groups: every retained group and its charge members\n")
        for sector in SECTORS
            println(io,"### $sector\n")
            sr=filter(r->r["sector"]==sector,rows)
            for k in sort!(unique([r["isospin_key"] for r in sr if r["open"]]))
                members=filter(r->r["isospin_key"]==k,sr);open_count=count(r->r["open"],members)
                println(io,"- **`$k`**: $open_count/$(length(members)) charge channels open.")
                for r in members
                    println(io,"  - $(r["open"] ? "OPEN" : "CLOSED") $(r["parent"]) → $(r["daughter"]) + $(r["emitted"]); ℓ=$(join(r["waves"],",")); Δ=",@sprintf("%+.3f",1000r["excess_GeV"])," MeV.")
                end
            end
        end
    end
end
main()
