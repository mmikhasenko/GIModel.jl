#!/usr/bin/env julia
# Precompute the radial amplitudes consumed by the paper renderer.
# Including this file defines the API without solving or writing anything.
include(joinpath(@__DIR__, "density_3d.jl"))
using TOML, SHA

const CANDIDATE_DIR = joinpath(@__DIR__, "data", "density_candidates")
const CANDIDATES = [
    (id="rho", name="ρ", sector=:light, label="1^3S_1"),
    (id="psi2s", name="ψ(2S)", sector=:charm, label="2^3S_1"),
    (id="upsilon3s", name="Υ(3S)", sector=:bottom, label="3^3S_1"),
]

"""
    save_density_state(spec, label; id, name=label, outdir=CANDIDATE_DIR,
                       rgrid, flavors=String[], provenance=Dict())

Export any state from a GIModel spectrum for later rendering. Signed amplitudes
are kept by (L,S), preserving radial and angular interference. Supply provenance
for custom solves; no model parameters are inferred from a label.
"""
function save_density_state(spec, label; id, name=label, outdir=CANDIDATE_DIR,
                            rgrid=collect(range(0.0,RMAX_INV_GEV;length=NRADIAL)),
                            flavors=String[], provenance=Dict())
    occursin(r"^[a-zA-Z0-9_-]+$",id) || throw(ArgumentError("id must be a plain filename stem"))
    grid = Float64.(rgrid)
    length(grid) >= 3 || throw(ArgumentError("rgrid needs at least three points"))
    h = grid[2]-grid[1]
    h > 0 && all(isapprox.(diff(grid),h;rtol=1e-10,atol=1e-12)) ||
        throw(ArgumentError("rgrid must be uniformly increasing"))
    J, channels = cloud_state(spec,label,grid)
    weights = radial_weights(channels)
    residuals = [check_angular_normalization(channels,J,m,weights;nx=4001) for m in 0:J]
    isotropy = check_mj_sum_isotropy(channels,J)
    @assert maximum(residuals) < 1e-6
    @assert isotropy < 1e-12
    symmetry = maximum(abs(angular_profile(channels,k,J,m,x) -
                           angular_profile(channels,k,J,-m,x))
                       for k in 1:20:length(grid), m in 0:J, x in (-0.8,0.0,0.7))
    @assert symmetry < 1e-12 * maximum(weights)
    data = Dict(
        "schema_version"=>1,"id"=>id,"name"=>name,"label"=>label,"J"=>J,
        "flavors"=>string.(flavors),
        "mass_GeV"=>spectrum_state(spec,label).mass_GeV,
        "r_rms_fm"=>sqrt(sum(weights .* grid.^2)/sum(weights))*HBARC_FM,
        "r_inv_GeV"=>grid,"hbarc_GeV_fm"=>HBARC_FM,
        "channels"=>[Dict("L"=>c.L,"S"=>c.S,"U"=>c.U,"weight"=>channel_weight(c,h)) for c in channels],
        "checks"=>Dict("angular_residuals_m0_to_J"=>residuals,"isotropy"=>isotropy,
                       "opposite_m_absolute_difference"=>symmetry),
        "provenance"=>merge(Dict("julia_version"=>string(VERSION)),provenance))
    mkpath(outdir)
    open(joinpath(outdir,id*".toml"),"w") do io
        TOML.print(io,data;sorted=true)
    end
    @printf("  %-12s M=%.4f GeV  rms=%.4f fm  norm=%.2e\n",id,data["mass_GeV"],data["r_rms_fm"],maximum(residuals))
    flush(stdout)
    return data
end

"""Solve each flavor sector once and save portable, signed radial amplitudes.
No images or Monte-Carlo samples enter this cache. Re-run explicitly after
changing model parameters, solver settings, or the source implementation."""
function precompute_candidates(; outdir=CANDIDATE_DIR)
    params,masses = load_parameters_and_quark_masses(PARAMS_PATH)
    solver = OscillatorSolver()
    source_dir = joinpath(dirname(@__DIR__),"src")
    sources = sort(filter(f->endswith(f,".jl"),readdir(source_dir;join=true)))
    source_hashes = Dict(relpath(f,dirname(@__DIR__))=>bytes2hex(sha256(read(f))) for f in sources)
    sectors = [(:light,:u,:d,2,[("S",3,1),("D",3,1)]),
               (:charm,:c,:c,2,[("S",3,1),("D",3,1)]),
               (:bottom,:b,:b,3,[("S",3,1),("D",3,1)])]
    for (sector,q,qb,nmax,keys) in sectors
        println("Solving $sector ..."); flush(stdout)
        levels = filter(b->(b.L_label,b.multiplicity,b.J) in keys,
                        spectrum_levels(nmax;L_labels=unique(first.(keys))))
        spec = compute_spectrum(params,Meson(masses,q,qb);levels,solver)
        provenance = Dict(
            "parameters"=>read(PARAMS_PATH,String),"parameters_sha256"=>bytes2hex(sha256(read(PARAMS_PATH))),
            "manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"Manifest.toml")))),
            "source_sha256"=>source_hashes,
            "solver"=>Dict(string(f)=>getfield(solver,f) for f in fieldnames(typeof(solver))),
            "nmax"=>nmax,"channels_in_solve"=>["$(k[1]), S=$((k[2]-1)÷2), J=$(k[3])" for k in keys],
            "status"=>"Provisional model parameters; names identify target states, masses are model predictions.")
        for c in filter(c->c.sector==sector,CANDIDATES)
            save_density_state(spec,c.label;id=c.id,name=c.name,flavors=[q,qb],provenance,outdir)
        end
    end
end

"""Load a saved state; rendering needs no spectrum solve."""
function load_candidate(id; datadir=CANDIDATE_DIR)
    d = TOML.parsefile(joinpath(datadir,id*".toml"))
    d["schema_version"] == 1 || error("Unsupported density cache schema")
    channels = [CloudChannel(c["L"],c["S"],Float64.(c["U"])) for c in d["channels"]]
    return (; data=d, channels, J=d["J"], rgrid=Float64.(d["r_inv_GeV"]))
end

function candidates_main(args)
    if args == ["precompute"]
        precompute_candidates()
    else
        println("Usage: julia examples/density_candidates.jl precompute")
        isempty(args) || error("Invalid arguments")
    end
end
abspath(PROGRAM_FILE) == abspath(@__FILE__) && candidates_main(ARGS)
