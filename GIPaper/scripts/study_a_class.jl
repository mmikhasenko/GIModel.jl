#!/usr/bin/env julia
# julia --project=GIPaper/scripts GIPaper/scripts/study_a_class.jl [--pilot]
using GIModel, QuarkModelTransitions, CSV, SHA, TOML, LinearAlgebra, Printf
include("a_class_kernels.jl")
using .AClassKernels
BLAS.set_num_threads(1)

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const OUTPUT = joinpath(ROOT,"GIPaper/docs/a_class_universality")
const Q_THRESHOLD = [0.12,0.06,0.03,0.015,0.0075]
const Q_SCAN = [0.15,0.30,0.60]
const BETA = 0.4
const PARAMETER_FILE = joinpath(ROOT,"data/parameters.provisional.toml")
const SOLVERS = (HO=OscillatorSolver(nbasis=48,max_nbasis=128,energy_tolerance_GeV=1e-6,nlevels_per_channel=1),
    FD1200=FiniteDifferenceSolver(ngrid=1200,rmax=48.0,nlevels_per_channel=1),
    FD1800=FiniteDifferenceSolver(ngrid=1800,rmax=48.0,nlevels_per_channel=1))

pairname(p) = join(string.(p),"-")
basisname(t) = "1^$(t.multiplicity)$(t.L)$(t.J)"
daughtername(t) = "1^$(t.daughter_spin)$(t.daughter_L)$(t.daughter_J)"
routeid(t,r) = join((r.sector,pairname(r.parent),pairname(r.daughter),pairname(r.emitted),t.id,
                     basisname(t),daughtername(t),"L$(t.relative_L)",r.topology),"|")
pointid(t,r,q,backend) = routeid(t,r)*"|q=$(q)|$(backend)"

function run_study(; pilot_only=false)
    mkpath(OUTPUT)
    params,mq=load_parameters_and_quark_masses(PARAMETER_FILE)
    parameter_hash=bytes2hex(sha256(read(PARAMETER_FILE)))
    waves=Dict{Any,Any}()
    integral_cache=Dict{Any,Any}()
    sho_waves=Dict(L=>OscillatorWave(L,BETA,[1.0]) for L in 0:2)
    wave_rows=NamedTuple[]
    sample_rows=NamedTuple[]
    function getwave(pair,L,multiplicity,backend)
        # u/d are exactly mass-degenerate; sharing the solve does not merge traces.
        masses=ConstituentMasses(mq[string(pair[1])],mq[string(pair[2])])
        key=(min(masses.m1_GeV,masses.m2_GeV),max(masses.m1_GeV,masses.m2_GeV),L,
             L=="S" ? multiplicity : 0,backend)
        get!(waves,key) do
            solver=getproperty(SOLVERS,backend)
            solution=L=="S" ? contact_hyperfine_nonperturbative_states(params,masses,L,multiplicity,1;solver) :
                channel_solution(params,masses,orbital_angular_momentum(L);solver,nlevels=1)
            isnothing(solution) && error("contact-resummed S waves unavailable")
            wave=radial_wave(solution,1)
            id=join(string.(key),"|")
            push!(wave_rows,(wave_id=id,backend=string(backend),L,multiplicity=key[4],mass1=key[1],mass2=key[2],
                eigenvalue_GeV=solution.eigenvalues_GeV[1],parameter_sha256=parameter_hash,
                prescription=L=="S" ? "central+resummed_contact" : "central",
                solver=repr(solver),convergence=repr(solution.convergence),
                wave_type=string(typeof(wave)),beta=wave isa OscillatorWave ? wave.beta : NaN,
                ho_coefficients=wave isa OscillatorWave ? join(wave.coefficients,";") : "mesh:wave_samples.csv"))
            if !(wave isa OscillatorWave)
                for (i,(r,u)) in enumerate(zip(wave.r,wave.u))
                    push!(sample_rows,(wave_id=id,index=i,r_GeV_inv=r,u=u))
                end
            end
            (wave=wave,mass=solution.eigenvalues_GeV[1],id=id)
        end
    end
    function evaluate(t,r,q,backend)
        native_backend=backend==:SHO ? :HO : backend
        p=getwave(r.parent,t.L,t.multiplicity,native_backend)
        d=getwave(r.daughter,t.daughter_L,t.daughter_spin,native_backend)
        wp=backend==:SHO ? sho_waves[orbital_angular_momentum(t.L)] : p.wave
        wd=backend==:SHO ? sho_waves[orbital_angular_momentum(t.daughter_L)] : d.wave
        c=coefficients(t,r,wp,wd,mq,q;parent_mass=p.mass,daughter_mass=d.mass,beta=BETA,integral_cache)
        return (c=c,parent=p,daughter=d)
    end

    # Gate the nine transitions with a resolved light flavor and both solvers
    # before expanding the census. The gate is on coefficients, not a fitted U.
    pilot_route=first(flavor_routes())
    pilot=NamedTuple[]
    for t in TRANSITIONS
        values=Dict(b=>evaluate(t,pilot_route,0.30,b).c for b in (:SHO,:HO,:FD1200,:FD1800))
        sholimit=threshold_limit(Q_THRESHOLD,[evaluate(t,pilot_route,q,:SHO).c for q in Q_THRESHOLD],t.relative_L)
        shoratio=sholimit.value.h/sholimit.value.g
        abs(shoratio-t.sho_hg)<1e-8 || error("equal-mass SHO A-class failure: $(t.id): $shoratio")
        for b in (:FD1200,:FD1800)
            cmp=compare_couplings(values[b],values[:HO])
            err=maximum(abs(abs(getproperty(values[b],p)/getproperty(values[:HO],p))-1) for p in (:g,:h))
            push!(pilot,(transition=t.id,backend=string(b),coefficient_error=err,
                         vector_distance=cmp.vector_distance,passed=err<0.03))
        end
    end
    CSV.write(joinpath(OUTPUT,"pilot.csv"),pilot)
    all(r.passed for r in pilot) || error("pilot failed; full survey not run (see pilot.csv)")
    println("Pilot passed; maximum coefficient error = ",maximum(r.coefficient_error for r in pilot)); flush(stdout)
    pilot_only && return

    raw=NamedTuple[]
    results=NamedTuple[]
    certificates=NamedTuple[]
    routes=flavor_routes()
    for (ri,r) in enumerate(routes), t in TRANSITIONS
        evaluated=Dict{Tuple{Symbol,Float64},Any}()
        for b in (:SHO,:HO,:FD1200,:FD1800), q in vcat(Q_THRESHOLD,Q_SCAN)
            ev=evaluate(t,r,q,b)
            evaluated[(b,q)]=ev.c
            id=pointid(t,r,q,b)
            for piece in (:g,:h)
                c=getproperty(ev.c,piece)
                push!(raw,(trace_id=id*"|"*string(piece),point_id=id,sector=r.sector,
                    parent_flavors=pairname(r.parent),daughter_flavors=pairname(r.daughter),
                    emitted_flavors=pairname(r.emitted),light_constituent=r.light_constituent,
                    transition=t.id,family=t.family,parent_basis=basisname(t),daughter_basis=daughtername(t),
                    relative_L=t.relative_L,topology=string(r.topology),emitter=string(r.emitter),
                    spectator=string(r.spectator),q_GeV=q,wave_backend=string(b),operator_piece=string(piece),
                    coefficient_re=real(c),coefficient_im=imag(c),reduced_re=real(c/q^t.relative_L),
                    reduced_im=imag(c/q^t.relative_L),parent_label="basis:"*basisname(t),
                    daughter_label="basis:"*daughtername(t),parent_mass_GeV=ev.parent.mass,
                    daughter_mass_GeV=ev.daughter.mass,emitted_mass_placeholder_GeV=0.14,
                    parent_mixing=1.0,daughter_mixing=1.0,emitted_mixing=1.0,
                    parent_wave_id=b==:SHO ? "SHO|$(t.L)|beta=$BETA" : ev.parent.id,
                    daughter_wave_id=b==:SHO ? "SHO|$(t.daughter_L)|beta=$BETA" : ev.daughter.id,
                    beta_comparator_GeV=BETA,parameter_sha256=parameter_hash,
                    kinematics="formal_common_q_not_on_shell",mixing_policy="unmixed_basis_components"))
            end
        end
        first_row = length(raw) - 2*4*length(vcat(Q_THRESHOLD,Q_SCAN)) + 1
        CSV.write(joinpath(OUTPUT,"traces.csv"), raw[first_row:end]; append=first_row>1)
        # Keep the complete raw evaluation before creating any derived reduction.
        # Checkpoint by route so an interrupted long run retains uncollapsed traces.
        limits=Dict(b=>threshold_limit(Q_THRESHOLD,[evaluated[(b,q)] for q in Q_THRESHOLD],t.relative_L)
                    for b in (:SHO,:HO,:FD1200,:FD1800))
        for q in vcat([0.0],Q_SCAN)
            val(b)=q==0 ? limits[b].value : evaluated[(b,q)]
            sho=val(:SHO)
            ho=val(:HO)
            fd=val(:FD1800)
            fdcoarse=val(:FD1200)
            relative(a,b)=maximum(abs(abs(getproperty(a,p)/getproperty(b,p))-1) for p in (:g,:h))
            solver_error=relative(fd,ho)
            refinement_error=relative(fd,fdcoarse)
            threshold_error=q==0 ? maximum(limits[b].error for b in keys(limits)) : 0.0
            certified=solver_error<0.03 && refinement_error<0.003 && threshold_error<0.002
            cid=pointid(t,r,q,:certificate)
            push!(certificates,(certificate_id=cid,route_id=routeid(t,r),q_GeV=q,
                ho_fd_error=solver_error,fd_refinement_error=refinement_error,
                threshold_error,certified))
            for b in (:HO,:FD1200,:FD1800)
                c=val(b)
                cmp=compare_couplings(c,sho)
                me,ms=mq[string(r.emitter)],mq[string(r.spectator)]
                source_ids(backend,p)=join([pointid(t,r,x,backend)*"|$p" for x in (q==0 ? Q_THRESHOLD : [q])],";")
                push!(results,(point_id=pointid(t,r,q,b),route_id=routeid(t,r),sector=r.sector,
                    parent_flavors=pairname(r.parent),daughter_flavors=pairname(r.daughter),
                    emitted_flavors=pairname(r.emitted),light_constituent=r.light_constituent,
                    transition=t.id,family=t.family,relative_L=t.relative_L,topology=string(r.topology),
                    q_GeV=q,regime=q==0 ? "threshold_extrapolation" : "common_q",backend=string(b),
                    Cg_re=real(c.g),Cg_im=imag(c.g),Ch_re=real(c.h),Ch_im=imag(c.h),
                    SHO_Cg_re=real(sho.g),SHO_Cg_im=imag(sho.g),SHO_Ch_re=real(sho.h),SHO_Ch_im=imag(sho.h),
                    SHO_hg_re=real(sho.h/sho.g),SHO_hg_im=imag(sho.h/sho.g),
                    mass_formula_applicable=t.family=="A",
                    requested_emitter_hg=t.family=="A" ? me/(2*(me+ms)) : NaN,
                    implemented_spectator_hg=t.family=="A" ? ms/(2*(me+ms)) : NaN,
                    equal_mass_SU3_hg=t.sho_hg,
                    requested_formula_pass=t.family=="A" ? abs(sho.h/sho.g-me/(2*(me+ms)))<1e-8 : missing,
                    Rg_re=real(cmp.Rg),Rg_im=imag(cmp.Rg),Rh_re=real(cmp.Rh),Rh_im=imag(cmp.Rh),
                    U_hg_re=real(cmp.U_hg),U_hg_im=imag(cmp.U_hg),
                    scalar_stable=cmp.scalar_stable,vector_valid=cmp.vector_valid,vector_distance=cmp.vector_distance,
                    certified,certificate_id=cid,
                    native_g_trace_ids=source_ids(b,:g),native_h_trace_ids=source_ids(b,:h),
                    sho_g_trace_ids=source_ids(:SHO,:g),sho_h_trace_ids=source_ids(:SHO,:h)))
            end
        end
        if t.id==last(TRANSITIONS).id
            println("Completed flavor route $ri/$(length(routes))"); flush(stdout)
        end
    end
    CSV.write(joinpath(OUTPUT,"waves.csv"),wave_rows)
    CSV.write(joinpath(OUTPUT,"wave_samples.csv"),sample_rows)
    CSV.write(joinpath(OUTPUT,"comparisons.csv"),results)
    CSV.write(joinpath(OUTPUT,"certificates.csv"),certificates)
    @assert length(unique(r.trace_id for r in raw))==length(raw)
    open(joinpath(OUTPUT,"run.toml"),"w") do io
        TOML.print(io,Dict("schema_version"=>2,"julia_version"=>string(VERSION),
            "source_revision"=>strip(read(`git -C $ROOT rev-parse HEAD`,String)),
            "parameter_sha256"=>parameter_hash,"beta_GeV"=>BETA,"threshold_q_GeV"=>Q_THRESHOLD,
            "common_q_GeV"=>Q_SCAN,"trace_count"=>length(raw),"flavor_route_count"=>length(routes),
            "transition_count"=>length(TRANSITIONS),"numerical_certificate_tolerance"=>0.03,"fd_refinement_tolerance"=>0.003,
            "threshold_tolerance"=>0.002,"scope"=>"unmixed n=1 spectroscopic kernels; no physical width prediction",
            "operator_convention"=>"existing Eq19 spectator coordinate fraction; emitter formula discrepancy retained"))
    end
    println("Wrote $(length(raw)) traces and $(length(results)) comparisons to $OUTPUT")
end

abspath(PROGRAM_FILE)==(@__FILE__) && run_study(pilot_only="--pilot" in ARGS)
