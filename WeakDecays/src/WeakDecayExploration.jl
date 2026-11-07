module WeakDecayExploration

using GIModel, GIModel.QuarkModelTransitions, LinearAlgebra
using QuadGK: quadgk
export CurrentModel, current_components, form_factors, differential_rate, integrated_rate,
    leptonic_constant, gauss, spectral_density, folded_rate, factorized_pp_width,
    recoil, GF, HBAR, WaveTable, prepare_wave

const GF = 1.1663787e-5 # GeV^-2
const HBAR = 6.582119569e-25 # GeV s

"""Gauss-Legendre nodes and weights on [a,b], independent of external quadrature packages."""
function gauss(n::Int, a::Real, b::Real)
    n >= 2 || throw(ArgumentError("at least two nodes required"))
    f = eigen(SymTridiagonal(zeros(n), [k/sqrt(4k*k-1) for k in 1:n-1]))
    return (a+b)/2 .+ (b-a)/2 .* f.values, (b-a).*f.vectors[1,:].^2
end

"""Explicit numerical tabulation for research integrals; original solved waves remain unchanged."""
struct WaveTable
    p::Vector{Float64}
    phi::Vector{Float64}
end
function (w::WaveTable)(p::Real)
    p < 0 && throw(ArgumentError("negative radial momentum"))
    p >= last(w.p) && return 0.0
    i = clamp(searchsortedlast(w.p, p), 1, length(w.p)-1)
    t = (p-w.p[i])/(w.p[i+1]-w.p[i])
    (1-t)*w.phi[i] + t*w.phi[i+1]
end
function prepare_wave(w::RadialWave; pmax=12.0, n=6001)
    if w isa GIModel.MeshWave
        mw = momentum_wave(w, 0; pmax, npoints=n)
        return WaveTable(mw.p, mw.phi)
    end
    mw = momentum_wave(w, 0)
    p = collect(range(0,pmax,length=n))
    WaveTable(p, mw.(p))
end

"""
Exploratory S-wave, active-quark spectator current. No Wilson matching, pair currents,
Wigner rotations or claim of Lorentz covariance. :galilean uses k'=k+ms/(mf+ms)P;
:boosted_spectator boosts the on-shell spectator and includes its scalar Jacobian,
but retains canonical unrotated spin states. It is a model-sensitivity diagnostic.
All dimensional inputs/outputs are GeV; widths are GeV (unlike the QMT public API).
"""
struct CurrentModel
    wi::WaveTable
    wf::WaveTable
    mi::Float64
    mf::Float64
    ms::Float64
    scheme::Symbol
    nr::Int
    nz::Int
    pmax::Float64
    function CurrentModel(wi,wf,mi,mf,ms; scheme=:galilean,nr=100,nz=40,pmax=10.)
        scheme in (:galilean,:boosted_spectator) || throw(ArgumentError("unknown recoil scheme"))
        all(x->isfinite(x)&&x>0,(mi,mf,ms,pmax)) || throw(ArgumentError("positive finite inputs required"))
        nr>=2 && nz>=2 || throw(ArgumentError("invalid quadrature"))
        new(wi,wf,mi,mf,ms,scheme,nr,nz,pmax)
    end
end
export CoherentCurrent, current_from_states
struct CoherentCurrent
    terms::Vector{Tuple{ComplexF64,CurrentModel}}
end
"""S-wave physical-state adapter. Changes the first (quark) flavor only;
antiquark-current conventions and non-S-wave components are deliberately rejected.
State mixing coefficients interfere coherently, and CKM remains outside the current.
"""
function current_from_states(final::PhysicalState,initial::PhysicalState,masses;
                             active::Pair{Symbol,Symbol},kwargs...)
    initial.J==0 && initial.parity==-1 || throw(ArgumentError("initial must be pseudoscalar"))
    final.J in (0,1) && final.parity==-1 || throw(ArgumentError("final must be P or V"))
    all(c.basis.L_label=="S" for st in (initial,final) for c in st.components) ||
        throw(ArgumentError("S-wave adapter cannot discard orbital components"))
    terms=Tuple{ComplexF64,CurrentModel}[]
    for ci in initial.components, cf in final.components
        fi,ff=ci.basis.flavors,cf.basis.flavors
        isnothing(fi) || isnothing(ff) ? throw(ArgumentError("explicit flavors required")) : nothing
        fi[1]==first(active) && ff[1]==last(active) && fi[2]==ff[2] || continue
        m=CurrentModel(prepare_wave(ci.wave),prepare_wave(cf.wave),
            masses[string(fi[1])],masses[string(ff[1])],masses[string(fi[2])];kwargs...)
        push!(terms,(ComplexF64(conj(cf.coefficient)*ci.coefficient),m))
    end
    isempty(terms) && throw(ArgumentError("no matching active-quark/spectator route"))
    CoherentCurrent(terms)
end
function current_components(c::CoherentCurrent,M,m,q2)
    hs=[(a,current_components(model,M,m,q2)) for (a,model) in c.terms]
    vals=[sum(a*getproperty(h,key) for (a,h) in hs) for key in (:V0,:Vz,:Ax,:Vy,:A0,:Az)]
    h=first(hs)[2]
    (;V0=vals[1],Vz=vals[2],Ax=vals[3],Vy=vals[4],A0=vals[5],Az=vals[6],P=h.P,E=h.E)
end

recoil(M,m,q2) = sqrt(max(0.,((M+m)^2-q2)*((M-m)^2-q2)))/(2M)

"""Return integrated V0,Vz, transverse Ax and |Im Vy|, longitudinal A0,Az.
Spin coefficients follow singlet epsilon/sqrt(2) and triplet sigma_a epsilon/sqrt(2).
Angular average is dcos(theta)/2, since Phi is a unit-normalized radial wave.
"""
function current_components(c::CurrentModel, M,m,q2)
    M>m>0 || throw(ArgumentError("require M > m > 0"))
    0<=q2<=(M-m)^2 || throw(ArgumentError("q2 outside physical range"))
    P = recoil(M,m,q2); E=sqrt(m*m+P*P)
    ps,wp=gauss(c.nr,0,c.pmax); zs,wz=gauss(c.nz,-1.,1.)
    out=zeros(6)
    for (k,wk) in zip(ps,wp), (z,w) in zip(zs,wz)
        kz=k*z; kt2=k*k*(1-z*z); pz=kz+P
        Ei=sqrt(c.mi^2+k*k); Ef=sqrt(c.mf^2+kt2+pz*pz)
        Di=Ei+c.mi; Df=Ef+c.mf; den=Di*Df
        N=sqrt(Di*Df/(4Ei*Ef))
        if c.scheme==:galilean
            kpz=kz+c.ms/(c.mf+c.ms)*P; jac=1.
        else
            Es=sqrt(c.ms^2+k*k); wboost=E/m
            kpz=wboost*kz+(P/m)*Es
            jac=(wboost*Es+(P/m)*kz)/Es
        end
        kf=sqrt(kt2+kpz*kpz)
        weight=wk*w/2*k*k*c.wi(k)*c.wf(kf)*sqrt(jac)*N
        dotp=k*k+kz*P
        out[1]+=weight*(1+dotp/den)
        out[2]+=weight*(kz/Di+pz/Df)
        out[3]+=weight*(1-(kz*kz+kz*P)/den)
        out[4]+=weight*(pz/Df-kz/Di)
        out[5]+=weight*(kz/Di+pz/Df)
        out[6]+=weight*(1+(kz*kz+kz*P-kt2)/den)
    end
    out .*= 2sqrt(M*E)
    (;V0=out[1],Vz=out[2],Ax=out[3],Vy=out[4],A0=out[5],Az=out[6],P,E)
end

function form_factors(c,M,m,q2; vector=false)
    # At exact zero recoil f+ and A2 cannot be separated from one current value.
    q2 < (M-m)^2 || throw(ArgumentError("use a one-sided limit at zero recoil"))
    h=current_components(c,M,m,q2); q0=M-h.E
    if !vector
        fp=(h.V0+q0*h.Vz/h.P)/(2M)
        fm=fp-h.Vz/h.P
        return (;fp, f0=fp+q2/(M*M-m*m)*fm, fm)
    end
    A1=h.Ax/(M+m)
    A0=(q0*h.A0+h.P*h.Az)/(2M*h.P)
    A2=((M+m)*A1*(h.P^2+q0*h.E)/m-(h.P*h.A0+q0*h.Az))*m*(M+m)/(2M*M*h.P^2)
    V=(M+m)*h.Vy/(2M*h.P)
    (;V,A0,A1,A2)
end

"""Finite-lepton-mass dGamma/dq2; no radiative corrections. Return GeV^-1."""
function differential_rate(c,M,m,q2,ml,ckm; vector=false)
    ml>=0 && isfinite(ml) && 0<=abs(ckm)<=1 || throw(ArgumentError("invalid lepton/CKM input"))
    M>m>0 || throw(ArgumentError("invalid masses"))
    (q2<ml^2 || (ml>0 && q2==ml^2) || q2>=(M-m)^2) && return (total=0.,longitudinal=0.,scalar=0.)
    h=current_components(c,M,m,q2); q0=M-h.E
    # sqrt(q2)*H0 and sqrt(q2)*Ht avoid removable 1/q2 singularities.
    if vector
        qH0=abs2(h.P*h.A0+q0*h.Az)
        qHt=abs2(q0*h.A0+h.P*h.Az)
        qHT=q2*(abs2(h.Ax-h.Vy)+abs2(h.Ax+h.Vy))
    else
        qH0=abs2(h.P*h.V0+q0*h.Vz)
        qHt=abs2(q0*h.V0+h.P*h.Vz)
        qHT=0.
    end
    delta=ml==0 ? 0. : ml^2/q2
    pref=GF^2*abs2(ckm)*h.P/(96pi^3*M*M)*(1-delta)^2
    scalar=pref*3delta/2*qHt
    longitudinal=pref*(1+delta/2)*qH0+scalar
    total=longitudinal+pref*(1+delta/2)*qHT
    (;total,longitudinal,scalar)
end
function integrated_rate(c,M,m,ml,ckm; vector=false,nq=48)
    M>m+ml || return (width=0.,FL=0.,scalar_fraction=0.)
    qs,ws=gauss(nq,ml^2,(M-m)^2)
    rates=[differential_rate(c,M,m,q,ml,ckm;vector) for q in qs]
    width=sum(w*r.total for (w,r) in zip(ws,rates))
    width==0 && return (width=0.,FL=0.,scalar_fraction=0.)
    (;width, FL=sum(w*r.longitudinal for (w,r) in zip(ws,rates))/width,
      scalar_fraction=sum(w*r.scalar for (w,r) in zip(ws,rates))/width)
end

"""Dimensionful fP in GeV from the pre-existing GI Appendix-D annihilation current."""
function leptonic_constant(state,masses)
    terms=[AnnihilationTerm(c.basis.flavors,2sqrt(3)) for c in state.components]
    # This benchmark takes a pure flavor sector; duplicate terms must not be summed.
    length(terms)==1 || throw(ArgumentError("pure single-component benchmark only"))
    real(state.mass_GeV*matrix_element(Vacuum(),LeptonicCurrent(:P_P,masses,terms),state).value)
end

"""P-wave running width spectral function, without finite-window renormalization.
Single elastic channel and phenomenological Blatt-Weisskopf radius R in GeV^-1.
"""
function spectral_density(s,mR,gamma,ma,mb; R=3.,running=true)
    gamma>0 && mR>ma+mb && ma>=0 && mb>=0 && R>=0 || throw(ArgumentError("invalid resonance"))
    s>(ma+mb)^2 || return 0.
    k=recoil(sqrt(s),ma,mb^2); k0=recoil(mR,ma,mb^2)
    g=running ? gamma*(k/k0)^3*mR/sqrt(s)*(1+(R*k0)^2)/(1+(R*k)^2) : gamma
    mR*g/(pi*((s-mR*mR)^2+(mR*g)^2))
end

"""Fold off-shell vector production over s. Fixed radial wave, phenomenological
mass continuation, all elastic charge channels; multiply a charge sub-BF separately.
Tangent map resolves arbitrarily narrow widths and supplies an independent NWA check.
"""
function folded_rate(c,M,mR,gamma,ma,mb,ml,ckm; ns=48,nq=32,R=3.,running=true)
    lo=(ma+mb)^2; hi=(M-ml)^2; hi>lo || return 0.
    scale=mR*gamma
    xs,ws=gauss(ns,atan((lo-mR*mR)/scale),atan((hi-mR*mR)/scale))
    sum(zip(xs,ws)) do (x,w)
        s=mR*mR+scale*tan(x)
        w*scale/cos(x)^2*spectral_density(s,mR,gamma,ma,mb;R,running)*
            integrated_rate(c,M,sqrt(s),ml,ckm;vector=true,nq).width
    end
end

function factorized_pp_width(c,M,m,mh,fh,ckm_product; a1=1.)
    M>m+mh || return 0.
    f0=form_factors(c,M,m,mh^2).f0
    amp=GF/sqrt(2)*ckm_product*a1*fh*(M*M-m*m)*f0
    recoil(M,m,mh^2)/(8pi*M*M)*abs2(amp)
end
include("LightFront.jl")

end
