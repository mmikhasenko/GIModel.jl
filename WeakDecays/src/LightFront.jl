# Included inside WeakDecayExploration. Only the zero-mode-free P->P plus
# current at q+=0 is implemented, at spacelike q2 and q2=0. This is NOT a full
# covariant light-front implementation or a timelike analytic continuation.
export lightfront_fp, lf_kinematics
function lf_kinematics(x,kt2,m1,m2)
    M0=sqrt((m1*m1+kt2)/x+(m2*m2+kt2)/(1-x))
    kz=(x-.5)*M0+(m2*m2-m1*m1)/(2M0)
    jac=M0/(4x*(1-x))*(1-((m1*m1-m2*m2)/M0^2)^2)
    (;k=sqrt(kt2+kz*kz),jac,A=(1-x)*m1+x*m2)
end
"""SLF P->P good-current overlap with Melosh spin factor; Chang et al.
1810.00296, Sec. 3.3. Arbitrary normalized GI radial functions replace a Gaussian.
Measure dx d²k/(4pi), with sqrt(dkz/dx) for each radial wave. At Q=0 the
azimuth is analytic. No external meson mass enters this particular overlap.
"""
function lightfront_fp(c::CurrentModel,Q2=0.;nx=100,nk=100,nphi=32,kmax=10.,spin=true)
    Q2>=0 || throw(ArgumentError("Q2=-q2 must be nonnegative; timelike LF not implemented"))
    xs,wx=gauss(nx,0.,1.); ks,wk=gauss(nk,0.,kmax)
    phis,wp=Q2==0 ? ([0.],[2pi]) : gauss(nphi,0.,2pi)
    Q=sqrt(Q2); result=0.
    for (x,ax) in zip(xs,wx),(k,ak) in zip(ks,wk)
        ki=lf_kinematics(x,k*k,c.mi,c.ms)
        for (phi,ap) in zip(phis,wp)
            kfx=k*cos(phi)-(1-x)*Q; kfy=k*sin(phi)
            ktf2=kfx*kfx+kfy*kfy; kdot=k*k-(1-x)*Q*k*cos(phi)
            kf=lf_kinematics(x,ktf2,c.mf,c.ms)
            melosh=spin ? (ki.A*kf.A+kdot)/sqrt((ki.A^2+k*k)*(kf.A^2+ktf2)) : 1.
            result+=ax*ak*ap*k/(4pi)*sqrt(ki.jac*kf.jac)*c.wi(ki.k)*c.wf(kf.k)*melosh
        end
    end
    result
end
