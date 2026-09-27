@testset "Annihilation amplitudes, interference, and widths" begin
    masses = QuarkMassTable("c"=>1.628)
    w = OscillatorWave(0, 0.5, [1.0])
    state(coefficient=1.0) = PhysicalState("eta_c", 2.98, [(
        basis=BasisState(1,"S",1,0; flavors=(:c,:c)), coefficient=coefficient, wave=w)])
    term = AnnihilationTerm((:c,:c), 4/9)
    op = TwoPhotonAnnihilation(masses, term)
    a = matrix_element(TwoPhotonChannel(), op, state())
    expected = two_photon_amplitude(:P,w,1.628,2.98,4/9)
    @test a.value ≈ expected
    @test decay_width(a) ≈ 1000expected^2
    @test matrix_element(TwoPhotonChannel(),op,state(im)).value ≈ im*a.value
    @test decay_width(TwoPhotonChannel(),op,state(im)) ≈ decay_width(a)
    mixed = PhysicalState("cancel",2.98,[
        (basis=BasisState(1,"S",1,0; flavors=(:c,:c)),coefficient=inv(sqrt(2)),wave=w),
        (basis=BasisState(2,"S",1,0; flavors=(:c,:c)),coefficient=-inv(sqrt(2)),wave=w)])
    @test abs(matrix_element(TwoPhotonChannel(),op,mixed).value) < 1e-14
    current = LeptonicCurrent(:P_P,masses,AnnihilationTerm((:c,:c),2sqrt(3)))
    f = matrix_element(Vacuum(),current,state()).value
    @test f ≈ 2sqrt(3)*leptonic_decay_factor(:P_P,w,1.628,1.628,2.98)
    final = LeptonNeutrinoChannel(0.105; ckm=0.8)
    @test decay_width(final,current,state()) ≈ 1000*0.8^2*leptonic_pseudoscalar_width(abs(f),2.98,0.105)
    g = GluonicAnnihilation(masses,0.3)
    S = wavefunction_origin_smearing(w,1.628)
    # Eq. (17) prefactors, evaluated directly in MeV.
    @test decay_width(TwoGluonChannel(),g,state()) ≈
        1000 * 8π * 0.3^2 / (3 * 1.628^2) * abs2(S) rtol=1e-12
    for (L,mult,J,prefactor,final) in (
            (0,3,1,40*(π^2-9)*0.3^3/(81*1.628^2),ThreeGluonChannel()),
            (1,3,0,8π*0.3^2/(3*1.628^2),TwoGluonChannel()),
            (1,3,2,32π*0.3^2/(45*1.628^2),TwoGluonChannel()))
        wave = OscillatorWave(L,0.5,[1.0])
        st = PhysicalState("test",3.5,[(basis=BasisState(1,L==0 ? "S" : "P",mult,J;
            flavors=(:c,:c)),coefficient=1.0,wave=wave)])
        overlap = wavefunction_origin_smearing(wave,1.628; L)
        @test decay_width(final,g,st) ≈ 1000prefactor*abs2(overlap) rtol=1e-12
        if J == 2
            @test matrix_element(TwoPhotonChannel(),op,st).value ≈
                two_photon_amplitude(:P2,wave,1.628,3.5,4/9)
        elseif L == 0
            vc = LeptonicCurrent(:V_V,masses,AnnihilationTerm((:c,:c),sqrt(16/3)))
            fv = matrix_element(Vacuum(),vc,st).value
            @test decay_width(MasslessLeptonPair(),vc,st) ≈
                1000dilepton_vector_width(abs(fv),3.5)
        end
    end
end

@testset "Derived electromagnetic current coefficients" begin
    # Analytic charge/flavor expectations independent of the implementation.
    for (L, flavor, expected) in ((0,:s,-sqrt(4/3)), (2,:s,-sqrt(8/27)),
                                  (0,:c,sqrt(16/3)), (2,:c,sqrt(32/27)),
                                  (0,:b,-sqrt(4/3)), (2,:b,-sqrt(8/27)))
        @test vector_current_prefactor(BasisState(1,L==0 ? "S" : "D",3,1;
            flavors=(flavor,flavor))) ≈ expected
    end
    masses = QuarkMassTable("u"=>0.22,"d"=>0.22,"c"=>1.628)
    op = LeptonicCurrent(:electromagnetic,masses)
    for (L,sign_d,expected) in ((0,-1,sqrt(6)), (0,1,sqrt(2/3)),
                               (2,-1,sqrt(4/3)), (2,1,sqrt(4/27)))
        wave = OscillatorWave(L,0.5,[1.0])
        st = PhysicalState("neutral light vector",1.5,[
            (basis=BasisState(1,L==0 ? "S" : "D",3,1;flavors=(f,f)),
             coefficient=a/sqrt(2),wave=wave) for (f,a) in ((:u,1),(:d,sign_d))])
        result = matrix_element(Vacuum(),op,st)
        V = leptonic_decay_factor(L==0 ? :V_V : :Vp_V,wave,0.22,0.22,1.5)
        @test result.value ≈ expected*V
        @test decay_width(MasslessLeptonPair(),op,st) ≈
            GEV_TO_MEV*dilepton_vector_width(abs(result.value),1.5)
    end
    mixed = PhysicalState("S-D mixing",3.8,[
        (basis=BasisState(1,L==0 ? "S" : "D",3,1;flavors=(:c,:c)),
         coefficient=a,wave=OscillatorWave(L,0.5,[1.0])) for (L,a) in ((0,0.8),(2,0.6im))])
    result = matrix_element(Vacuum(),op,mixed)
    # Explicit S/D coherent sum with the analytic charge/angular coefficients.
    s_wave = OscillatorWave(0,0.5,[1.0])
    d_wave = OscillatorWave(2,0.5,[1.0])
    expected = 0.8sqrt(16/3)*leptonic_decay_factor(:V_V,s_wave,1.628,1.628,3.8) +
               0.6im*sqrt(32/27)*leptonic_decay_factor(:Vp_V,d_wave,1.628,1.628,3.8)
    @test result.value ≈ expected rtol=1e-12
end
