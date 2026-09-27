@testset "Explicit mass corrections at fixed wavefunctions" begin
    masses = QuarkMassTable("c"=>1.628,"q"=>0.22)
    state(L, mult, J, M; coefficients=[1.0], flavors=(:c,:c)) = PhysicalState("test", M, [(
        basis=BasisState(1, ("S","P","D")[L+1], mult, J; flavors),
        coefficient=1.0, wave=OscillatorWave(L,0.5,coefficients))])
    for (kind,L,mult,J,power) in ((:P_P,0,1,0,1), (:V_V,0,3,1,2),
                                (:Vp_V,2,3,1,2), (:Pp_A1,1,3,1,2))
        initial = state(L,mult,J,3.5)
        op = LeptonicCurrent(kind,masses,AnnihilationTerm((:c,:c),1.0))
        a = matrix_element(Vacuum(),op,initial)
        correction = mass_correction_factor(Vacuum(),op,initial; target_mass=4.0)
        @test correction ≈ (3.5/4.0)^power
        @test correction*a.value ≈ matrix_element(Vacuum(),op,state(L,mult,J,4.0)).value
        @test mass_correction_factor(Vacuum(),op,initial; target_mass=3.5) ≈ 1
    end
    mixed = PhysicalState("S-D",3.8,[(basis=c.basis, coefficient=a, wave=c.wave)
        for (L,a) in ((0,0.8),(2,0.6im)) for c in state(L,3,1,3.8).components])
    em = LeptonicCurrent(:electromagnetic,masses)
    @test mass_correction_factor(Vacuum(),em,mixed; target_mass=4.1) ≈ (3.8/4.1)^2

    gg = TwoPhotonAnnihilation(masses,AnnihilationTerm((:c,:c),4/9))
    for (L,mult,J) in ((0,1,0),(1,3,2))
        initial = state(L,mult,J,3.5)
        correction = mass_correction_factor(TwoPhotonChannel(),gg,initial; target_mass=4.0)
        @test correction ≈ (4/3.5)^1.5
        @test correction*matrix_element(TwoPhotonChannel(),gg,initial).value ≈
            matrix_element(TwoPhotonChannel(),gg,state(L,mult,J,4.0)).value
        @test abs2(correction)*decay_width(TwoPhotonChannel(),gg,initial) ≈
            decay_width(TwoPhotonChannel(),gg,state(L,mult,J,4.0))
    end
    running = M -> 0.8/M
    for (L,mult,J,final,power) in ((0,1,0,TwoGluonChannel(),2),
                                 (1,3,0,TwoGluonChannel(),2),
                                 (1,3,2,TwoGluonChannel(),2),
                                 (0,3,1,ThreeGluonChannel(),3))
        initial = state(L,mult,J,3.5)
        for coupling in (0.3,running)
            op = GluonicAnnihilation(masses,coupling)
            correction = mass_correction_factor(final,op,initial; target_mass=4.0)
            expected = coupling isa Real ? 1.0 : (running(4)/running(3.5))^power
            @test correction ≈ expected
            @test correction*decay_width(final,op,initial) ≈
                decay_width(final,op,state(L,mult,J,4.0))
        end
    end

    vector, pseudoscalar, tensor = state(0,3,1,3.5), state(0,1,0,3.0), state(1,3,2,4.0)
    # Use an allowed charged M2 channel. Neutral ccbar tensor-to-pseudoscalar
    # emission cancels between the two constituent currents.
    charged_p = state(0,1,0,3.0; flavors=(:u,:d))
    charged_t = state(1,3,2,4.0; flavors=(:u,:d))
    for (final,initial,power) in ((vector,tensor,1.5),
                                 (charged_p,charged_t,2.5),
                                 (pseudoscalar,vector,0.0))
        op = PhotonEmission(masses)
        a = matrix_element(final,op,initial)
        q = a.momentum_GeV
        @test mass_correction_factor(final,op,initial; target_momentum=q) ≈ 1
        @test mass_correction_factor(final,op,initial; target_momentum=2q) ≈ 2^power
    end
    # Independent low-level recoil calculation: not just a q power.
    excited = state(0,3,1,4.0; coefficients=[0.8,0.6])
    op = PhotonEmission(masses; recoil_order=2,recoil_form_factor=true)
    a = matrix_element(pseudoscalar,op,excited)
    qt = 0.3
    corrected = a.value*mass_correction_factor(pseudoscalar,op,excited; target_momentum=qt)
    expected = m1_recoil_moment(only(pseudoscalar.components).wave,
        only(excited.components).wave, masses["c"],4/3,qt)*photon_recoil_form_factor(qt)
    @test corrected ≈ expected
end
