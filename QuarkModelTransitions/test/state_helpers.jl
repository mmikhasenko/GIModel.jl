@testset "Coherent superpositions and unsupported photon kernels" begin
    wave = OscillatorWave(0, 0.5, [1.0])
    up = PhysicalState("u", 0.8, [(basis=BasisState(1,"S",3,1; flavors=(:u,:u)), coefficient=1.0, wave)])
    down = PhysicalState("d", 0.8, [(basis=BasisState(1,"S",3,1; flavors=(:d,:d)), coefficient=1.0, wave)])
    omega = superpose([up, down], [1, 1]; mass_GeV=0.8)
    rho = superpose([up, down], [1, -1]; mass_GeV=0.8)
    _, mq = load_parameters_and_quark_masses(default_parameters_path())
    em = LeptonicCurrent(:electromagnetic, mq)
    @test matrix_element(Vacuum(), em, rho).value ≈ 3matrix_element(Vacuum(), em, omega).value
    repeated = superpose([up, up], [1, 1]; mass_GeV=0.8)
    @test matrix_element(Vacuum(), em, repeated).value ≈ matrix_element(Vacuum(), em, up).value
    @test_throws ArgumentError superpose([up, up], [1, -1]; mass_GeV=0.8)
    @test_throws ArgumentError superpose([up], [1, 2]; mass_GeV=0.8)
    mixed = PhysicalState("vector", 3.1, [
        (basis=BasisState(1,"S",3,1; flavors=(:c,:c)), coefficient=sqrt(0.99), wave),
        (basis=BasisState(1,"D",3,1; flavors=(:c,:c)), coefficient=0.1, wave=OscillatorWave(2,0.5,[1.0])),
    ])
    initial = PhysicalState("P", 3.5, [
        (basis=BasisState(1,"P",3,1; flavors=(:c,:c)), coefficient=1.0, wave=OscillatorWave(1,0.5,[1.0])),
    ])
    exception = try
        matrix_element(mixed, PhotonEmission(mq), initial)
        nothing
    catch err
        err
    end
    @test exception isa ErrorException
    @test exception.msg == "Not implemented yet. Please submit issue if needed, and/or PR with implementation."
end
