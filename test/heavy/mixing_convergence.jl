@testset "Unequal-mass P-wave mixing agrees across radial solvers" begin
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    levels = spectrum_levels(1; L_labels = ("P",))
    angle(spec) = begin
        components = spectrum_state(spec, "1^1P_1").mixings[1].components
        rad2deg(atan(components[2], components[1]))
    end
    for flavors in ((:u, :s), (:c, :u), (:b, :u))
        meson = Meson(mq, flavors...)
        fd = compute_spectrum(params, meson; levels,
                              solver = FiniteDifferenceSolver(ngrid = 900))
        ho = compute_spectrum(params, meson; levels, solver = OscillatorSolver())
        common = compute_spectrum(params, meson; levels,
                                  solver = OscillatorSolver(beta_grid = [0.6]))
        # Energy and wave truncation still limit the solved states. The former
        # 0.7--2 degree discrepancy was an independent operator truncation.
        @test abs(angle(fd) - angle(ho)) < 0.1
        @test abs(angle(common) - angle(ho)) < 0.1
    end
end
