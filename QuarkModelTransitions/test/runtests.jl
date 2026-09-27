using Test
using LinearAlgebra
using GIModel
using QuarkModelTransitions
# Internal implementation imports used only by the reference/kernel audit.
using QuarkModelTransitions: ALPHA_EM,
    GEV_TO_MEV,
    GLUONIC_CHANNELS,
    HBARC_FM2,
    ReferenceState,
    allowed_partial_waves,
    axial_tau_width,
    calibrate_strong_decay_model,
    decay_momentum,
    dilepton_vector_width,
    e1_angular_coefficient,
    e1_transition_amplitude,
    gluonic_annihilation_amplitude,
    gluonic_annihilation_width,
    leptonic_decay_factor,
    leptonic_pseudoscalar_width,
    m1_radiative_width,
    m1_recoil_moment,
    m1_transition_moment,
    mock_mean_energy,
    mock_meson_mass,
    mock_meson_overlap,
    mock_meson_radial_moment,
    mock_wave_mass,
    neutral_m1_charge,
    partial_wave_projection,
    photon_momentum,
    photon_recoil_form_factor,
    reduced_decay_amplitude,
    spin_flip_photon_amplitude,
    strong_decay_amplitude,
    two_photon_amplitude,
    vector_current_prefactor,
    wavefunction_origin_smearing
const REPOSITORY_ROOT = dirname(dirname(@__DIR__))
const root = REPOSITORY_ROOT

@testset "QuarkModelTransitions" begin
    include("flavor_algebra.jl")
    include("spin_algebra.jl")
    include("algebraic_decomposition.jl")
    include("pseudoscalar_emission.jl")
    include("observables.jl")
    include("annihilation.jl")
    include("mass_correction.jl")
    include("radiative_decays.jl")
    include("strong_decays.jl")
    include("transition_amplitudes.jl")
end

include("../scripts/audit_documentation.jl")
audit_documentation()
