module QuarkModelTransitions

using ..GIModel
using PartialWaveFunctions: CG
using Printf

# Supported API. Undeclared implementation names are internal.
export PhysicalState,
    physical_state,
    TwoMesonChannel,
    PartialWave,
    matrix_element,
    partial_waves,
    mass_correction_factor,
    LeptonicCurrent,
    TwoPhotonAnnihilation,
    GluonicAnnihilation,
    Vacuum,
    TwoPhotonChannel,
    TwoGluonChannel,
    ThreeGluonChannel,
    LeptonNeutrinoChannel,
    MasslessLeptonPair,
    AnnihilationTerm,
    PseudoscalarEmission,
    decay_width,
    charge_radius_squared,
    PhotonEmission

include("observable_inputs.jl")

include("flavor_algebra.jl")
include("spin_algebra.jl")
include("algebraic_decomposition.jl")

include("transition_amplitudes.jl")

include("mass_correction.jl")
include("annihilation_channels.jl")
include("annihilation_term.jl")
include("annihilation_amplitude.jl")

include("pseudoscalar_emission.jl")

include("strong_decays.jl")

# =============================================================================
# Table VII annihilation widths: gluonic QQ̄ → gluons
# =============================================================================

include("annihilation_widths.jl")
include("leptonic_current.jl")
include("two_photon_annihilation.jl")
include("gluonic_annihilation.jl")

include("mock_meson_overlaps.jl")

include("radiative_decays.jl")

include("photon_emission.jl")

end
