# Paper prescription for isoscalar annihilation mixing: digitized Fig. 5
# pseudoscalar targets and the Table III channel-amplitude table. The numeric
# Eq. (16)-(18) machinery lives in GIModel (pseudoscalar_annihilation.jl and
# the Spectrum bridges in flavor_mixing.jl); this file only encodes which
# channel gets which amplitude and the digitized calibration data.
#
# Public API (exported from GIPaper.jl): GI_PSEUDOSCALAR_FIG5_TARGETS_GEV,
#   table_iii_amplitude

"""Digitized GI Fig. 5 isoscalar pseudoscalar masses (η, η′, η(2S), η′(2S)) in GeV."""
const GI_PSEUDOSCALAR_FIG5_TARGETS_GEV = (0.520, 0.960, 1.440, 1.630)

"""
    table_iii_amplitude(params, L_label, multiplicity, J) -> Union{Float64,Nothing}

Table III general annihilation amplitudes: only `^3S_1` and `^3P_2` carry a
nonzero `A`; the pseudoscalars use P1/P2, and every other channel is ideally
mixed ("Other states ... have been assumed for now to be ideally mixed").
"""
function table_iii_amplitude(
    params::GIParameters,
    L_label::AbstractString,
    multiplicity::Integer,
    J::Integer,
)
    multiplicity == 3 || return nothing
    L_label == "S" && J == 1 && return params.annihilation.s1_A
    L_label == "P" && J == 2 && return params.annihilation.a_3p2
    return nothing
end
