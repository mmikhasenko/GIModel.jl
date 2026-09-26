# Public API (exported from GIModel.jl):
#   alpha_s_q, alpha_s_r

# Godfrey–Isgur parameterizes α_s(Q²) as a sum of Gaussians. The coordinate-space
# erf profile is the corresponding Coulomb kernel representation.
function erf_prime(x::AbstractFloat)
    2 / sqrt(π) * exp(-x^2)
end

function erf_second(x::AbstractFloat)
    -4x / sqrt(π) * exp(-x^2)
end

"""
    alpha_s_r(r) -> Float64

The coordinate-space running coupling ``\\alpha_s(r) = \\sum_k \\alpha_k\\,
\\mathrm{erf}(\\gamma_k r)``, the Coulomb-kernel representation of
[`alpha_s_q`](@ref). Used by the static potential ``G(r) = -4\\alpha_s(r)/3r``.
"""
alpha_s_r(r::Real) = sum(a * erf(g * r) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

"""
    alpha_s_q(Q) -> Float64

The Godfrey–Isgur running coupling ``\\alpha_s(Q^2) = \\sum_k \\alpha_k\\,
e^{-Q^2/4\\gamma_k^2}`` (Eq. 12), `Q` in GeV — three fixed Gaussians with **no
quark-mass or flavor dependence at all**: there is no ``n_f`` threshold, so
the historical `Lambda_MeV` reference value is not a runtime input. This is the precise
sense in which the model's medium is flavor-blind; quark mass enters the model
only through [`contact_smearing_sigma`](@ref) and the relativistic weight.
"""
alpha_s_q(Q::Real) = sum(a * exp(-(Q^2) / (4g^2)) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

function central_potential(r::Real, params::GIParameters)
    params.potential.b * r - (4 / 3) * alpha_s_r(r) / r + params.potential.c
end

"""Coulomb piece G(r) = -4 α_s / (3 r) from the text; independent of the linear + constant term S(r) = b r + c."""
function static_coulomb_G(r::Real, params::GIParameters)
    ri = max(r, 1.0e-12)
    -(4 / 3) * alpha_s_r(ri) / ri
end

function static_confinement_S(r::Real, params::GIParameters)
    params.potential.b * r + params.potential.c
end
