# Public API (exported from GIModel.jl): (none — use `GIModel.fn` in tests/scripts)

# Godfrey–Isgur parameterizes α_s(Q²) as a sum of Gaussians. The coordinate-space
# erf profile is the corresponding Coulomb kernel representation.
function erf_prime(x::AbstractFloat)
    2 / sqrt(π) * exp(-x^2)
end

function erf_second(x::AbstractFloat)
    -4x / sqrt(π) * exp(-x^2)
end

alpha_s_r(r::Real) = sum(a * erf(g * r) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

alpha_s_q(Q::Real) = sum(a * exp(-(Q^2) / (4g^2)) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

function central_potential(r::Real, params::GIParameters)
    params.b * r - (4 / 3) * alpha_s_r(r) / r + params.c
end

"""Coulomb piece G(r) = -4 α_s / (3 r) from the text; independent of the linear + constant term S(r) = b r + c."""
function static_coulomb_G(r::Real, params::GIParameters)
    ri = max(r, 1.0e-12)
    -(4 / 3) * alpha_s_r(ri) / ri
end

function static_confinement_S(r::Real, params::GIParameters)
    params.b * r + params.c
end
