# Public API (exported from GIModel.jl): (none — use `GIModel.fn` in tests/scripts)

# Godfrey-Isgur parameterizes α_s(r) as a sum of error functions; use the
# library erf while keeping these thin wrappers to centralize derivative formulas.
function gi_erf(x::Real)
    erf(float(x))
end

function gi_erf_prime(x::Real)
    z = float(x)
    2 / sqrt(π) * exp(-z^2)
end

function gi_erf_second(x::Real)
    z = float(x)
    -4z / sqrt(π) * exp(-z^2)
end

alpha_s_r(r::Real) = sum(a * gi_erf(g * r) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

function central_potential(r::Real, params::GIParameters)
    params.b * r - (4 / 3) * alpha_s_r(r) / r + params.c
end

"""Coulomb piece G(r) = -4 α_s / (3 r) from the text; independent of the linear + constant term S(r) = b r + c."""
function static_coulomb_G(r::Real, params::GIParameters)
    ri = max(float(r), 1.0e-12)
    -(4 / 3) * alpha_s_r(ri) / ri
end

function static_confinement_S(r::Real, params::GIParameters)
    params.b * float(r) + params.c
end
