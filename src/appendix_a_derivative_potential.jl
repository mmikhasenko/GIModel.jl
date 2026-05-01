# Public API (exported from GIModel.jl): (none — use `GIModel.fn` in tests/scripts)

# Appendix-A-oriented finite-difference derivative proxy for the central potential.
#
# The paper's final spin-independent implementation uses expanded effective
# forms (A12)-(A13) in an HO basis. Until the exact coefficients are audited from
# the PDF scan, this path implements the first Gaussian-smearing derivative term
# for the Coulomb block G only:
#
#   exp(∇² / (4σ²)) G ≈ G + ∇²G / (4σ²)
#
# with S(r)=br+c left pointwise, matching the simplification noted in
# docs/appendix_a_from_paper.md. This is intentionally a named comparator, not
# the production default.

function radial_laplacian_values(v::AbstractVector{<:Real}, r::AbstractVector{<:Real})
    n = length(v)
    n == length(r) || throw(ArgumentError("v and r must have the same length"))
    n == 0 && return Float64[]
    n < 3 && return zeros(Float64, n)
    rf = collect(float.(r))
    h = rf[2] - rf[1]
    h > 0 || throw(ArgumentError("r grid must be increasing"))
    for i = 2:(n-1)
        isapprox(rf[i+1] - rf[i], h; rtol = 1e-8, atol = 1e-12) ||
            throw(ArgumentError("radial_laplacian_values expects a uniform grid"))
    end
    vf = collect(float.(v))
    out = zeros(Float64, n)
    for i = 2:(n-1)
        d1 = (vf[i+1] - vf[i-1]) / (2h)
        d2 = (vf[i+1] - 2vf[i] + vf[i-1]) / h^2
        out[i] = d2 + 2d1 / max(rf[i], 1.0e-12)
    end
    out[1] = out[2]
    out[n] = out[n-1]
    return out
end

"""
    appendix_a_derivative_central_values(params, m1, m2, r; order=1)

Finite-difference derivative-expansion proxy for Appendix A: first-order
Gaussian-smearing correction to the Coulomb block G(r), with S(r) kept pointwise.
The default `order=1` corresponds to `G + ∇²G/(4σ²)`.
"""
function appendix_a_derivative_central_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector{<:Real};
    order::Integer = 1,
)
    rf = collect(float.(r))
    isempty(rf) && return Float64[]
    g0 = [static_coulomb_G(ri, params) for ri in rf]
    s0 = [static_confinement_S(ri, params) for ri in rf]
    order <= 0 && return g0 .+ s0
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    lap_g = radial_laplacian_values(g0, rf)
    g_eff = g0 .+ lap_g ./ (4σ^2)
    return g_eff .+ s0
end
