# 1D radial Gaussian renormalization of the pointwise Coulomb G(r) on a uniform grid.
# This is a *reduced* smearing: same σ family as the contact (A9), no 3D volume factors.
# Not a transcription of the full (A12) derivative expansion (paper uses HO basis there).

function convolve_1d_gaussian_same_length(
    v::Vector{Float64},
    r::Vector{Float64},
    _h::Real,
    σ::Real,
)
    n = length(v)
    n == 0 && return v
    σf = max(float(σ), 1.0e-9)
    out = similar(v)
    for i in 1:n
        ri = r[i]
        s = 0.0
        wsum = 0.0
        for j in 1:n
            w = exp(-0.5 * ((ri - r[j]) / σf)^2)
            s += w * v[j]
            wsum += w
        end
        out[i] = wsum > 0.0 ? s / wsum : v[i]
    end
    return out
end

"""Pointwise S(r) + 1D Gaussian-smeared Coulomb G(r) (same Table II σ as contact (A9))."""
function coulomb_1d_smeared_central_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector{<:Real},
)
    r = collect(float.(r))
    h = r[2] - r[1]
    σ = contact_smearing_sigma(params, m1, m2)
    g0 = [static_coulomb_G(ri, params) for ri in r]
    g1 = convolve_1d_gaussian_same_length(g0, r, h, σ)
    s0 = [static_confinement_S(ri, params) for ri in r]
    return g1 .+ s0
end
