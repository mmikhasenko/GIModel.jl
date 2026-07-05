# Public API (exported from GIModel.jl): (none — use `GIModel.fn` in tests/scripts)

#
# 3D isotropic Gaussian smearing of a spherically symmetric radial function V (|r|),
# Appendix A, Eqs. (A7)–(A8), PDF p. 36. Same σ as contact_smearing_sigma (A9), Table II.
# R > 0: one-dimensional form from the angle-integrated convolution (e.g. difference of
# Gaussians). R → 0: direct radial integral with isotropic 3D Gaussian.
#
function smear_3d_radial(v::AbstractVector{<:Real}, r::AbstractVector{<:Real}, σ::Real)
    n = length(r)
    n == 0 && return eltype(r)[]
    n == 1 && return v
    h = r[2] - r[1]
    w = fill(h, n)
    w[1] = h / 2
    w[n] = h / 2
    σf = max(float(σ), 1.0e-12)
    R0 = 0.25 * h
    out = similar(r, Float64)
    for i in eachindex(r)
        R = r[i]
        s = 0.0
        if R < R0
            for j in eachindex(r)
                rp = r[j]
                # ρ(r) = σ^3 / π^(3/2) exp(-σ^2 r^2), with σ in GeV and r in GeV^-1.
                ρ = σf^3 / (π^(3 / 2)) * exp(-(σf * rp)^2)
                s += 4 * π * rp^2 * w[j] * ρ * v[j]
            end
        else
            for j in eachindex(r)
                rp = r[j]
                # Angle-integrated convolution of a 3D isotropic Gaussian:
                # f̃(R) = (σ / (√π R)) ∫ dr' r' [e^{-σ^2 (R-r')^2} - e^{-σ^2 (R+r')^2}] f(r')
                pre = σf / (sqrt(π) * R)
                s +=
                    w[j] *
                    pre *
                    rp *
                    (exp(-(σf * (R - rp))^2) - exp(-(σf * (R + rp))^2)) *
                    v[j]
            end
        end
        out[i] = s
    end
    return out
end

# Experimental. Not the GI (A12)–(A13) smeared potential used in the paper: (A7)–(A8)
# with the Table II width applied to the pointwise G and S (Eqs. (11)–(13) orient.)
# is numerically uncontrolled on a fixed radial line when $\sigma$ is O(1): the 3D
# convolution weights the large-$r$ region by volume and can remove the $1/r$ well.
# Enable only for research; production defaults keep `appendix_a_smearing = false`.
function smeared_central_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector{<:Real},
)
    σ = contact_smearing_sigma(params, m1, m2)
    n = length(r)
    n < 2 && return [central_potential(ri, params) for ri in r]
    h = r[2] - r[1]
    rmax0 = r[end]
    # Kernel tail: exp(-(σ Δr)^2) at Δr = 8/σ gives exp(-64), effectively zero.
    n_tail = σ > 0 ? max(0, Int(ceil(8 / (σ * h)))) : 0
    r_ext =
        n_tail > 0 ? vcat(r, collect(range(rmax0 + h, rmax0 + n_tail * h; step = h))) : r
    g0 = [static_coulomb_G(ri, params) for ri in r_ext]
    s0 = [static_confinement_S(ri, params) for ri in r_ext]
    vsum = smear_3d_radial(g0, r_ext, σ) .+ smear_3d_radial(s0, r_ext, σ)
    return vsum[1:n]
end

function smeared_coulomb_G_closed(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = float(r)
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    τs = map(γ -> 1 / sqrt(1 / σ^2 + 1 / γ^2), ALPHA_GAMMAS)
    if abs(ri) < 1.0e-8
        return -sum(8 * α * τ / (3 * sqrt(π)) for (α, τ) in zip(ALPHA_COEFFS, τs))
    end
    -sum(4 * α * erf(τ * ri) / (3 * ri) for (α, τ) in zip(ALPHA_COEFFS, τs))
end

function smeared_confinement_S_closed(params::GIParameters, m1::Real, m2::Real, r::Real)
    ri = float(r)
    σ = max(contact_smearing_sigma(params, m1, m2), 1.0e-12)
    if abs(ri) < 1.0e-8
        return 2 * params.potential.b / (sqrt(π) * σ) + params.potential.c
    end
    z = σ * ri
    bracket = exp(-z^2) / (sqrt(π) * z) + (1 + 1 / (2 * z^2)) * erf(z)
    params.potential.b * ri * bracket + params.potential.c
end

function appendix_a_closed_central_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector{<:Real},
)
    [
        smeared_coulomb_G_closed(params, m1, m2, ri) +
        smeared_confinement_S_closed(params, m1, m2, ri) for ri in r
    ]
end
