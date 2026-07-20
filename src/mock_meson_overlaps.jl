# Appendix-D mock-meson transition matrix elements (page 40), promoted from the
# Table VI photon-decay audit into src with tests.
#
# A mock meson is the weakly-bound |q q̄⟩ state whose momentum-space radial wave
# Φ_L(p) is the spherical-Bessel transform of the reduced radial wave u(r), and
# whose mock mass is M̃ = ⟨E₁⟩ + ⟨E₂⟩. The two kernels are:
#
#   I_i(x,y)   = √(4 M̃ₓ M̃_y)/(M̃ₓ + M̃_y) · ∫ dp p² Φₓ(p) Φ_y(p) (1/mᵢ)(mᵢ/Eᵢ)^f
#   Eₙⁱ(x,y)   = |mᵢ / √(⟨Eᵢ⟩ₓ ⟨Eᵢ⟩_y)|^g · ∫ dr uₓ(r) u_y(r) rⁿ
#
# with Eᵢ = √(mᵢ² + p²) the relativistic quark energy and f, g the paper's fitted
# exponents (0.7 for I_i, 0.5 for Eₙⁱ). I_i drives the M1 magnetic moments; Eₙⁱ
# the E1/M2 multipole amplitudes and the hindered-transition recoil term.

"""A mock meson's normalized momentum-space radial wave Φ_L(p) on a p-grid."""
struct MockMomentumWave
    p::Vector{Float64}
    phi::Vector{Float64}
end

_mm_trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i = 1:length(x)-1)

"""
    mock_momentum_wave(radial, L; pmax=30.0, npoints=1501) -> MockMomentumWave

Spherical-Bessel transform of the reduced radial wave to momentum space,
normalized so `∫ p² Φ² dp = 1` (i.e. `Φ_L(p) = ∫ dr u(r) j_L(pr) √(2/π) p`, the
[`_momentum_radial_wave`](@ref) kernel).
"""
function mock_momentum_wave(radial::RadialWaveOnUniformMesh, L::Integer;
                            pmax::Real = 30.0, npoints::Integer = 1501)
    p = collect(range(0.0, float(pmax); length = npoints))
    phi = [_momentum_radial_wave(radial, pk, L) for pk in p]
    nrm = sqrt(_mm_trapz(p, p .^ 2 .* phi .^ 2))
    nrm > 0 || throw(ArgumentError("mock_momentum_wave: zero-norm wave"))
    return MockMomentumWave(p, phi ./ nrm)
end

"""Mean relativistic quark energy `⟨E⟩ = ∫ p² Φ² √(m²+p²) dp` over a mock wave."""
mock_mean_energy(mw::MockMomentumWave, m::Real) =
    _mm_trapz(mw.p, mw.p .^ 2 .* mw.phi .^ 2 .* sqrt.(float(m)^2 .+ mw.p .^ 2))

"""Mock mass `M̃ = ⟨E₁⟩ + ⟨E₂⟩` of a mock wave with constituent masses `m1, m2`."""
mock_wave_mass(mw::MockMomentumWave, m1::Real, m2::Real) =
    mock_mean_energy(mw, m1) + mock_mean_energy(mw, m2)

"""
    mock_meson_overlap(mwx, mwy, m_emit; Mx, My, exponent=0.7) -> I_i (GeV⁻¹)

Appendix-D M1 transition overlap `I_i(x,y)` for emitting quark mass `m_emit`,
with `Mx`, `My` the mock masses of the two states (see [`mock_wave_mass`](@ref)):

    I_i = √(4 Mx My)/(Mx + My) · ∫ dp p² Φₓ Φ_y (1/m_emit)(m_emit/E)^exponent .
"""
function mock_meson_overlap(mwx::MockMomentumWave, mwy::MockMomentumWave, m_emit::Real;
                            Mx::Real, My::Real, exponent::Real = 0.7)
    m = float(m_emit)
    pref = sqrt(4 * Mx * My) / (Mx + My)
    kern = @. mwx.p^2 * mwx.phi * mwy.phi * (1 / m) * (m / sqrt(m^2 + mwx.p^2))^exponent
    return pref * _mm_trapz(mwx.p, kern)
end

"""
    mock_meson_radial_moment(wx, wy, Ex, Ey, m_emit; n=1, exponent=0.5) -> Eₙⁱ (GeV⁻ⁿ)

Appendix-D E1/M2 radial moment `Eₙⁱ(x,y)` on the position-space reduced waves,
with `Ex`, `Ey` the mean quark energies (see [`mock_mean_energy`](@ref)):

    Eₙⁱ = |m_emit / √(Ex Ey)|^exponent · ∫ dr uₓ(r) u_y(r) rⁿ .
"""
function mock_meson_radial_moment(wx::RadialWaveOnUniformMesh, wy::RadialWaveOnUniformMesh,
                                  Ex::Real, Ey::Real, m_emit::Real;
                                  n::Integer = 1, exponent::Real = 0.5)
    wx.r == wy.r || throw(ArgumentError("mock_meson_radial_moment: meshes differ"))
    radial = sum(@. wx.u * wy.u * wx.r^n) * wx.h
    return abs(float(m_emit) / sqrt(Ex * Ey))^exponent * radial
end
