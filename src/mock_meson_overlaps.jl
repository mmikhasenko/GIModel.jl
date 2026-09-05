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

"""
    abstract type MomentumWave

A radial wave in momentum space. Separate from [`RadialWave`](@ref) because the
two spaces are different objects — *except* in the oscillator basis, where the
Fourier–Bessel transform maps the family onto itself (scale `β → 1/β`), so an
oscillator wave transformed is still an oscillator wave. A mesh wave transformed
is samples on a different grid, hence [`MeshMomentumWave`](@ref).
"""
abstract type MomentumWave end

"""A normalized momentum-space radial wave Φ_L(p) sampled on a p-grid."""
struct MeshMomentumWave <: MomentumWave
    p::Vector{Float64}
    phi::Vector{Float64}
end

"""Exact momentum-space view of a native [`OscillatorWave`](@ref)."""
struct OscillatorMomentumWave <: MomentumWave
    source::OscillatorWave
end

_mm_trapz(x, y) = sum(0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i]) for i = 1:length(x)-1)

function _spherical_bessel_j(L::Integer, x::Real)
    j0 = abs(x) < 1.0e-8 ? 1.0 - x^2 / 6 : sin(x) / x
    L == 0 && return j0
    abs(x) < 1.0e-4 && return x^L / prod(1:2:(2L+1))
    jm1 = j0
    j = sin(x) / x^2 - cos(x) / x
    for l in 1:(L-1)
        jp1 = (2l + 1) / x * j - jm1
        jm1, j = j, jp1
    end
    return j
end

function _momentum_radial_wave(radial::MeshWave, p::Real, L::Integer)
    accum = sum(
        radial.r[k] * radial.u[k] * _spherical_bessel_j(L, p * radial.r[k]) for
        k in eachindex(radial.r)
    )
    return sqrt(2 / π) * accum * radial.h
end

"""
    mock_momentum_wave(radial, L; pmax=30.0, npoints=1501) -> MeshMomentumWave

Spherical-Bessel transform of the reduced radial wave to momentum space,
normalized so `∫ p² Φ² dp = 1` (i.e. `Φ_L(p) = ∫ dr u(r) j_L(pr) √(2/π) p`, the
[`_momentum_radial_wave`](@ref) kernel).
"""
function mock_momentum_wave(radial::MeshWave, L::Integer;
                            pmax::Real = 30.0, npoints::Integer = 1501)
    p = collect(range(0.0, float(pmax); length = npoints))
    phi = [_momentum_radial_wave(radial, pk, L) for pk in p]
    nrm = sqrt(_mm_trapz(p, p .^ 2 .* phi .^ 2))
    nrm > 0 || throw(ArgumentError("mock_momentum_wave: zero-norm wave"))
    return MeshMomentumWave(p, phi ./ nrm)
end

"""Mean relativistic quark energy `⟨E⟩ = ∫ p² Φ² √(m²+p²) dp` over a mock wave."""
mock_mean_energy(mw::MomentumWave, m::Real) =
    momentum_expect(mw, p -> sqrt(float(m)^2 + p^2))

"""Mock mass `M̃ = ⟨E₁⟩ + ⟨E₂⟩` of a mock wave with constituent masses `m1, m2`."""
mock_wave_mass(mw::MomentumWave, m1::Real, m2::Real) =
    mock_mean_energy(mw, m1) + mock_mean_energy(mw, m2)

"""
    mock_meson_overlap(mwx, mwy, m_emit; Mx, My, exponent=0.7) -> I_i (GeV⁻¹)

Appendix-D M1 transition overlap `I_i(x,y)` for emitting quark mass `m_emit`,
with `Mx`, `My` the mock masses of the two states (see [`mock_wave_mass`](@ref)):

    I_i = √(4 Mx My)/(Mx + My) · ∫ dp p² Φₓ Φ_y (1/m_emit)(m_emit/E)^exponent .
"""
function mock_meson_overlap(mwx::MomentumWave, mwy::MomentumWave, m_emit::Real;
                            Mx::Real, My::Real, exponent::Real = 0.7)
    m = float(m_emit)
    pref = sqrt(4 * Mx * My) / (Mx + My)
    radial = momentum_overlap(
        mwx,
        mwy,
        p -> (1 / m) * (m / sqrt(m^2 + p^2))^exponent,
    )
    return pref * radial
end

"""
    mock_meson_radial_moment(wx, wy, Ex, Ey, m_emit; n=1, exponent=0.5) -> Eₙⁱ (GeV⁻ⁿ)

Appendix-D E1/M2 radial moment `Eₙⁱ(x,y)` on the position-space reduced waves,
with `Ex`, `Ey` the mean quark energies (see [`mock_mean_energy`](@ref)):

    Eₙⁱ = |m_emit / √(Ex Ey)|^exponent · ∫ dr uₓ(r) u_y(r) rⁿ .
"""
function mock_meson_radial_moment(wx::RadialWave, wy::RadialWave,
                                  Ex::Real, Ey::Real, m_emit::Real;
                                  n::Integer = 1, exponent::Real = 0.5)
    # This kernel IS `radial_overlap` with f(r) = r^n, times the Appendix-D
    # energy prefactor. Going through the interface normalizes both waves, which
    # this used to leave to the caller without saying so -- the same unstated
    # requirement that `_rel_momentum_average` carried. For the normalized waves
    # every solve now returns, the divisor is 1 and no number moves.
    radial = radial_overlap(wx, wy, x -> x^n)
    return abs(float(m_emit) / sqrt(Ex * Ey))^exponent * radial
end

# --- Radiative transition assembly (Eq. 22, Table VI) ------------------------
# The kernels above are the Appendix-D matrix elements; these three assemble
# them into the multipole amplitudes the paper tabulates. Promoted out of the
# Table VI audit script, where they sat next to the paper-comparison rows.

"""Proton mass in GeV — the unit of the M1 moments, which Table VI prints as μ/μ_N."""
const NUCLEON_MASS_GEV = 0.93827

"""
    photon_momentum(M_parent_GeV, M_child_GeV) -> Float64

Photon momentum `q = (M² - M'²) / 2M` for the radiative transition
`parent → child + γ`, in GeV.
"""
photon_momentum(M_parent_GeV::Real, M_child_GeV::Real) =
    (float(M_parent_GeV)^2 - float(M_child_GeV)^2) / (2 * float(M_parent_GeV))

"""
    m1_transition_moment(singlet, triplet, m1_GeV, m2_GeV, terms) -> Float64

The M1 transition moment `μ/μ_N` for `V → P γ` between two mock mesons given by
their momentum-space waves, in Table VI's units of `e/2`.

`terms` is an iterable of `(coefficient, m_i)` pairs — one per contributing
quark line — summed as `Σ c·I_i`. The coefficients follow the transition-moment
rule `μ = e_q I_q - e_q̄ I_q̄` (the antiquark charge enters flipped), so
charmonium is `+4/3 I_c` and bottomonium `-2/3 I_b`.
"""
function m1_transition_moment(
    singlet::MomentumWave,
    triplet::MomentumWave,
    m1_GeV::Real,
    m2_GeV::Real,
    terms,
)
    Mx = mock_wave_mass(singlet, m1_GeV, m2_GeV)
    My = mock_wave_mass(triplet, m1_GeV, m2_GeV)
    return sum(
        c * mock_meson_overlap(singlet, triplet, m_i; Mx = Mx, My = My) for (c, m_i) in terms
    ) * NUCLEON_MASS_GEV
end

"""
    e1_transition_amplitude(wave_S, mom_S, wave_P, mom_P, m_i, coeff_of_q,
                            M_parent_GeV, M_child_GeV; q = nothing) -> Float64

The E1 (and M2) amplitude `coeff(q) · E₁ⁱ · √(α q)` in `MeV^(1/2)`, with `q` in
MeV inside the square root — Table VI's amplitude convention.

`coeff_of_q` is the row's angular/charge factor as a function of `q` in GeV, so
the paper's printed q-dependence stays visible at the call site. Pass an
explicit `q` (GeV) to override the photon momentum implied by the two masses,
e.g. to use model rather than measured masses.
"""
function e1_transition_amplitude(
    wave_S::RadialWave,
    mom_S::MomentumWave,
    wave_P::RadialWave,
    mom_P::MomentumWave,
    m_i::Real,
    coeff_of_q,
    M_parent_GeV::Real,
    M_child_GeV::Real;
    q = nothing,
)
    q_GeV = isnothing(q) ? photon_momentum(M_parent_GeV, M_child_GeV) : float(q)
    E1 = mock_meson_radial_moment(
        wave_S, wave_P,
        mock_mean_energy(mom_S, m_i), mock_mean_energy(mom_P, m_i), m_i; n = 1,
    )
    return coeff_of_q(q_GeV) * E1 * sqrt(ALPHA_EM * 1000 * q_GeV)
end


# --- RadialWave interface: the momentum-space and origin operations ----------
# MeshWave implements these by numerical transform on the mesh. The oscillator
# implementation will not: for oscillator functions the Fourier-Bessel transform
# is exact (scale beta -> 1/beta) and the origin value is closed form.

"""
    momentum_wave(w::RadialWave, L; pmax=30.0, npoints=1501) -> MeshMomentumWave

The momentum-space radial wave `Phi(p)`, normalized so `integral p^2 Phi^2 dp = 1`.
For a [`MeshWave`](@ref) this is a numerical spherical-Bessel transform.
"""
momentum_wave(w::MeshWave, L::Integer; pmax::Real = 30.0, npoints::Integer = 1501) =
    mock_momentum_wave(w, L; pmax = pmax, npoints = npoints)

function momentum_wave(w::OscillatorWave, L::Integer = w.L; kwargs...)
    L == w.L || throw(ArgumentError(
        "momentum_wave: requested L=$L for OscillatorWave with L=$(w.L)",
    ))
    return OscillatorMomentumWave(w)
end

function _oscillator_momentum_value(mw::OscillatorMomentumWave, p::Real)
    w = mw.source
    pf = float(p)
    if pf == 0
        w.L > 0 && return 0.0
        pf = sqrt(eps(Float64)) * w.beta
    end
    return sum(
        w.coefficients[n + 1] * (-1.0)^n *
        ho_reduced_radial(n, w.L, inv(w.beta), pf) / pf for
        n in 0:(length(w.coefficients)-1)
    )
end

"""
    momentum_expect(mw::MeshMomentumWave, g) -> Float64

`integral p^2 Phi(p)^2 g(p) dp`. `g` is called as `g(p)`; `g = p -> sqrt(m^2+p^2)`
gives the mean relativistic quark energy.
"""
momentum_expect(mw::MeshMomentumWave, g) =
    _mm_trapz(mw.p, mw.p .^ 2 .* mw.phi .^ 2 .* map(g, mw.p))

function momentum_expect(mw::OscillatorMomentumWave, g)
    w = mw.source
    p2 = Symmetric(Matrix(ho_p2_matrix(w.L, w.beta, length(w.coefficients))))
    fact = eigen(p2)
    values = map(x -> g(sqrt(max(x, 0.0))), fact.values)
    op = fact.vectors * Diagonal(values) * fact.vectors'
    return dot(w.coefficients, op * w.coefficients) / wave_norm(w)
end

"""`integral p^2 Phi_x(p) Phi_y(p) g(p) dp`, with both waves normalized."""
function momentum_overlap(left::MeshMomentumWave, right::MeshMomentumWave, g)
    length(left.p) == length(right.p) || throw(ArgumentError(
        "momentum_overlap: momentum grids have different sizes",
    ))
    all(isapprox.(left.p, right.p; rtol = 1e-12, atol = 1e-14)) ||
        throw(ArgumentError("momentum_overlap: momentum grids differ"))
    nl = momentum_expect(left, _ -> 1.0)
    nr = momentum_expect(right, _ -> 1.0)
    (nl > 0 && nr > 0) || throw(ArgumentError("momentum_overlap: zero-norm wave"))
    value = _mm_trapz(
        left.p,
        left.p .^ 2 .* left.phi .* right.phi .* map(g, left.p),
    )
    return value / sqrt(nl * nr)
end


function momentum_overlap(
    left::OscillatorMomentumWave,
    right::OscillatorMomentumWave,
    g,
)
    wl, wr = left.source, right.source
    if wl.L == wr.L && wl.beta == wr.beta
        n = max(length(wl.coefficients), length(wr.coefficients))
        cl = vcat(wl.coefficients, zeros(n - length(wl.coefficients)))
        cr = vcat(wr.coefficients, zeros(n - length(wr.coefficients)))
        p2 = Symmetric(Matrix(ho_p2_matrix(wl.L, wl.beta, n)))
        fact = eigen(p2)
        values = map(x -> g(sqrt(max(x, 0.0))), fact.values)
        op = fact.vectors * Diagonal(values) * fact.vectors'
        return dot(cl, op * cr) / sqrt(wave_norm(wl) * wave_norm(wr))
    end
    pmax = max(
        _oscillator_momentum_cutoff(wl),
        _oscillator_momentum_cutoff(wr),
    )
    value, _ = quadgk(
        p -> p^2 * _oscillator_momentum_value(left, p) *
             _oscillator_momentum_value(right, p) * g(p),
        0.0,
        pmax;
        rtol = 1e-9,
    )
    return value / sqrt(wave_norm(wl) * wave_norm(wr))
end

"""
    momentum_functional(mw::MomentumWave, K) -> Float64

`integral p^2 Phi(p) K(p) dp` — **linear** in `Phi`, unlike
[`momentum_expect`](@ref) which is quadratic. `K` is supplied by the physics
that needs it, not by the wave: the annihilation amplitude `S_L` uses
`K(p) = (m/E) (p/E)^L`, which is why the mass belongs to the kernel and not to
the wave interface.
"""
momentum_functional(mw::MeshMomentumWave, K) =
    _mm_trapz(mw.p, mw.p .^ 2 .* mw.phi .* map(K, mw.p))

function momentum_functional(mw::OscillatorMomentumWave, K)
    pmax = _oscillator_momentum_cutoff(mw.source)
    value, _ = quadgk(
        p -> p^2 * _oscillator_momentum_value(mw, p) * K(p),
        0.0,
        pmax;
        rtol = 1e-9,
    )
    return value
end
