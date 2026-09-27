# Generic momentum-wave representations and overlap operations.

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

"""
    wave_mean_squares(wave, L; momentum_kwargs...) -> (r2, p2)
    wave_mean_squares(wave::OscillatorWave; momentum_kwargs...) -> (r2, p2)

Return the mean-square relative radius `r2 = <r^2>` in `GeV^-2` and
mean-square relative momentum `p2 = <p^2>` in `GeV^2` for one radial wave.
The calculation composes the representation-independent [`radial_expect`](@ref),
[`momentum_wave`](@ref), and [`momentum_expect`](@ref) operations.

`L` is required for a [`MeshWave`](@ref), whose sampled reduced radial function
does not itself store the orbital label. An [`OscillatorWave`](@ref) already
stores `L`, so the one-argument form is available. Keywords such as `pmax` and
`npoints` are forwarded to the numerical momentum transform of a mesh wave.

## Example

```julia
moments = wave_mean_squares(wave, 0)
moments.r2
moments.p2
```
"""
function wave_mean_squares(wave::RadialWave, L::Integer; momentum_kwargs...)
    L >= 0 || throw(ArgumentError("wave_mean_squares: L must be non-negative"))
    radial = radial_expect(wave, r -> r^2)
    momentum = momentum_expect(
        momentum_wave(wave, L; momentum_kwargs...),
        p -> p^2,
    )
    return (r2 = radial, p2 = momentum)
end

wave_mean_squares(wave::OscillatorWave; momentum_kwargs...) =
    wave_mean_squares(wave, wave.L; momentum_kwargs...)

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
