# Photon emission from the nonrelativistic quark current, exact in the photon
# momentum and resolved into every multipole.
#
# Per constituent i (charge ê, mass m, position r_i = R + s_i w_i r, with
# s = +1, -1 and w = m_other/M for quark and antiquark) the emission operator is
#
#   -(ê/m) e^{-iq.r_i} eps*.p_i + (iê/2m) sigma_i.(q x eps*) e^{-iq.r_i},
#
# the convection and magnetization currents of H_I = -(ê/m) A.p - (ê/2m) sigma.B.
# With q = q z-hat the plane wave is expanded in spherical Bessel functions and
# every angular integral is a Clebsch-Gordan product, so no multipole is
# truncated. The width is Gamma = 2 alpha q (2J_i+1)^-1 sum |M_lambda|^2.

"""
    MultipolePhotonEmission(quark_masses; electric_form = :siegert)

Photon emission evaluated from the quark convection and magnetization
currents, exact in the photon momentum and resolved into all multipoles
(E1, M1, E2, M2, E3, ...). It applies to any pair of `PhysicalState`s with
shared flavors, including D and F waves and mixed states. Unlike
[`PhotonEmission`](@ref), which evaluates the leading multipole in the 1985
mock-meson prescription with its `(m/E)` smearing, this operator has no
fitted factors.

`electric_form = :siegert` (default) takes electric multipoles from the
charge density through Siegert's theorem with `omega = q`, as GI and later
GI-model papers do for E1; their leading E1 is reproduced exactly. The
magnetization current's electric multipoles (the spin-flip E1, and the
O(q/m) spin-orbit-like correction to allowed E1) are always in current form.
`electric_form = :current` uses the convection current for every multipole;
for relativized wavefunctions it differs from the Siegert form, which is a
measure of the gauge ambiguity of the model.

`matrix_element` returns a result whose `multipoles` field lists
`(multipole, k, amplitude)` with amplitudes in MeV^(1/2), so that
`decay_width` is their sum of squares. Use [`multipole_fractions`](@ref) for
the normalized amplitudes measured in charmonium.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
masses = QuarkMassTable("c" => 1.628)
state(L, S, J, wave, M) = PhysicalState("cc", M, [(
    basis = BasisState(1, L, 2S + 1, J; flavors = (:c, :c)),
    coefficient = 1.0, wave = wave,
)])
chi_c2 = state("P", 1, 2, OscillatorWave(1, 0.6, [1.0]), 3.556)
psi = state("S", 1, 1, OscillatorWave(0, 0.7, [1.0]), 3.097)
amplitude = matrix_element(psi, MultipolePhotonEmission(masses), chi_c2)
[m.multipole for m in amplitude.multipoles]      # [:E1, :M2, :E3]
@assert decay_width(amplitude) > 0               # MeV
```

## Related

- [`multipole_fractions`](@ref) — normalized multipole amplitudes.
- [`PhotonEmission`](@ref) — leading multipole in the GI 1985 prescription.
- [`decay_width`](@ref) — width in MeV.
- [`matrix_element`](@ref) — evaluate an operator between states.
"""
struct MultipolePhotonEmission <: TransitionOperator
    quark_masses::QuarkMassTable
    electric_form::Symbol
    function MultipolePhotonEmission(quark_masses::QuarkMassTable; electric_form::Symbol = :siegert)
        electric_form in (:siegert, :current) || throw(ArgumentError(
            "electric_form must be :siegert or :current",
        ))
        all(isfinite(m) && m > 0 for m in values(quark_masses)) ||
            throw(ArgumentError("constituent masses must be finite and positive"))
        return new(copy(quark_masses), electric_form)
    end
end

"""Multipole-resolved photon amplitude; `multipoles` amplitudes are in MeV^(1/2)."""
struct PhotonMultipoleAmplitude{O,I,F,P}
    operator::O
    initial::I
    final::F
    momentum_GeV::Float64
    multipoles::Vector{NamedTuple{(:multipole, :k, :amplitude),Tuple{Symbol,Int,Float64}}}
    provenance::P
end

"""
    multipole_fractions(amplitude)

Normalized multipole amplitudes `a_k = A_k / sqrt(sum A^2)`, keyed by
multipole (`:E1`, `:M2`, ...), with the sign convention that the lowest
multipole is positive. For a decay `J_i -> J_f gamma` this is the convention
of Karl, Meshkov and Rosner and of the CLEO and BESIII analyses of
`chi_cJ -> J/psi gamma` (`a_2` is the M2 fraction). Their amplitudes `b_k`
for `psi(2S) -> gamma chi_cJ` put chi_cJ, the final state, in the angular-momentum
slot and equal `(-1)^(k-1)` times these.

## Example

```julia
using GIModel, GIModel.QuarkModelTransitions
masses = QuarkMassTable("c" => 1.628)
state(L, S, J, wave, M) = PhysicalState("cc", M, [(
    basis = BasisState(1, L, 2S + 1, J; flavors = (:c, :c)),
    coefficient = 1.0, wave = wave,
)])
chi_c1 = state("P", 1, 1, OscillatorWave(1, 0.6, [1.0]), 3.511)
psi = state("S", 1, 1, OscillatorWave(0, 0.7, [1.0]), 3.097)
fractions = multipole_fractions(matrix_element(psi, MultipolePhotonEmission(masses), chi_c1))
fractions[:M2]       # about -E_gamma/(4 m_c)
@assert fractions[:M2] < 0
```

## Related

- [`MultipolePhotonEmission`](@ref) — the multipole-resolved photon operator.
"""
function multipole_fractions(amplitude::PhotonMultipoleAmplitude)
    entries = amplitude.multipoles
    norm = sqrt(sum(abs2(m.amplitude) for m in entries))
    norm > 0 || throw(ArgumentError("the transition has a vanishing amplitude"))
    reference = first(m for m in entries if !iszero(m.amplitude))
    sgn = sign(reference.amplitude)
    return Dict(m.multipole => sgn * m.amplitude / norm for m in entries)
end

decay_width(amplitude::PhotonMultipoleAmplitude) = sum(abs2(m.amplitude) for m in amplitude.multipoles)

function decay_width(final::PhysicalState, operator::MultipolePhotonEmission, initial::PhysicalState)
    initial.mass_GeV <= final.mass_GeV && return 0.0
    return decay_width(matrix_element(final, operator, initial))
end

function Base.show(io::IO, amplitude::PhotonMultipoleAmplitude)
    parts = join(("$(m.multipole)=$(round(m.amplitude; sigdigits = 4))" for m in amplitude.multipoles), ", ")
    print(io, "PhotonMultipoleAmplitude(", amplitude.initial.label, " -> ",
        amplitude.final.label, " + gamma, ", parts, " MeV^1/2)")
end

# --- angular algebra -----------------------------------------------------------

function _cg_or_zero(j1, m1, j2, m2, J, M)
    (abs(m1) <= j1 && abs(m2) <= j2 && abs(M) <= J && m1 + m2 == M &&
     abs(j1 - j2) <= J <= j1 + j2) || return 0.0
    return Float64(CG(j1, m1, j2, m2, J, M))
end

# Integral of conj(Y_{Lf mf}) Y_{l 0} Y_{L m} over the sphere.
function _gaunt_z(Lf, mf, l, L, m)
    mf == m || return 0.0
    return sqrt((2l + 1) * (2L + 1) / (4pi * (2Lf + 1))) *
           _cg_or_zero(L, 0, l, 0, Lf, 0) * _cg_or_zero(L, m, l, 0, Lf, mf)
end

# Spherical Bessel function for a signed argument: j_l(-x) = (-1)^l j_l(x).
_signed_bessel_j(l, x) = x < 0 ? (-1)^l * _spherical_bessel_j(l, -x) : _spherical_bessel_j(l, x)

# --- radial integrals ∫ u_f(r) K(r, u_i, u_i') dr --------------------------------

_photon_wave_value(w::OscillatorWave, r) = GIModel._oscillator_radial_value(w, r)
_photon_wave_derivative(w::OscillatorWave, r) = GIModel._oscillator_radial_derivative(w, r)

function _photon_radial(final::OscillatorWave, initial::OscillatorWave, kernel)
    rmax = max(GIModel._oscillator_coordinate_cutoff(final),
               GIModel._oscillator_coordinate_cutoff(initial))
    value, _ = GIModel.quadgk(0.0, rmax; rtol = 1e-10) do r
        _photon_wave_value(final, r) *
            kernel(r, _photon_wave_value(initial, r), _photon_wave_derivative(initial, r))
    end
    return value / sqrt(wave_norm(final) * wave_norm(initial))
end

function _photon_radial(final::MeshWave, initial::MeshWave, kernel)
    length(final.r) == length(initial.r) && isapprox(final.h, initial.h; rtol = 1e-10) &&
        all(isapprox.(final.r, initial.r; rtol = 1e-10, atol = 1e-12)) ||
        throw(ArgumentError("photon multipoles need both waves on the same mesh"))
    u, h, n = initial.u, initial.h, length(initial.u)
    total = 0.0
    # Dirichlet boundaries u = 0 one step beyond each end of the interior mesh.
    for i in 1:n
        left = i == 1 ? 0.0 : u[i - 1]
        right = i == n ? 0.0 : u[i + 1]
        total += final.u[i] * kernel(initial.r[i], u[i], (right - left) / (2h))
    end
    return total * h / sqrt(wave_norm(final) * wave_norm(initial))
end

_photon_radial(::RadialWave, ::RadialWave, _) = throw(ArgumentError(
    "photon multipoles need both waves in the same representation (oscillator or mesh)",
))

# Radial integrals of one component pair, memoized by (kind, l, kappa).
struct _PhotonRadialCache{F,I}
    final::F
    initial::I
    L_initial::Int
    values::Dict{Tuple{Symbol,Int,Float64},Float64}
end
_PhotonRadialCache(final, initial, L) = _PhotonRadialCache(final, initial, L, Dict{Tuple{Symbol,Int,Float64},Float64}())

function _radial_term(cache::_PhotonRadialCache, kind::Symbol, l::Int, kappa::Float64)
    get!(cache.values, (kind, l, kappa)) do
        L = cache.L_initial
        kernel = if kind === :plain
            (r, u, du) -> u * _signed_bessel_j(l, kappa * r)
        elseif kind === :raise   # d_{L+1} = u' - (L+1) u / r
            (r, u, du) -> _signed_bessel_j(l, kappa * r) * (du - (L + 1) * u / r)
        else                     # d_{L-1} = u' + L u / r
            (r, u, du) -> _signed_bessel_j(l, kappa * r) * (du + L * u / r)
        end
        _photon_radial(cache.final, cache.initial, kernel)
    end
end

# <Lf mf| e^{-i kappa z} |L m>
function _plane_wave_element(cache, Lf, mf, L, m, kappa)
    mf == m || return 0.0im
    total = 0.0im
    for l in abs(Lf - L):(Lf + L)
        g = _gaunt_z(Lf, mf, l, L, m)
        iszero(g) && continue
        total += (-im)^l * sqrt(4pi * (2l + 1)) * g * _radial_term(cache, :plain, l, kappa)
    end
    return total
end

# <Lf mf| e^{-i kappa z} (-i nabla_mu) |L m>, gradient formula of Varshalovich 7.3.
function _convection_element(cache, Lf, mf, L, m, mu, kappa)
    total = 0.0im
    for (Lp, c, kind) in ((L + 1, sqrt((L + 1) / (2L + 3)), :raise),
                          (L - 1, L == 0 ? 0.0 : -sqrt(L / (2L - 1)), :lower))
        (Lp < 0 || iszero(c)) && continue
        coupling = _cg_or_zero(L, m, 1, mu, Lp, m + mu)
        iszero(coupling) && continue
        for l in abs(Lf - Lp):(Lf + Lp)
            g = _gaunt_z(Lf, mf, l, Lp, m + mu)
            iszero(g) && continue
            total += -im * c * coupling * (-im)^l * sqrt(4pi * (2l + 1)) * g *
                     _radial_term(cache, kind, l, kappa)
        end
    end
    return total
end

# --- helicity amplitudes ----------------------------------------------------------

function _multipole_mass(operator::MultipolePhotonEmission, flavor::Symbol)
    key = String(flavor in (:u, :d, :n) ? :q : flavor)
    haskey(operator.quark_masses, key) || throw(ArgumentError(
        "MultipolePhotonEmission has no constituent mass for flavor :$flavor",
    ))
    return operator.quark_masses[key]
end

# Convection and magnetization parts of M_{+1}(Mi) and the density element
# rho(Mi) = <f Mi| sum_i ê_i e^{-i kappa_i z} |i Mi>, summed over component pairs.
function _photon_component_sums(final, operator, initial, q)
    Ji, Jf = initial.J, final.J
    conv = zeros(ComplexF64, 2Ji + 1)
    mag = zeros(ComplexF64, 2Ji + 1)
    rho = zeros(ComplexF64, 2Ji + 1)
    matched = false
    for parent in _electromagnetic_components(initial), daughter in _electromagnetic_components(final)
        flavors = parent.basis.flavors
        daughter.basis.flavors == flavors || continue
        matched = true
        m = (_multipole_mass(operator, flavors[1]), _multipole_mass(operator, flavors[2]))
        M = m[1] + m[2]
        charge = (_quark_charge(flavors[1]), -_quark_charge(flavors[2]))
        sgn = (1, -1)
        kappa = (m[2] / M * q, -m[1] / M * q)
        Li, Lf = orbital_angular_momentum(parent.basis.L_label), orbital_angular_momentum(daughter.basis.L_label)
        Si, Sf = (parent.basis.multiplicity - 1) ÷ 2, (daughter.basis.multiplicity - 1) ÷ 2
        mixing = parent.coefficient * conj(daughter.coefficient)
        cache = _PhotonRadialCache(daughter.wave, parent.wave, Li)
        for (index, Mi) in enumerate(-Ji:Ji)
            for mLi in -Li:Li, mSi in -Si:Si
                a = _cg_or_zero(Li, mLi, Si, mSi, Ji, Mi)
                iszero(a) && continue
                spin_i = _coupled_spin(Si, mSi)
                # photon helicity +1: M_f = Mi - 1
                for mLf in -Lf:Lf, mSf in -Sf:Sf
                    b = _cg_or_zero(Lf, mLf, Sf, mSf, Jf, Mi - 1)
                    iszero(b) && continue
                    pref = -mixing * a * b          # (-1)^lambda with lambda = +1
                    spin_f = _coupled_spin(Sf, mSf)
                    for i in 1:2
                        if Sf == Si && mSf == mSi
                            conv[index] += pref * (-charge[i] * sgn[i] / m[i]) *
                                _convection_element(cache, Lf, mLf, Li, mLi, -1, kappa[i])
                        end
                        topology = i == 1 ? _QuarkEmission() : _AntiquarkEmission()
                        s = _spin_factor(topology, spin_f, spin_i, -1)
                        iszero(s) && continue
                        mag[index] += pref * (-charge[i] * q / (2m[i])) * s *
                            _plane_wave_element(cache, Lf, mLf, Li, mLi, kappa[i])
                    end
                end
                # density, photon projection 0: M_f = Mi, spin untouched
                Sf == Si || continue
                b = _cg_or_zero(Lf, mLi, Sf, mSi, Jf, Mi)
                iszero(b) && continue
                for i in 1:2
                    rho[index] += mixing * a * b * charge[i] *
                        _plane_wave_element(cache, Lf, mLi, Li, mLi, kappa[i])
                end
            end
        end
    end
    matched || throw(ArgumentError("initial and final states share no flavor component"))
    return conv, mag, rho
end

# Wigner-Eckart projection of f(Mi) on multipoles k for photon projection mu:
# f(Mi) = sum_k c_k sqrt((2k+1)/(2Ji+1)) <k mu; Jf Mi-mu | Ji Mi>.
function _multipole_projection(values, Ji, Jf, mu)
    ks = max(abs(Ji - Jf), abs(mu)):(Ji + Jf)
    basis(k, Mi) = sqrt((2k + 1) / (2Ji + 1)) * _cg_or_zero(k, mu, Jf, Mi - mu, Ji, Mi)
    Ms = [Mi for Mi in -Ji:Ji if abs(Mi - mu) <= Jf]
    out = Dict(k => sum(values[Mi + Ji + 1] * basis(k, Mi) for Mi in Ms; init = 0.0im) for k in ks)
    residual = maximum((abs(values[Mi + Ji + 1] - sum(out[k] * basis(k, Mi) for k in ks)) for Mi in Ms); init = 0.0)
    scale = maximum(abs, values; init = 0.0)
    residual <= 1e-9 * max(scale, 1e-300) || error("multipole projection is not complete")
    return out
end

# Siegert: transverse electric multipole = -sqrt((k+1)/(2k)) (omega/q) Coulomb multipole.
_siegert_factor(k) = -sqrt((k + 1) / (2k))

function matrix_element(final::PhysicalState, operator::MultipolePhotonEmission, initial::PhysicalState)
    q = _on_shell_photon_momentum(final, initial)
    return _multipole_matrix_element(final, operator, initial, q)
end

function _multipole_matrix_element(final, operator, initial, q)
    Ji, Jf = initial.J, final.J
    (Ji == 0 && Jf == 0) && throw(ArgumentError("J = 0 -> J = 0 photon emission is forbidden"))
    conv, mag, rho = _photon_component_sums(final, operator, initial, q)
    cm = _multipole_projection(conv, Ji, Jf, 1)
    mm = _multipole_projection(mag, Ji, Jf, 1)
    dm = operator.electric_form === :siegert ? _multipole_projection(rho, Ji, Jf, 0) : nothing
    parity_change = initial.parity * final.parity == -1
    # All multipoles of one transition share the phase (-i)^l of the orbital
    # rank: real for parity-conserving, -i times real for parity-changing. The
    # overall sign is a convention; it is chosen so that P -> S E1 and S -> S M1
    # carry the sign of the leading `PhotonEmission` kernels.
    phase = parity_change ? -im : -1.0
    width_scale = sqrt(1000 * 4 * ALPHA_EM * q / (2Ji + 1))   # MeV^(1/2)
    entries = NamedTuple{(:multipole, :k, :amplitude),Tuple{Symbol,Int,Float64}}[]
    for k in sort(collect(keys(cm)))
        k == 0 && continue
        electric = (initial.parity * final.parity) == (-1)^k
        convection = (electric && !isnothing(dm)) ? _siegert_factor(k) * get(dm, k, 0.0im) : cm[k]
        value = (convection + mm[k]) * phase
        abs(imag(value)) <= 1e-8 * max(abs(value), 1e-300) + 1e-15 ||
            error("multipoles of one transition must share a common phase")
        push!(entries, (multipole = Symbol(electric ? "E" : "M", k), k = k,
                        amplitude = real(value) * width_scale))
    end
    # Multipoles forbidden for these components come out at rounding level.
    largest = maximum(abs(e.amplitude) for e in entries; init = 0.0)
    entries = [abs(e.amplitude) <= 1e-12 * largest ? (e..., amplitude = 0.0) : e for e in entries]
    provenance = (
        backend = :nonrelativistic_current,
        electric_form = operator.electric_form,
        omega = :photon_momentum,
        amplitude_units = :MeV_sqrt,
        width_units = :MeV,
    )
    return PhotonMultipoleAmplitude(operator, initial, final, q, entries, provenance)
end
