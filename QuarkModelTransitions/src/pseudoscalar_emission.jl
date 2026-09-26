# =============================================================================
# Solver-native GI pseudoscalar emission (Eq. 19): spatial integral layer
# =============================================================================

"""
    PseudoscalarEmission(g, h, quark_masses)

Godfrey--Isgur elementary pseudoscalar-emission operator of Eq. (19).
`g` and `h` are the dimensionless direct and recoil couplings.
`quark_masses` is a [`QuarkMassTable`](@ref) and is copied into the operator so
the constituent-coordinate momentum routing is explicit and reproducible.

This is distinct from the legacy fitted `A`/`S0` Table-IV model. The relation
`A = (g + h/4) beta` holds only for the paper's equal-mass, equal-`beta`,
single-oscillator reduction.
"""
struct PseudoscalarEmission <: StrongDecayOperator
    g::Float64
    h::Float64
    quark_masses::QuarkMassTable
    function PseudoscalarEmission(
        g::Real,
        h::Real,
        quark_masses::QuarkMassTable,
    )
        gf, hf = Float64(g), Float64(h)
        isfinite(gf) || throw(ArgumentError("pseudoscalar coupling g must be finite"))
        isfinite(hf) || throw(ArgumentError("pseudoscalar coupling h must be finite"))
        all(isfinite(mass) && mass > 0 for mass in values(quark_masses)) ||
            throw(ArgumentError("constituent masses must be finite and positive"))
        return new(gf, hf, copy(quark_masses))
    end
end

function _constituent_mass(operator::PseudoscalarEmission, flavor::Symbol)
    canonical = flavor === :n ? :q : flavor
    key = String(canonical)
    haskey(operator.quark_masses, key) || throw(ArgumentError(
        "PseudoscalarEmission has no constituent mass for flavor :$flavor",
    ))
    return operator.quark_masses[key]
end

"""Positive coefficient alpha in exp(sign*i*alpha*q.r) for one emitter."""
function _eq19_momentum_fraction(
    operator::PseudoscalarEmission,
    topology::_EmissionTopology,
    parent_flavors::Tuple{Symbol,Symbol},
)
    m_quark = _constituent_mass(operator, parent_flavors[1])
    m_antiquark = _constituent_mass(operator, parent_flavors[2])
    denominator = m_quark + m_antiquark
    # The spectator mass multiplies the relative coordinate of the emitter.
    return topology isa _QuarkEmission ?
           m_antiquark / denominator : m_quark / denominator
end

function _spherical_bessel_j0(x::Number)
    abs(x) < 1e-4 && return one(x) - x^2 / 6 + x^4 / 120 - x^6 / 5040
    return sin(x) / x
end

function _spherical_bessel_j1(x::Number)
    abs(x) < 1e-3 && return x / 3 - x^3 / 30 + x^5 / 840 - x^7 / 45360
    return sin(x) / x^2 - cos(x) / x
end

function _odd_double_factorial(n::Integer)
    n >= -1 || throw(ArgumentError("double-factorial argument must be >= -1"))
    result = 1
    for value in 1:2:n
        result *= value
    end
    return result
end

"""Spherical Bessel `j_l(x)` for real or complex `x` and integer `l >= 0`."""
function _spherical_bessel_j(order::Integer, x::Number)
    ell = Int(order)
    ell >= 0 || throw(ArgumentError("spherical-Bessel order must be non-negative"))
    ell == 0 && return _spherical_bessel_j0(x)
    ell == 1 && return _spherical_bessel_j1(x)
    iszero(x) && return zero(x)
    if abs(x) < ell + 1
        term = x^ell / _odd_double_factorial(2ell + 1)
        result = term
        for k in 0:255
            term *= -x^2 / (2 * (k + 1) * (2ell + 2k + 3))
            result_next = result + term
            abs(term) <= 8eps(Float64) * max(abs(result_next), 1.0) &&
                return result_next
            result = result_next
        end
        throw(ErrorException("spherical-Bessel series failed to converge"))
    end
    previous = _spherical_bessel_j0(x)
    current = _spherical_bessel_j1(x)
    for l in 1:(ell-1)
        previous, current = current, (2l + 1) * current / x - previous
    end
    return current
end

_spherical_phase(index::Integer) = isodd(abs(index)) ? -1 : 1

"""
Coefficient of one `j_order(alpha*q*r)` after expanding
`exp(sign*i*alpha*q*z)` between two Condon--Shortley spherical harmonics.
"""
function _plane_wave_angular_coefficient(
    target_L::Integer,
    target_m::Integer,
    source_L::Integer,
    source_m::Integer,
    order::Integer,
    sign::Integer,
)
    Lt, mt = Int(target_L), Int(target_m)
    Ls, ms = Int(source_L), Int(source_m)
    ell = Int(order)
    sign in (-1, 1) || throw(ArgumentError("plane-wave sign must be +/-1"))
    (Lt >= 0 && Ls >= 0 && ell >= 0) || return 0.0 + 0.0im
    (abs(mt) <= Lt && abs(ms) <= Ls && mt == ms) || return 0.0 + 0.0im
    abs(Lt - Ls) <= ell <= Lt + Ls || return 0.0 + 0.0im
    c0 = ComplexF64(CG(ell, 0, Ls, 0, Lt, 0))
    cm = ComplexF64(CG(ell, 0, Ls, ms, Lt, mt))
    return (sign * im)^ell * (2ell + 1) * sqrt((2Ls + 1) / (2Lt + 1)) *
           c0 * cm
end

function _bessel_overlap(
    final_wave::RadialWave,
    initial_wave::RadialWave,
    order::Int,
    momentum,
)
    return radial_overlap(
        final_wave,
        initial_wave,
        r -> _spherical_bessel_j(order, momentum * r),
    )
end

function _bessel_over_r(order::Int, momentum, r::Real)
    iszero(r) && return zero(momentum)
    return _spherical_bessel_j(order, momentum * r) / r
end

function _gradient_radial_overlap(
    final_wave::RadialWave,
    initial_wave::RadialWave,
    order::Int,
    momentum,
    inverse_r_coefficient::Int,
)
    bessel(r) = _spherical_bessel_j(order, momentum * r)
    derivative = GIModel.radial_derivative_overlap(final_wave, initial_wave, bessel)
    iszero(inverse_r_coefficient) && return derivative
    inverse_r = radial_overlap(
        final_wave,
        initial_wave,
        r -> _bessel_over_r(order, momentum, r),
    )
    return derivative + inverse_r_coefficient * inverse_r
end

function _validate_orbital_wave(wave::OscillatorWave, L::Int, role::AbstractString)
    wave.L == L || throw(ArgumentError(
        "$role OscillatorWave has L=$(wave.L), but the orbital label requires L=$L",
    ))
    return nothing
end

_validate_orbital_wave(::RadialWave, ::Int, ::AbstractString) = nothing

function _direct_orbital_integral(
    operator::PseudoscalarEmission,
    label::_SpectroscopicOrbitalLabel,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    q::Number,
    alpha::Real,
)
    elementary = label.elementary
    iszero(elementary.vector_component) || return zero(complex(float(q)))
    total = zero(complex(float(q)))
    for order in abs(label.final_L-label.initial_L):(label.final_L+label.initial_L)
        angular = _plane_wave_angular_coefficient(
            label.final_L,
            label.final_mL,
            label.initial_L,
            label.initial_mL,
            order,
            elementary.plane_wave_sign,
        )
        iszero(angular) && continue
        total += angular * _bessel_overlap(
            final_wave, initial_wave, order, alpha * q,
        )
    end
    return operator.g * q * total
end

function _recoil_gradient_branch(
    label::_SpectroscopicOrbitalLabel,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    q::Number,
    alpha::Real,
    gradient_L::Int,
)
    Lf, mf = label.final_L, label.final_mL
    Li, mi = label.initial_L, label.initial_mL
    component = label.elementary.vector_component
    gradient_component = -component
    gradient_m = mf + gradient_component
    abs(gradient_m) <= gradient_L || return zero(complex(float(q)))

    if gradient_L == Lf + 1
        gradient_coefficient = sqrt((Lf + 1) / (2Lf + 3)) * ComplexF64(CG(
            Lf, mf, 1, gradient_component, gradient_L, gradient_m,
        ))
        inverse_r_coefficient = -(Lf + 1)
    elseif Lf > 0 && gradient_L == Lf - 1
        gradient_coefficient = -sqrt(Lf / (2Lf - 1)) * ComplexF64(CG(
            Lf, mf, 1, gradient_component, gradient_L, gradient_m,
        ))
        inverse_r_coefficient = Lf
    else
        return zero(complex(float(q)))
    end
    iszero(gradient_coefficient) && return zero(complex(float(q)))

    total = zero(complex(float(q)))
    for order in abs(gradient_L-Li):(gradient_L+Li)
        angular = _plane_wave_angular_coefficient(
            gradient_L,
            gradient_m,
            Li,
            mi,
            order,
            label.elementary.plane_wave_sign,
        )
        iszero(angular) && continue
        radial = _gradient_radial_overlap(
            final_wave,
            initial_wave,
            order,
            alpha * q,
            inverse_r_coefficient,
        )
        total += angular * radial
    end
    return gradient_coefficient * total
end

function _recoil_orbital_integral(
    operator::PseudoscalarEmission,
    label::_SpectroscopicOrbitalLabel,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    q::Number,
    alpha::Real,
)
    component = label.elementary.vector_component
    component in -1:1 || throw(ArgumentError("vector component must be -1, 0, or +1"))
    upper = _recoil_gradient_branch(
        label, final_wave, initial_wave, q, alpha, label.final_L + 1,
    )
    lower = label.final_L > 0 ? _recoil_gradient_branch(
        label, final_wave, initial_wave, q, alpha, label.final_L - 1,
    ) : zero(upper)
    # This is the spherical component of -i*gradient acting to the left. The
    # (-1)^component is the Hermitian spherical-tensor conjugation phase.
    return operator.h * (-im) * _spherical_phase(component) * (upper + lower)
end

"""
Evaluate an S-wave-to-S-wave Eq. (19) spatial column in the helicity frame.

The direct column is `g*q*<j0(alpha*q*r)>`. The recoil column evaluates the
radial part of `sigma.p'`, including the `d/dr-u/r` gradient of the final
S-wave and the sign of the constituent plane wave. Couplings are included;
spin, flavor, topology, and spherical-contraction coefficients are not.
"""
function _eq19_s_wave_spatial_integral(
    operator::PseudoscalarEmission,
    label::_Eq19OrbitalLabel,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    momentum_GeV::Number,
    parent_flavors::Tuple{Symbol,Symbol},
)
    spectroscopic = _SpectroscopicOrbitalLabel(label, 0, 0, 0, 0)
    return _eq19_spatial_integral(
        operator,
        spectroscopic,
        final_wave,
        initial_wave,
        momentum_GeV,
        parent_flavors,
    )
end

function _eq19_spatial_integral(
    operator::PseudoscalarEmission,
    label::_SpectroscopicOrbitalLabel,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    momentum_GeV::Number,
    parent_flavors::Tuple{Symbol,Symbol},
)
    Li, mi = label.initial_L, label.initial_mL
    Lf, mf = label.final_L, label.final_mL
    (Li >= 0 && abs(mi) <= Li) || throw(ArgumentError(
        "invalid initial orbital labels L=$Li, m=$mi",
    ))
    (Lf >= 0 && abs(mf) <= Lf) || throw(ArgumentError(
        "invalid final orbital labels L=$Lf, m=$mf",
    ))
    _validate_orbital_wave(initial_wave, Li, "initial")
    _validate_orbital_wave(final_wave, Lf, "final")
    alpha = _eq19_momentum_fraction(
        operator, label.elementary.topology, parent_flavors,
    )
    if label.elementary.piece isa _DirectPseudoscalarPiece
        return _direct_orbital_integral(
            operator, label, final_wave, initial_wave, momentum_GeV, alpha,
        )
    end
    return _recoil_orbital_integral(
        operator, label, final_wave, initial_wave, momentum_GeV, alpha,
    )
end

"""Numerical values over exactly the shared columns produced in Phase III."""
struct _Eq19WaveAmplitude{I<:Tuple,H<:Tuple,W<:Tuple}
    spatial_integrals::I
    helicity::H
    partial_waves::W
end

function _eq19_wave_amplitude(
    operator::PseudoscalarEmission,
    decomposition::_AngularCoefficientDecomposition,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    momentum_GeV::Number,
    parent_flavors::Tuple{Symbol,Symbol},
)
    integrals = Tuple(
        _eq19_spatial_integral(
            operator,
            label,
            final_wave,
            initial_wave,
            momentum_GeV,
            parent_flavors,
        ) for label in decomposition.integral_labels
    )
    values = collect(integrals)
    helicity_values = decomposition.helicity_coefficients * values
    partial_wave_values = decomposition.partial_wave_coefficients * values
    helicity = Tuple(
        label => value for (label, value) in
        zip(decomposition.helicities, helicity_values)
    )
    partial_waves = Tuple(
        label => value for (label, value) in
        zip(decomposition.partial_waves, partial_wave_values)
    )
    return _Eq19WaveAmplitude(integrals, helicity, partial_waves)
end
