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
    label.vector_component == 0 || throw(ArgumentError(
        "S-to-S Eq. (19) spatial integrals vanish outside vector component zero",
    ))
    q = momentum_GeV
    alpha = _eq19_momentum_fraction(operator, label.topology, parent_flavors)
    kernel0(r) = _spherical_bessel_j0(alpha * q * r)
    if label.piece isa _DirectPseudoscalarPiece
        return operator.g * q * radial_overlap(final_wave, initial_wave, kernel0)
    end
    kernel1(r) = _spherical_bessel_j1(alpha * q * r)
    derivative = GIModel.radial_derivative_overlap(
        final_wave, initial_wave, kernel1,
    )
    angular = radial_overlap(
        final_wave, initial_wave, r -> kernel1(r) / r,
    )
    return operator.h * label.plane_wave_sign * (derivative - angular)
end

function _eq19_spatial_integral(
    operator::PseudoscalarEmission,
    label::_SpectroscopicOrbitalLabel,
    final_wave::RadialWave,
    initial_wave::RadialWave,
    momentum_GeV::Number,
    parent_flavors::Tuple{Symbol,Symbol},
)
    (label.initial_L, label.final_L) == (0, 0) || throw(ArgumentError(
        "the first native Eq. (19) spatial slice supports S-to-S waves; " *
        "L=$(label.initial_L) to L=$(label.final_L) is not implemented yet",
    ))
    (label.initial_mL, label.final_mL) == (0, 0) || throw(ArgumentError(
        "an S-wave orbital projection must have mL=0",
    ))
    return _eq19_s_wave_spatial_integral(
        operator,
        label.elementary,
        final_wave,
        initial_wave,
        momentum_GeV,
        parent_flavors,
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
