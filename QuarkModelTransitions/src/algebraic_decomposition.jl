# =============================================================================
# Algebraic coefficient / spatial-integral boundary for GI Eq. (19)
# =============================================================================

abstract type _ElementaryEmissionPiece end
struct _DirectPseudoscalarPiece <: _ElementaryEmissionPiece end
struct _RecoilPseudoscalarPiece <: _ElementaryEmissionPiece end

"""
Label one spatial tensor component required by the elementary-emission algebra.

`vector_component` is the spherical component of `q` or `p'`, while
`plane_wave_sign` records the sign in `exp(sign * i q' dot r)`.  Wave functions,
momentum, couplings, and numerical controls deliberately do not enter this
algebraic label.
"""
struct _Eq19OrbitalLabel{P<:_ElementaryEmissionPiece,T<:_EmissionTopology}
    piece::P
    topology::T
    vector_component::Int
    plane_wave_sign::Int
end

function orbital_decomposition(
    piece::_ElementaryEmissionPiece,
    topology::_EmissionTopology,
    spin_component::Integer,
)
    mu = Int(spin_component)
    mu in -1:1 || throw(ArgumentError("spherical spin component must be -1, 0, or +1"))
    route = topology isa _QuarkEmission ? -1 : 1
    return _Eq19OrbitalLabel(piece, topology, -mu, route)
end

# After factoring out the common `i` in GI Eq. (19), the lower signs give an
# antiquark minus for the direct g sigma.q term and two cancelling minuses for
# the recoil h sigma.p' term.
_eq19_topology_phase(::_DirectPseudoscalarPiece, ::_QuarkEmission) = 1
_eq19_topology_phase(::_DirectPseudoscalarPiece, ::_AntiquarkEmission) = -1
_eq19_topology_phase(::_RecoilPseudoscalarPiece, ::_QuarkEmission) = 1
_eq19_topology_phase(::_RecoilPseudoscalarPiece, ::_AntiquarkEmission) = 1

_eq19_coupling(::_DirectPseudoscalarPiece) = :g
_eq19_coupling(::_RecoilPseudoscalarPiece) = :h
_eq19_component_selection(::_DirectPseudoscalarPiece, mu::Int) = iszero(mu) ? 1 : 0
_eq19_component_selection(::_RecoilPseudoscalarPiece, ::Int) = 1

"""One inspectable algebraic multiplier of a named Eq. (19) spatial integral."""
struct _AlgebraicCoefficientTerm{P,T,F,S,C,L,R}
    piece::P
    topology::T
    spin_component::Int
    flavor_factor::F
    spin_factor::S
    topology_phase::Int
    spherical_phase::Int
    component_selection::Int
    coefficient::C
    orbital_label::L
    provenance::R
end

"""
    _eq19_coefficient_decomposition(piece, daughter_spin, parent_spin,
                                    daughter_flavor, emitted_flavor,
                                    parent_flavor)

Build the wave-function-independent part of one GI Eq. (19) operator piece.
The returned tuple retains both quark-line topologies and all three spherical
components, including exact zeros.  Each coefficient multiplies the spatial
integral named by `orbital_label`; no radial wave is evaluated here.
"""
function _eq19_coefficient_decomposition(
    piece::_ElementaryEmissionPiece,
    daughter_spin::_SpinState,
    parent_spin::_SpinState,
    daughter_flavor::_FlavorState,
    emitted_flavor::_FlavorState,
    parent_flavor::_FlavorState,
)
    topologies = (_QuarkEmission(), _AntiquarkEmission())
    return Tuple(
        begin
            flavor = _flavor_factor(
                topology, daughter_flavor, emitted_flavor, parent_flavor,
            )
            spin = _spin_factor(topology, daughter_spin, parent_spin, mu)
            topology_phase = _eq19_topology_phase(piece, topology)
            spherical_phase = isodd(abs(mu)) ? -1 : 1
            component_selection = _eq19_component_selection(piece, mu)
            coefficient = topology_phase * spherical_phase * component_selection *
                          flavor * spin
            label = orbital_decomposition(piece, topology, mu)
            provenance = (
                source = :GI1985_Eq19,
                coupling = _eq19_coupling(piece),
                contraction = :spherical_scalar_product,
                spin_component = mu,
                vector_component = -mu,
            )
            _AlgebraicCoefficientTerm(
                piece,
                topology,
                mu,
                flavor,
                spin,
                topology_phase,
                spherical_phase,
                component_selection,
                coefficient,
                label,
                provenance,
            )
        end for topology in topologies for mu in -1:1
    )
end
