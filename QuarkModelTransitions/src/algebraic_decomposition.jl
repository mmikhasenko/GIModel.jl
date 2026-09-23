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

# =============================================================================
# Spectroscopic recoupling and one shared integral basis
# =============================================================================

"""A spatial label after resolving the external `L dot S` coupled states."""
struct _SpectroscopicOrbitalLabel{L}
    elementary::L
    initial_L::Int
    initial_mL::Int
    final_L::Int
    final_mL::Int
end

"""One recoupled contribution before equal spatial labels are collected."""
struct _SpectroscopicCoefficientTerm{B,L,C}
    elementary::B
    orbital_label::L
    initial_cg::Float64
    final_cg::Float64
    coefficient::C
end

function _basis_spin(basis::BasisState)
    basis.multiplicity in (1, 3) || throw(ArgumentError(
        "elementary meson emission currently supports singlet/triplet q-qbar states",
    ))
    return (basis.multiplicity - 1) ÷ 2
end

function _spectroscopic_decomposition(
    piece::_ElementaryEmissionPiece,
    daughter::BasisState,
    parent::BasisState,
    daughter_flavor::_FlavorState,
    emitted_flavor::_FlavorState,
    parent_flavor::_FlavorState,
    helicity::Integer,
)
    Mi = Mf = Int(helicity)
    abs(Mi) <= parent.J || throw(ArgumentError("parent helicity is outside its J"))
    abs(Mf) <= daughter.J || throw(ArgumentError("daughter helicity is outside its J"))
    Li = orbital_angular_momentum(parent.L_label)
    Lf = orbital_angular_momentum(daughter.L_label)
    Si = _basis_spin(parent)
    Sf = _basis_spin(daughter)
    terms = _SpectroscopicCoefficientTerm[]
    for mLi in -Li:Li, mSi in -Si:Si
        mLi + mSi == Mi || continue
        initial_cg = Float64(CG(Li, mLi, Si, mSi, parent.J, Mi))
        iszero(initial_cg) && continue
        initial_spin = _coupled_spin(Si, mSi)
        for mLf in -Lf:Lf, mSf in -Sf:Sf
            mLf + mSf == Mf || continue
            final_cg = Float64(CG(Lf, mLf, Sf, mSf, daughter.J, Mf))
            iszero(final_cg) && continue
            final_spin = _coupled_spin(Sf, mSf)
            elementary = _eq19_coefficient_decomposition(
                piece,
                final_spin,
                initial_spin,
                daughter_flavor,
                emitted_flavor,
                parent_flavor,
            )
            for term in elementary
                iszero(term.coefficient) && continue
                label = _SpectroscopicOrbitalLabel(
                    term.orbital_label, Li, mLi, Lf, mLf,
                )
                coefficient = conj(final_cg) * initial_cg * term.coefficient
                push!(terms, _SpectroscopicCoefficientTerm(
                    term, label, initial_cg, final_cg, coefficient,
                ))
            end
        end
    end
    return Tuple(terms)
end

"""
Symbolic Eq. (19) helicities and partial waves over one common integral basis.

Columns name spatial integrals. `helicity_coefficients` converts those integrals
to the compressed Appendix-C helicity vector, and `partial_wave_coefficients`
converts the same columns directly to every allowed partial wave.
"""
struct _AngularCoefficientDecomposition{H,W,L,M,P,R}
    helicities::H
    partial_waves::W
    integral_labels::L
    helicity_coefficients::M
    partial_wave_coefficients::P
    provenance::R
end

function _find_or_append_label!(labels, label)
    index = findfirst(existing -> isequal(existing, label), labels)
    isnothing(index) || return index
    push!(labels, label)
    return length(labels)
end

function _eq19_angular_decomposition(
    daughter::BasisState,
    parent::BasisState,
    daughter_flavor::_FlavorState,
    emitted_flavor::_FlavorState,
    parent_flavor::_FlavorState,
)
    daughter_state = ReferenceState(
        "daughter", 1.0; J = daughter.J,
        parity = iseven(orbital_angular_momentum(daughter.L_label)) ? -1 : 1,
    )
    emitted_state = ReferenceState("emitted-P", 0.1; J = 0, parity = -1)
    parent_state = ReferenceState(
        "parent", 2.0; J = parent.J,
        parity = iseven(orbital_angular_momentum(parent.L_label)) ? -1 : 1,
    )
    channel = TwoMesonChannel(daughter_state, emitted_state)
    projection = partial_wave_projection(channel, parent_state)
    labels = _SpectroscopicOrbitalLabel[]
    rows = Vector{Vector{Tuple{Int,ComplexF64}}}()
    for helicity_label in projection.helicities
        twice_m = channel.first.label == "daughter" ?
            helicity_label.twice_lambda1 : helicity_label.twice_lambda2
        iseven(twice_m) || throw(ArgumentError("meson helicity must be integral"))
        m = twice_m ÷ 2
        entries = Tuple{Int,ComplexF64}[]
        scale = iszero(m) ? 1.0 : sqrt(2.0)
        for piece in (_DirectPseudoscalarPiece(), _RecoilPseudoscalarPiece())
            for term in _spectroscopic_decomposition(
                piece,
                daughter,
                parent,
                daughter_flavor,
                emitted_flavor,
                parent_flavor,
                m,
            )
                column = _find_or_append_label!(labels, term.orbital_label)
                push!(entries, (column, scale * ComplexF64(term.coefficient)))
            end
        end
        push!(rows, entries)
    end
    helicity_coefficients = zeros(ComplexF64, length(rows), length(labels))
    for (row, entries) in pairs(rows), (column, coefficient) in entries
        helicity_coefficients[row, column] += coefficient
    end
    partial_wave_coefficients = projection.matrix * helicity_coefficients
    provenance = (
        source = (:GI1985_Eq19, :GI1985_AppendixB, :GI1985_AppendixC),
        coupling_order = :L_dot_S,
        helicity_compression = :C3,
        numerical_spatial_integrals = false,
    )
    return _AngularCoefficientDecomposition(
        projection.helicities,
        projection.partial_waves,
        Tuple(labels),
        helicity_coefficients,
        partial_wave_coefficients,
        provenance,
    )
end
