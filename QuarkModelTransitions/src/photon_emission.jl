# High-level photon-emission operator for the common matrix_element API.

"""
    PhotonEmitter(flavors, constituent, coefficient)

One explicitly resolved term in the electromagnetic current. `flavors` is the
ordered `(quark, antiquark)` pair to which the term applies, `constituent` is
`1` or `2`, and `coefficient` contains the charge and flavor/isospin factor.

Keeping this information explicit is intentional: a coarse `(:q, :q)` radial
component does not say whether a physical transition contains `e_u + e_d` or
`e_u - e_d`.
"""
struct PhotonEmitter
    flavors::Tuple{Symbol,Symbol}
    constituent::Int
    coefficient::Float64
    function PhotonEmitter(
        flavors::Tuple{Symbol,Symbol},
        constituent::Integer,
        coefficient::Real,
    )
        constituent in (1, 2) || throw(ArgumentError(
            "photon-emitting constituent must be 1 (quark) or 2 (antiquark)",
        ))
        c = Float64(coefficient)
        isfinite(c) || throw(ArgumentError("photon-current coefficient must be finite"))
        return new(flavors, Int(constituent), c)
    end
end

"""
    PhotonEmission(multipole, quark_masses, emitters;
                   recoil=false, recoil_form_factor=false,
                   magnetic_exponent=0.7, electric_exponent=0.5,
                   recoil_beta_GeV=0.40)

The Godfrey--Isgur photon-emission operator in the published Appendix-D
mock-meson prescription. `multipole` is `:M1`, `:E1`, or `:M2`; `emitters` is
an iterable of [`PhotonEmitter`](@ref) terms.

This is the paper's hybrid realization of the one-body current introduced by
Eq. (22), not a literal unsmeared evaluation of Eq. (22). The two `m/E`
exponents and optional recoil form factor are stored in the operator and
therefore retained in result provenance.
"""
struct PhotonEmission{E<:Tuple} <: TransitionOperator
    multipole::Symbol
    quark_masses::QuarkMassTable
    emitters::E
    recoil::Bool
    recoil_form_factor::Bool
    magnetic_exponent::Float64
    electric_exponent::Float64
    recoil_beta_GeV::Float64
end

function PhotonEmission(
    multipole::Symbol,
    quark_masses::QuarkMassTable,
    emitters;
    recoil::Bool = false,
    recoil_form_factor::Bool = false,
    magnetic_exponent::Real = ELECTROMAGNETIC_DEFAULTS.magnetic_exponent,
    electric_exponent::Real = ELECTROMAGNETIC_DEFAULTS.electric_exponent,
    recoil_beta_GeV::Real = ELECTROMAGNETIC_DEFAULTS.recoil_beta_GeV,
)
    multipole in (:M1, :E1, :M2) || throw(ArgumentError(
        "photon multipole must be :M1, :E1, or :M2",
    ))
    resolved = Tuple(emitters)
    isempty(resolved) && throw(ArgumentError("PhotonEmission needs at least one emitter term"))
    all(term -> term isa PhotonEmitter, resolved) || throw(ArgumentError(
        "every photon-current term must be a PhotonEmitter",
    ))
    all(isfinite(mass) && mass > 0 for mass in values(quark_masses)) ||
        throw(ArgumentError("constituent masses must be finite and positive"))
    fm, fe, beta = Float64(magnetic_exponent), Float64(electric_exponent),
                   Float64(recoil_beta_GeV)
    isfinite(fm) && fm >= 0 || throw(ArgumentError(
        "magnetic m/E exponent must be finite and non-negative",
    ))
    isfinite(fe) && fe >= 0 || throw(ArgumentError(
        "electric m/E exponent must be finite and non-negative",
    ))
    isfinite(beta) && beta > 0 || throw(ArgumentError(
        "photon recoil beta must be finite and positive",
    ))
    recoil && multipole != :M1 && throw(ArgumentError(
        "the implemented hindered-transition recoil term applies only to M1",
    ))
    return PhotonEmission{typeof(resolved)}(
        multipole, copy(quark_masses), resolved, recoil, recoil_form_factor,
        fm, fe, beta,
    )
end

PhotonEmission(
    multipole::Symbol,
    quark_masses::QuarkMassTable,
    emitter::PhotonEmitter;
    kwargs...,
) = PhotonEmission(multipole, quark_masses, (emitter,); kwargs...)

"""A photon matrix element and its complete coherent component decomposition."""
struct RadiativeAmplitude{O,I,F,K,V<:Number,T<:Tuple,P}
    operator::O
    initial::I
    final::F
    kinematics::K
    multipole::Symbol
    value::V
    terms::T
    provenance::P
end

function _photon_mass(operator::PhotonEmission, flavor::Symbol)
    canonical = flavor === :n ? :q : flavor
    key = String(canonical)
    haskey(operator.quark_masses, key) || throw(ArgumentError(
        "PhotonEmission has no constituent mass for flavor :$flavor",
    ))
    return operator.quark_masses[key]
end

function _on_shell_photon_momentum(final::TransitionState, initial::TransitionState)
    initial.mass_GeV > final.mass_GeV || throw(ArgumentError(
        "on-shell photon channel $(initial.label) -> $(final.label) + gamma is closed: " *
        "parent mass $(initial.mass_GeV) GeV is not above daughter mass $(final.mass_GeV) GeV",
    ))
    return photon_momentum(initial.mass_GeV, final.mass_GeV)
end

_resolve_kinematics(final::TransitionState, initial::TransitionState, ::OnShell) =
    CMKinematics(_on_shell_photon_momentum(final, initial))

function _photon_component_kind(
    multipole::Symbol,
    daughter::StateComponent,
    parent::StateComponent,
)
    Ld = orbital_angular_momentum(daughter.basis.L_label)
    Lp = orbital_angular_momentum(parent.basis.L_label)
    if multipole == :M1
        Ld == Lp || throw(ArgumentError("M1 requires equal internal orbital angular momentum"))
        daughter.basis.multiplicity != parent.basis.multiplicity || throw(ArgumentError(
            "the implemented M1 kernel requires a singlet--triplet transition",
        ))
        Set((daughter.basis.multiplicity, parent.basis.multiplicity)) == Set((1, 3)) ||
            throw(ArgumentError("M1 requires spin multiplicities 1 and 3"))
        Ld == 0 || throw(ArgumentError("the implemented M1 mock-meson kernel is S-wave only"))
        return :magnetic
    end
    Set((Ld, Lp)) == Set((0, 1)) || throw(ArgumentError(
        "$multipole requires an S--P transition in the implemented Table-VI kernel",
    ))
    sw = Ld == 0 ? daughter : parent
    pw = Ld == 1 ? daughter : parent
    if sw.basis.multiplicity == pw.basis.multiplicity
        multipole == :E1 || throw(ArgumentError("spin-conserving S--P emission is E1"))
        return :electric
    end
    (sw.basis.multiplicity, pw.basis.multiplicity) == (1, 3) || throw(ArgumentError(
        "the implemented spin-flip kernel requires 1S0 and 3P_J states",
    ))
    expected = pw.basis.J == 2 ? :M2 : pw.basis.J == 1 ? :E1 : nothing
    multipole == expected || throw(ArgumentError(
        "3P$(pw.basis.J) -> 1S0 spin flip requires $(something(expected, :unsupported))",
    ))
    Lp == 1 || throw(ArgumentError(
        "the implemented spin-flip Table-VI kernel is for P-wave parent emission",
    ))
    return :spin_flip
end

function _photon_pure_component(
    operator::PhotonEmission,
    emitter::PhotonEmitter,
    daughter::StateComponent,
    parent::StateComponent,
    q::Real,
)
    flavors = parent.basis.flavors
    daughter.basis.flavors == flavors || return nothing
    emitter.flavors == flavors || return nothing
    kind = _photon_component_kind(operator.multipole, daughter, parent)
    masses = (_photon_mass(operator, flavors[1]), _photon_mass(operator, flavors[2]))
    m_emit = masses[emitter.constituent]
    coefficient = emitter.coefficient
    value = if kind == :magnetic
        singlet, triplet = parent.basis.multiplicity == 1 ?
            (parent, daughter) : (daughter, parent)
        if operator.recoil
            masses[1] == masses[2] || throw(ArgumentError(
                "the published hindered-M1 recoil prescription is equal-flavor only",
            ))
            m1_recoil_moment(
                singlet.wave, triplet.wave, m_emit, coefficient, q;
                magnetic_exponent = operator.magnetic_exponent,
                electric_exponent = operator.electric_exponent,
            )
        else
            m1_transition_moment(
                momentum_wave(singlet.wave, 0), momentum_wave(triplet.wave, 0),
                masses[1], masses[2], [(coefficient, m_emit)];
                exponent = operator.magnetic_exponent,
            )
        end
    else
        Ld = orbital_angular_momentum(daughter.basis.L_label)
        sw, pw = Ld == 0 ? (daughter, parent) : (parent, daughter)
        if kind == :spin_flip
            spin_flip_photon_amplitude(
                sw.wave, pw.wave, [(coefficient, m_emit)], pw.basis.J, q;
                exponent = operator.electric_exponent,
            )
        else
            angular = e1_angular_coefficient(
                pw.basis.J;
                singlet = sw.basis.multiplicity == 1,
                parent_is_S = orbital_angular_momentum(parent.basis.L_label) == 0,
            )
            e1_transition_amplitude(
                sw.wave, momentum_wave(sw.wave, 0),
                pw.wave, momentum_wave(pw.wave, 1), m_emit,
                qvalue -> coefficient * angular * qvalue,
                1.0, 0.0; q = q, exponent = operator.electric_exponent,
            )
        end
    end
    operator.recoil_form_factor &&
        (value *= photon_recoil_form_factor(q; beta = operator.recoil_beta_GeV))
    return value
end

"""
    matrix_element(final, operator::PhotonEmission, initial; kinematics=OnShell())

Evaluate a radiative matrix element between resolved physical mesons. The
result coherently composes their pure flavor components and records each
emitter contribution. M1 values are magnetic moments in nuclear magnetons;
E1/M2 values are Table-VI amplitudes in `MeV^(1/2)`.
"""
function matrix_element(
    final::PhysicalState,
    operator::PhotonEmission,
    initial::PhysicalState;
    kinematics::TransitionKinematics = OnShell(),
)
    resolved = _resolve_kinematics(final, initial, kinematics)
    q = resolved.momentum_GeV
    q isa Real && q >= 0 || throw(ArgumentError(
        "PhotonEmission requires a real non-negative photon momentum",
    ))
    terms = TransitionTerm[]
    total = 0.0 + 0.0im
    matched = false
    for parent in initial.components, daughter in final.components
        mixing = parent.coefficient * conj(daughter.coefficient)
        for emitter in operator.emitters
            pure = _photon_pure_component(operator, emitter, daughter, parent, q)
            isnothing(pure) && continue
            matched = true
            contribution = mixing * pure
            total += contribution
            T = ComplexF64
            push!(terms, TransitionTerm(
                "$(parent.basis.label) -> $(daughter.basis.label), " *
                "emitter=$(emitter.constituent)",
                T(mixing), one(T), T(pure), T(contribution),
                (
                    source = :GI1985_AppendixD,
                    multipole = operator.multipole,
                    flavors = emitter.flavors,
                    constituent = emitter.constituent,
                    charge_coefficient = emitter.coefficient,
                ),
            ))
        end
    end
    matched || throw(ArgumentError(
        "no PhotonEmitter term matches a shared pure-flavor component of " *
        "$(initial.label) and $(final.label)",
    ))
    units = operator.multipole == :M1 ? :nuclear_magnetons : :MeV_sqrt
    provenance = (
        backend = :mock_meson_appendix_d,
        source = (:GI1985_Eq22, :GI1985_AppendixD, :GI1985_TableVI),
        prescription = :hybrid_mock_meson,
        magnetic_exponent = operator.magnetic_exponent,
        electric_exponent = operator.electric_exponent,
        hindered_m1_recoil = operator.recoil,
        recoil_form_factor = operator.recoil_form_factor,
        recoil_beta_GeV = operator.recoil_beta_GeV,
        amplitude_units = units,
        width_units = :MeV,
    )
    return RadiativeAmplitude(
        operator, initial, final, resolved, operator.multipole,
        total, Tuple(terms), provenance,
    )
end

function matrix_element(
    final::TransitionState,
    ::PhotonEmission,
    initial::TransitionState;
    kinematics::TransitionKinematics = OnShell(),
)
    throw(ArgumentError(
        "PhotonEmission requires resolved PhysicalState inputs; got " *
        "$(nameof(typeof(initial))) -> $(nameof(typeof(final)))",
    ))
end

function decay_width(amplitude::RadiativeAmplitude)
    q = amplitude.kinematics.momentum_GeV
    if amplitude.multipole == :M1
        return 1000 * m1_radiative_width(
            abs(amplitude.value), q; parent_spin = amplitude.initial.J,
        )
    end
    return abs2(amplitude.value)
end

function decay_width(
    final::PhysicalState,
    operator::PhotonEmission,
    initial::PhysicalState,
)
    initial.mass_GeV <= final.mass_GeV && return 0.0
    return decay_width(matrix_element(final, operator, initial; kinematics = OnShell()))
end

function Base.show(io::IO, amplitude::RadiativeAmplitude)
    print(io, "RadiativeAmplitude(", amplitude.initial.label, " -> ",
        amplitude.final.label, " + gamma, ", amplitude.multipole, ", ",
        amplitude.value, ")")
end
