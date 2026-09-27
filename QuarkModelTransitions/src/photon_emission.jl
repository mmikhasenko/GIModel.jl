# High-level photon-emission operator for the common matrix_element API.

"""
A spectroscopically selected photon-transition kernel.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = PhotonEmission(masses)
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.PhotonTransitionClass
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
abstract type PhotonTransitionClass end

"""
Leading M1 transition between states with the same radial quantum number.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
initial_wave = OscillatorWave(0, 0.5, [1.0])
initial = PhysicalState("parent", 3.5, [(
    basis=BasisState(1, "S", 3, 1; flavors=(:c, :c)),
    coefficient=1.0, wave=initial_wave,
)])
final = PhysicalState("daughter", 3.0, [(
    basis=BasisState(1, "S", 1, 0; flavors=(:c, :c)),
    coefficient=1.0, wave=OscillatorWave(0, 0.5, [1.0]),
)])
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.DirectM1
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
struct DirectM1 <: PhotonTransitionClass end

"""
M1 transition between states with different radial quantum numbers.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
initial_wave = OscillatorWave(0, 0.5, [0.0, 1.0])
initial = PhysicalState("parent", 3.5, [(
    basis=BasisState(2, "S", 3, 1; flavors=(:c, :c)),
    coefficient=1.0, wave=initial_wave,
)])
final = PhysicalState("daughter", 3.0, [(
    basis=BasisState(1, "S", 1, 0; flavors=(:c, :c)),
    coefficient=1.0, wave=OscillatorWave(0, 0.5, [1.0]),
)])
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.HinderedM1
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
struct HinderedM1 <: PhotonTransitionClass end

"""
Spin-conserving electric-dipole transition between S and P states.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
initial_wave = OscillatorWave(1, 0.5, [1.0])
initial = PhysicalState("parent", 3.5, [(
    basis=BasisState(1, "P", 3, 1; flavors=(:c, :c)),
    coefficient=1.0, wave=initial_wave,
)])
final = PhysicalState("daughter", 3.0, [(
    basis=BasisState(1, "S", 3, 1; flavors=(:c, :c)),
    coefficient=1.0, wave=OscillatorWave(0, 0.5, [1.0]),
)])
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.AllowedE1
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
struct AllowedE1 <: PhotonTransitionClass end

"""
Spin-flip electric-dipole transition from a triplet P1 to a singlet S0 state.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
initial_wave = OscillatorWave(1, 0.5, [1.0])
initial = PhysicalState("parent", 3.5, [(
    basis=BasisState(1, "P", 3, 1; flavors=(:c, :c)),
    coefficient=1.0, wave=initial_wave,
)])
final = PhysicalState("daughter", 3.0, [(
    basis=BasisState(1, "S", 1, 0; flavors=(:c, :c)),
    coefficient=1.0, wave=OscillatorWave(0, 0.5, [1.0]),
)])
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.SpinFlipE1
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
struct SpinFlipE1 <: PhotonTransitionClass end

"""
Spin-flip magnetic-quadrupole transition from a triplet P2 to a singlet S0 state.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
initial_wave = OscillatorWave(1, 0.5, [1.0])
initial = PhysicalState("parent", 3.5, [(
    basis=BasisState(1, "P", 3, 2; flavors=(:c, :c)),
    coefficient=1.0, wave=initial_wave,
)])
final = PhysicalState("daughter", 3.0, [(
    basis=BasisState(1, "S", 1, 0; flavors=(:c, :c)),
    coefficient=1.0, wave=OscillatorWave(0, 0.5, [1.0]),
)])
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.SpinFlipM2
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
"""
struct SpinFlipM2 <: PhotonTransitionClass end

_photon_multipole(::Union{DirectM1,HinderedM1}) = :M1
_photon_multipole(::Union{AllowedE1,SpinFlipE1}) = :E1
_photon_multipole(::SpinFlipM2) = :M2

abstract type PhotonCurrent end

"""Derive the electromagnetic current from explicit quark flavors."""
struct StandardPhotonCurrent <: PhotonCurrent end

"""
PhotonEmitter(flavors, constituent, coefficient)

One resolved term in an electromagnetic current. This technical adapter is
used by paper-reproduction code whose legacy `:q` states do not retain enough
flavor information to derive isospin charges automatically.

## Example

```julia
using QuarkModelTransitions
import QuarkModelTransitions as QMT
emitter = QMT.PhotonEmitter((:c, :c), 1, 4/3)
@assert emitter.coefficient == 4/3
```

## Related

[`PhotonEmission`](@ref).
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

"""A pre-resolved current for reproducing a source with coarse flavor states."""
struct ResolvedPhotonCurrent{E<:Tuple} <: PhotonCurrent
    emitters::E
    function ResolvedPhotonCurrent(emitters)
        resolved = Tuple(emitters)
        isempty(resolved) && throw(ArgumentError("resolved photon current cannot be empty"))
        all(term -> term isa PhotonEmitter, resolved) || throw(ArgumentError(
            "every resolved photon-current term must be a PhotonEmitter",
        ))
        return new{typeof(resolved)}(resolved)
    end
end

"""
    PhotonEmission(quark_masses; recoil_order=0, ...)
    PhotonEmission(multipole, quark_masses, emitters; recoil=false, ...)

The Godfrey--Isgur photon-emission operator in the Appendix-D mock-meson
prescription. The transition class and multipole are inferred from the final
and initial states by `photon_transition_class`.

`recoil_order=0` selects the class-specific leading expression. Order `2`
adds the relative `(qr)^2` correction where that correction is implemented;
currently this is the published M1 `E_2` term. Unsupported combinations fail
explicitly rather than silently dropping a requested correction.

The explicit-current constructor accepts `:M1`, `:E1`, or `:M2` and
`PhotonEmitter` terms; the supplied multipole must match the states.
Use [`mass_correction_factor`](@ref) for comparisons at another photon momentum.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = PhotonEmission(masses)
@assert operator.recoil_order == 0
```

## Related

- `AllowedE1` — spin-conserving S–P classification.
- `DirectM1` — same-radial-level M1 classification.
- `ELECTROMAGNETIC_DEFAULTS` — phenomenological electromagnetic inputs.
- `HinderedM1` — different-radial-level M1 classification.
- `PhotonTransitionClass` — common photon classification type.
- `RadiativeAmplitude` — photon result and coherent terms.
- `SpinFlipE1` — triplet P1 to singlet S0 classification.
- `SpinFlipM2` — triplet P2 to singlet S0 classification.
- `e1_transition_amplitude` — assemble an E1 amplitude.
- `m1_transition_moment` — assemble a moment in nuclear magnetons.
- [`matrix_element`](@ref) — evaluate an operator between states.
- `photon_momentum` — two-body photon momentum in GeV.
- `photon_transition_class` — infer the kernel from spectroscopy.
- `spin_flip_photon_amplitude` — assemble a spin-flip photon amplitude.
"""
struct PhotonEmission{C<:PhotonCurrent} <: TransitionOperator
    quark_masses::QuarkMassTable
    recoil_order::Int
    current::C
    recoil_form_factor::Bool
    magnetic_exponent::Float64
    electric_exponent::Float64
    recoil_beta_GeV::Float64
end

function PhotonEmission(
    quark_masses::QuarkMassTable;
    recoil_order::Integer = 0,
    current::PhotonCurrent = StandardPhotonCurrent(),
    recoil_form_factor::Bool = false,
    magnetic_exponent::Real = ELECTROMAGNETIC_DEFAULTS.magnetic_exponent,
    electric_exponent::Real = ELECTROMAGNETIC_DEFAULTS.electric_exponent,
    recoil_beta_GeV::Real = ELECTROMAGNETIC_DEFAULTS.recoil_beta_GeV,
)
    recoil_order in (0, 2) || throw(ArgumentError(
        "recoil_order must be 0 (leading expression) or 2 (relative (qr)^2 correction)",
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
    return PhotonEmission{typeof(current)}(
        copy(quark_masses), Int(recoil_order), current, recoil_form_factor,
        fm, fe, beta,
    )
end

"""
A photon matrix element and its complete coherent component decomposition.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = PhotonEmission(masses)
amplitude = matrix_element(final, operator, initial)
@assert amplitude isa QuarkModelTransitions.RadiativeAmplitude
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
- [`decay_width`](@ref) — convert a transition result to MeV.
- `photon_transition_class` — infer the kernel from spectroscopy.
"""
struct RadiativeAmplitude{O,I,F,K,C<:PhotonTransitionClass,V<:Number,T<:Tuple,P}
    operator::O
    initial::I
    final::F
    kinematics::K
    transition_class::C
    multipole::Symbol
    value::V
    terms::T
    provenance::P
end

function _photon_mass(operator::PhotonEmission, flavor::Symbol)
    canonical = flavor in (:u, :d, :n) ? :q : flavor
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

function _component_transition_class(daughter::StateComponent, parent::StateComponent)
    Ld = orbital_angular_momentum(daughter.basis.L_label)
    Lp = orbital_angular_momentum(parent.basis.L_label)
    md, mp = daughter.basis.multiplicity, parent.basis.multiplicity
    if Ld == 0 && Lp == 0 && Set((md, mp)) == Set((1, 3))
        return daughter.basis.n == parent.basis.n ? DirectM1() : HinderedM1()
    end
    Set((Ld, Lp)) == Set((0, 1)) || throw(ArgumentError(
        "photon kernel supports S--S spin flips and S--P transitions; got L=$Lp -> L=$Ld",
    ))
    sw, pw = Ld == 0 ? (daughter, parent) : (parent, daughter)
    if sw.basis.multiplicity == pw.basis.multiplicity
        return AllowedE1()
    end
    (sw.basis.multiplicity, pw.basis.multiplicity) == (1, 3) || throw(ArgumentError(
        "spin-flip photon kernel requires 3P_J -> 1S0",
    ))
    Lp == 1 || throw(ArgumentError(
        "the implemented spin-flip photon kernel requires the P wave to be the parent",
    ))
    pw.basis.J == 1 && return SpinFlipE1()
    pw.basis.J == 2 && return SpinFlipM2()
    throw(ArgumentError("3P$(pw.basis.J) -> 1S0 has no implemented photon kernel"))
end

_same_photon_family(a::PhotonTransitionClass, b::PhotonTransitionClass) =
    (a isa Union{DirectM1,HinderedM1} && b isa Union{DirectM1,HinderedM1}) ||
    typeof(a) == typeof(b)

"""
    photon_transition_class(final, initial)

Infer the photon kernel from the states' orbital angular momentum, spin,
total angular momentum, and radial quantum numbers. All coherently mixed
components must select one compatible class.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = PhotonEmission(masses)
@assert QuarkModelTransitions.photon_transition_class(final, initial) isa QMT.DirectM1
```

## Related

- [`PhotonEmission`](@ref) — infer and evaluate a photon transition.
- `RadiativeAmplitude` — inspect the selected class on a result.
"""
function photon_transition_class(final::PhysicalState, initial::PhysicalState)
    classes = PhotonTransitionClass[]
    for parent in initial.components, daughter in final.components
        daughter.basis.flavors == parent.basis.flavors || continue
        push!(classes, _component_transition_class(daughter, parent))
    end
    isempty(classes) && throw(ArgumentError(
        "initial and final states have no shared pure-flavor component",
    ))
    selected = first(classes)
    all(candidate -> _same_photon_family(selected, candidate), classes) ||
        throw(ArgumentError("mixed-state components imply incompatible photon-transition classes"))
    if selected isa Union{DirectM1,HinderedM1}
        return any(candidate -> candidate isa HinderedM1, classes) ? HinderedM1() : DirectM1()
    end
    return selected
end

function _quark_charge(flavor::Symbol)
    flavor in (:u, :c) && return 2 / 3
    flavor in (:d, :s, :b) && return -1 / 3
    flavor in (:q, :n) && throw(ArgumentError(
        "flavor :$flavor does not determine an electromagnetic current; " *
        "use explicit :u/:d flavor components",
    ))
    throw(ArgumentError("unsupported electromagnetic flavor :$flavor"))
end

function _photon_emitters(
    ::StandardPhotonCurrent,
    class::PhotonTransitionClass,
    flavors::Tuple{Symbol,Symbol},
)
    f1, f2 = flavors
    c1, c2 = _quark_charge(f1), _quark_charge(f2)
    class isa Union{SpinFlipE1,SpinFlipM2} && (c2 = -c2)
    return (PhotonEmitter(flavors, 1, c1), PhotonEmitter(flavors, 2, c2))
end

function _photon_emitters(
    current::ResolvedPhotonCurrent,
    ::PhotonTransitionClass,
    flavors::Tuple{Symbol,Symbol},
)
    return filter(emitter -> emitter.flavors == flavors, current.emitters)
end

function _photon_pure_component(
    operator::PhotonEmission,
    class::Union{DirectM1,HinderedM1},
    emitter::PhotonEmitter,
    daughter::StateComponent,
    parent::StateComponent,
    q::Real,
)
    masses = (_photon_mass(operator, parent.basis.flavors[1]),
              _photon_mass(operator, parent.basis.flavors[2]))
    singlet, triplet = parent.basis.multiplicity == 1 ?
        (parent, daughter) : (daughter, parent)
    if operator.recoil_order == 2
        masses[1] == masses[2] || throw(ArgumentError(
            "the published order-2 M1 recoil prescription is equal-mass only",
        ))
        return m1_recoil_moment(
            singlet.wave, triplet.wave, masses[emitter.constituent],
            emitter.coefficient, q;
            magnetic_exponent = operator.magnetic_exponent,
            electric_exponent = operator.electric_exponent,
        )
    end
    return m1_transition_moment(
        momentum_wave(singlet.wave, 0), momentum_wave(triplet.wave, 0),
        masses[1], masses[2],
        [(emitter.coefficient, masses[emitter.constituent])];
        exponent = operator.magnetic_exponent,
    )
end

function _photon_pure_component(
    operator::PhotonEmission,
    class::AllowedE1,
    emitter::PhotonEmitter,
    daughter::StateComponent,
    parent::StateComponent,
    q::Real,
)
    operator.recoil_order == 0 || throw(ArgumentError(
        "recoil_order=2 is not implemented for AllowedE1",
    ))
    Ld = orbital_angular_momentum(daughter.basis.L_label)
    sw, pw = Ld == 0 ? (daughter, parent) : (parent, daughter)
    angular = e1_angular_coefficient(
        pw.basis.J;
        singlet = sw.basis.multiplicity == 1,
        parent_is_S = orbital_angular_momentum(parent.basis.L_label) == 0,
    )
    m_emit = _photon_mass(operator, parent.basis.flavors[emitter.constituent])
    return e1_transition_amplitude(
        sw.wave, momentum_wave(sw.wave, 0), pw.wave, momentum_wave(pw.wave, 1),
        m_emit, qvalue -> emitter.coefficient * angular * qvalue,
        1.0, 0.0; q = q, exponent = operator.electric_exponent,
    )
end

function _photon_pure_component(
    operator::PhotonEmission,
    class::Union{SpinFlipE1,SpinFlipM2},
    emitter::PhotonEmitter,
    daughter::StateComponent,
    parent::StateComponent,
    q::Real,
)
    operator.recoil_order == 0 || throw(ArgumentError(
        "recoil_order=2 is not implemented for $(nameof(typeof(class)))",
    ))
    Ld = orbital_angular_momentum(daughter.basis.L_label)
    sw, pw = Ld == 0 ? (daughter, parent) : (parent, daughter)
    m_emit = _photon_mass(operator, parent.basis.flavors[emitter.constituent])
    return spin_flip_photon_amplitude(
        sw.wave, pw.wave, [(emitter.coefficient, m_emit)], pw.basis.J, q;
        exponent = operator.electric_exponent,
    )
end

function _assemble_photon_amplitude(
    class::PhotonTransitionClass,
    final::PhysicalState,
    operator::PhotonEmission,
    initial::PhysicalState,
    resolved::CMKinematics,
)
    q = resolved.momentum_GeV
    terms = TransitionTerm[]
    total = 0.0 + 0.0im
    matched = false
    for parent in initial.components, daughter in final.components
        flavors = parent.basis.flavors
        daughter.basis.flavors == flavors || continue
        component_class = _component_transition_class(daughter, parent)
        _same_photon_family(class, component_class) || continue
        mixing = parent.coefficient * conj(daughter.coefficient)
        for emitter in _photon_emitters(operator.current, class, flavors)
            matched = true
            pure = _photon_pure_component(operator, component_class, emitter, daughter, parent, q)
            operator.recoil_form_factor &&
                (pure *= photon_recoil_form_factor(q; beta = operator.recoil_beta_GeV))
            contribution = mixing * pure
            total += contribution
            T = ComplexF64
            push!(terms, TransitionTerm(
                "$(parent.basis.label) -> $(daughter.basis.label), emitter=$(emitter.constituent)",
                T(mixing), one(T), T(pure), T(contribution),
                (
                    source = :GI1985_AppendixD,
                    transition_class = nameof(typeof(component_class)),
                    multipole = _photon_multipole(component_class),
                    recoil_order = operator.recoil_order,
                    flavors = emitter.flavors,
                    constituent = emitter.constituent,
                    charge_coefficient = emitter.coefficient,
                ),
            ))
        end
    end
    matched || throw(ArgumentError(
        "the photon current has no term matching a shared pure-flavor component",
    ))
    multipole = _photon_multipole(class)
    units = multipole == :M1 ? :nuclear_magnetons : :MeV_sqrt
    provenance = (
        backend = :mock_meson_appendix_d,
        source = (:GI1985_Eq22, :GI1985_AppendixD),
        prescription = :hybrid_mock_meson,
        transition_class = nameof(typeof(class)),
        recoil_order = operator.recoil_order,
        magnetic_exponent = operator.magnetic_exponent,
        electric_exponent = operator.electric_exponent,
        recoil_form_factor = operator.recoil_form_factor,
        recoil_beta_GeV = operator.recoil_beta_GeV,
        amplitude_units = units,
        width_units = :MeV,
    )
    return RadiativeAmplitude(
        operator, initial, final, resolved, class, multipole,
        total, Tuple(terms), provenance,
    )
end

photon_emission_matrix_element(
    class::Union{DirectM1,HinderedM1}, final, operator, initial, resolved,
) = _assemble_photon_amplitude(class, final, operator, initial, resolved)

photon_emission_matrix_element(
    class::AllowedE1, final, operator, initial, resolved,
) = _assemble_photon_amplitude(class, final, operator, initial, resolved)

photon_emission_matrix_element(
    class::Union{SpinFlipE1,SpinFlipM2}, final, operator, initial, resolved,
) = _assemble_photon_amplitude(class, final, operator, initial, resolved)


function matrix_element(
    final::PhysicalState,
    operator::PhotonEmission,
    initial::PhysicalState;
    verbose::Bool = false,
)
    return _photon_matrix_element(final, operator, initial,
        _on_shell_photon_momentum(final, initial); verbose)
end

function mass_correction_factor(final::PhysicalState, operator::PhotonEmission,
                                initial::PhysicalState; target_momentum::Real)
    reference = matrix_element(final, operator, initial)
    target = _photon_matrix_element(final, operator, initial, target_momentum)
    return _correction_ratio(target.value, reference.value)
end

function _photon_matrix_element(final::PhysicalState, operator::PhotonEmission,
                                initial::PhysicalState, momentum::Real; verbose::Bool=false)
    resolved = CMKinematics(momentum)
    q = resolved.momentum_GeV
    q isa Real && q >= 0 || throw(ArgumentError(
        "PhotonEmission requires a real non-negative photon momentum",
    ))
    class = photon_transition_class(final, initial)
    if verbose
        @info "photon transition selected" initial=initial.label final=final.label transition_class=nameof(typeof(class)) multipole=_photon_multipole(class) recoil_order=operator.recoil_order q_GeV=q
    end
    return photon_emission_matrix_element(class, final, operator, initial, resolved)
end

function matrix_element(
    final::TransitionState,
    ::PhotonEmission,
    initial::TransitionState;
    verbose::Bool = false,
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
    return decay_width(matrix_element(final, operator, initial))
end

function Base.show(io::IO, amplitude::RadiativeAmplitude)
    print(io, "RadiativeAmplitude(", amplitude.initial.label, " -> ",
        amplitude.final.label, " + gamma, ", nameof(typeof(amplitude.transition_class)),
        ", recoil_order=", amplitude.operator.recoil_order, ", ", amplitude.value, ")")
end

# Explicit-current compatibility for callers of the independently developed API.
struct SpecifiedPhotonCurrent{C<:PhotonCurrent} <: PhotonCurrent
    multipole::Symbol
    resolved::C
end
function PhotonEmission(multipole::Symbol, masses::QuarkMassTable, emitters;
                        recoil::Bool=false, kwargs...)
    multipole in (:M1, :E1, :M2) || throw(ArgumentError("unsupported photon multipole"))
    recoil && multipole != :M1 && throw(ArgumentError("recoil is implemented only for M1"))
    terms = emitters isa PhotonEmitter ? (emitters,) : Tuple(emitters)
    current = SpecifiedPhotonCurrent(multipole, ResolvedPhotonCurrent(terms))
    return PhotonEmission(masses; current, recoil_order=recoil ? 2 : 0, kwargs...)
end
function _photon_emitters(current::SpecifiedPhotonCurrent, class::PhotonTransitionClass,
                         flavors::Tuple{Symbol,Symbol})
    _photon_multipole(class) == current.multipole ||
        throw(ArgumentError("specified multipole is incompatible with the states"))
    return _photon_emitters(current.resolved, class, flavors)
end
