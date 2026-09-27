# =============================================================================
# Strong decays  M* -> M + P   (Sec. IV, Tables IV-V, Appendices B-C)
# =============================================================================
#
# PROVENANCE LEGEND
#   [PAPER]   transcribed / coded directly from Godfrey-Isgur 1985. NOT derived
#             here: per-row flavor-spin coefficients (App. B), the Table IV class
#             algebra and its k-values, the form-factor shape, beta, and the two
#             numeric fit targets.
#   [DERIVED] computed analytically in this module: the two-body breakup momentum
#             (Kallen), the harmonic-oscillator momentum-space spatial overlap,
#             and the numeric solve that turns the two fit rows into (A, S0).
#
# A Table V row factorizes as
#
#   amp  =  c  *  X(qbar)  *  [ qbar^L * sqrt(q/2pi) * exp(-q^2 / 16 beta^2) ]
#           ^^^   ^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#     coefficient  reduced             spatial_overlap
#       c [PAPER]  X [Table IV]        [DERIVED]  (qbar = q/beta)
#
#   * c    : signed flavor-spin-color factor, one per decay row      [PAPER]
#   * X    : Table IV reduced partial-wave amplitude; its FORM is    [PAPER]
#            the paper's, but the numeric strengths A, S0 are fit    [DERIVED]
#            here from two rows.
#   * spatial_overlap : the single-beta SHO momentum overlap.        [DERIVED]
#
# The result is in MeV^(1/2); the partial width is its square (`decay_width`).
# `StrongDecayAmplitude` returns all three factors so a computed number is
# self-documenting. See docs/observable_ledger.md for the full equation ledger.
#
# STRUCTURE-DEPENDENT CONVENTION (TableIVPolynomial vs LeadingS0)
#   Table IV writes the structure-dependent classes as S0 - k A qbar^2 with
#   k = 1/2, 3/10, 3/4 for S, D, P. That is `TableIVPolynomial()`. But the paper's NUMERIC
#   column was computed with only the LEADING constant S0 = 3 h beta (the
#   -k A qbar^2 polynomial dropped); that is `LeadingS0()`, and it is what
#   reproduces both the reported S0 ~ 3.27 fit and the tabulated D/P numbers to
#   ~1%. `TableIVPolynomial()` is the default (faithful to the printed formula);
#   the reproduction harness uses `LeadingS0()`.
#
# Qualified reference API (public in QuarkModelTransitions.jl):
#   StrongDecayModel, decay_momentum, reduced_decay_amplitude, spatial_overlap,
#   strong_decay_amplitude, calibrate_strong_decay_model,
#   DecayChannel, StrongDecayAmplitude, decay_amplitude, reduced_matrix_element,
#   decay_width, MesonMasses, mass

"""
    StrongDecayModel(A, S0, beta_GeV)

Two-parameter Table IV/V decay model: `A` is the structure-independent reduced
amplitude (fit to `rho -> pi pi`), `S0 = 3 h beta` the structure-dependent
strength (fit to `B -> [omega pi]_S`), with oscillator scale `beta` (0.40 GeV
in the paper).

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
@assert model.beta_GeV == 0.4
```

## Related

- `TableVReference` — evaluate a frozen paper row.
- `calibrate_strong_decay_model` — fit the two reference strengths.
- `decay_amplitude` — evaluate a reference row.
"""
struct StrongDecayModel
    A::Float64
    S0::Float64
    beta_GeV::Float64
end

"""Convention for evaluating the structure-dependent Table IV amplitudes."""
abstract type ReducedAmplitudeConvention end

"""
Use the full polynomial printed in Table IV.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
@assert QMT.reduced_decay_amplitude(model, :S, 1.0; convention=QMT.TableIVPolynomial()) ≈ 2.77
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
struct TableIVPolynomial <: ReducedAmplitudeConvention end

"""
Use only the leading `S0` term employed for the paper's numeric column.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
@assert QMT.reduced_decay_amplitude(model, :S, 1.0; convention=QMT.LeadingS0()) == model.S0
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
struct LeadingS0 <: ReducedAmplitudeConvention end

"""
    decay_momentum(M, m1, m2)

[DERIVED] Two-body breakup momentum of `M -> m1 + m2` (GeV); zero below
threshold. Pure Kallen kinematics, no paper input.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
@assert QMT.decay_momentum(0.77, 0.14, 0.14) > 0
@assert QMT.decay_momentum(0.2, 0.14, 0.14) == 0
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
function decay_momentum(M::Real, m1::Real, m2::Real)
    M <= m1 + m2 && return 0.0
    return sqrt((M^2 - (m1 + m2)^2) * (M^2 - (m1 - m2)^2)) / (2M)
end

"""
    reduced_decay_amplitude(model, class, qbar;
                            convention=TableIVPolynomial())

[PAPER form, DERIVED strengths] Table IV reduced partial-wave amplitude for
`class` at `qbar = q/beta`. Structure-independent classes (`:A`, `:Aprime`,
`:Adoubleprime`, `:A0`, `:A_c`) return the fitted `A` (the paper's
`A' ~ A'' ~ A0 ~ A` convention).

Structure-dependent classes (`:S`, `:S_c`, `:D`, `:P`) depend on `convention`:
- `TableIVPolynomial()` (default) returns `S0 - k A qbar^2` — the formula printed in
  Table IV, with `k = 3/10` for `:D` and `3/4` for `:P`.
- `LeadingS0()` returns the constant `S0` — the convention the paper actually used
  for the numeric column (drops the `-k A qbar^2` polynomial), reproducing the
  tabulated D/P amplitudes to ~1%.

For the S family, Table IV (page-image 13) prints

    S   = [3h - (1/2) (g + h/4) q^2/beta^2] beta
    S_c = [3h - (m_c A_c/((m_d+m_c) beta)) q^2/beta_c^2] beta_c

so `k` is `1/2` for `:S` and `r = m_c/(m_d+m_c)` for `:S_c` — and at equal
constituent masses `r = 1/2`, i.e. the light `S` *is* the equal-mass case of
`S_c`, exactly as the light Gaussian is the `r = 1/2` case of the charmed form
factor. Both therefore take `k = heavy_fraction`, which defaults to `0.5`.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
@assert QMT.reduced_decay_amplitude(model, :A, 0.9) == model.A
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
function reduced_decay_amplitude(
    model::StrongDecayModel, class::Symbol, qbar::Real;
    convention::ReducedAmplitudeConvention = TableIVPolynomial(),
    heavy_fraction::Real = 0.5,
)
    return reduced_decay_amplitude(convention, model, class, qbar; heavy_fraction)
end

function reduced_decay_amplitude(
    ::LeadingS0, model::StrongDecayModel, class::Symbol, qbar::Real;
    heavy_fraction::Real = 0.5,
)
    class in (:A, :Aprime, :Adoubleprime, :A0, :A_c) && return model.A
    class in (:S, :S_c, :D, :P) && return model.S0
    throw(ArgumentError("unknown reduced-amplitude class `$class`"))
end

function reduced_decay_amplitude(
    ::TableIVPolynomial, model::StrongDecayModel, class::Symbol, qbar::Real;
    heavy_fraction::Real = 0.5,
)
    class in (:A, :Aprime, :Adoubleprime, :A0, :A_c) && return model.A
    class in (:S, :S_c) && return model.S0 - heavy_fraction * model.A * qbar^2
    class === :D && return model.S0 - 0.3 * model.A * qbar^2
    class === :P && return model.S0 - 0.75 * model.A * qbar^2
    throw(ArgumentError("unknown reduced-amplitude class `$class`"))
end

"""
    spatial_overlap(q_GeV, L, beta_GeV; heavy_fraction=0.5, recoil=false,
                    beta_c_GeV=beta_GeV)

[DERIVED] The harmonic-oscillator momentum-space overlap factor

    qbar^L * sqrt(q/2pi) * exp(-q^2 / 16 beta^2),   qbar = q/beta,

in `MeV^(1/2)` (q enters the square root in MeV). This is the single-beta SHO
integral the paper suppresses from the Table V formula column; it carries the
whole `q`-dependence of an amplitude. Returns 0 at or below threshold.

`heavy_fraction` is `r = m_Q/(m_Q + m_q)`, the heavy quark's share of the
constituent mass, and it is the *only* thing that distinguishes an unequal-mass
row from a light one:

    form factor = exp[-(1/4) r^2 q^2 / beta_c^2]

At `r = 1/2` (equal constituent masses) this is exactly the light Gaussian
`exp(-q^2/16 beta^2)` — there is no separate light branch, and the two agree
bitwise. Charmed rows (Table V footnote d) pass `r = m_c/(m_c+m_d) ≈ 0.881`;
`b`-flavored rows would simply pass a larger `r`. The A_c P-wave rows
additionally carry the unequal-mass recoil multiplier
`m_c beta / ((m_c+m_d) beta_c)` when `recoil=true`.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
@assert QMT.spatial_overlap(0.36, 1, 0.4) > 0
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
function spatial_overlap(
    q_GeV::Real, L::Integer, beta_GeV::Real;
    heavy_fraction::Real = 0.5, recoil::Bool = false,
    beta_c_GeV::Real = beta_GeV,
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / beta_c_GeV
    form_factor = exp(-0.25 * heavy_fraction^2 * q_GeV^2 / beta_c_GeV^2)
    recoil_mult = recoil ? heavy_fraction * beta_GeV / beta_c_GeV : 1.0
    return qbar^L * sqrt(1000q_GeV / (2π)) * form_factor * recoil_mult
end

# Back-compat private aliases (used by the scalar amplitude functions below).
_suppressed_factor(q_GeV::Real, beta_GeV::Real) = spatial_overlap(q_GeV, 0, beta_GeV)

"""
    strong_decay_amplitude(model, coefficient, class, qbar_power, q_GeV;
                           heavy_fraction=0.5, recoil=false,
                           beta_c_GeV=model.beta_GeV,
                           convention=TableIVPolynomial())

Scalar (compat) Table V amplitude in `MeV^(1/2)`: `coefficient` is the signed
flavor/spin factor, `class` the reduced-amplitude class, `qbar_power` the
explicit `qbar^L` power, `q_GeV` the breakup momentum. Prefer the row-oriented
`decay_amplitude` for new code; this returns only the product.

`heavy_fraction` is `r = m_Q/(m_Q + m_q)`; the default `0.5` is the equal-mass
(light) case. Charmed rows pass `r = m_c/(m_c+m_d)` and, on the A_c P-waves,
`recoil=true` for the footnote-d multiplier `r * beta/beta_c`. There is no
separate charm entry point — the mass ratio is the only difference.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
@assert QMT.strong_decay_amplitude(model, sqrt(4/3), :A, 1, 0.36) > 0
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
function strong_decay_amplitude(
    model::StrongDecayModel,
    coefficient::Real,
    class::Symbol,
    qbar_power::Integer,
    q_GeV::Real;
    heavy_fraction::Real = 0.5,
    recoil::Bool = false,
    beta_c_GeV::Real = model.beta_GeV,
    convention::ReducedAmplitudeConvention = TableIVPolynomial(),
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / beta_c_GeV
    return coefficient *
           reduced_decay_amplitude(model, class, qbar; convention, heavy_fraction) *
           spatial_overlap(q_GeV, qbar_power, model.beta_GeV;
               heavy_fraction = heavy_fraction, recoil = recoil,
               beta_c_GeV = beta_c_GeV)
end

"""
    calibrate_strong_decay_model(rho_q_GeV, B_q_GeV; rho_amplitude, B_amplitude,
                                 beta_GeV,
                                 convention=TableIVPolynomial())

Fix `A` from `rho -> pi pi` (`+(4/3)^(1/2) A qbar`, paper `+12.4 MeV^(1/2)`)
and then `S0` from `B -> [omega pi]_S` (`-(2/9)^(1/2) S(qbar)`, paper `-11`),
given the breakup momenta of the two fit decays. `A` is convention-independent.
`S0` differs: `TableIVPolynomial()` back-solves the full
`S0 - (1/2) A qbar_B^2`, while `LeadingS0()` sets `S0` directly to the value at
`q_B` (the paper's numeric convention, giving `S0 ~ 3.27`).

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.calibrate_strong_decay_model(0.36, 0.35; convention=QMT.LeadingS0())
@assert model.A > 0
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
function calibrate_strong_decay_model(
    rho_q_GeV::Real,
    B_q_GeV::Real;
    rho_amplitude::Real = STRONG_DECAY_DEFAULTS.rho_amplitude,
    B_amplitude::Real = STRONG_DECAY_DEFAULTS.B_amplitude,
    beta_GeV::Real = STRONG_DECAY_DEFAULTS.beta_GeV,
    convention::ReducedAmplitudeConvention = TableIVPolynomial(),
)
    qbar_rho = rho_q_GeV / beta_GeV
    A = rho_amplitude /
        (sqrt(4 / 3) * qbar_rho * _suppressed_factor(rho_q_GeV, beta_GeV))
    qbar_B = B_q_GeV / beta_GeV
    S_at_B = B_amplitude / (-sqrt(2 / 9) * _suppressed_factor(B_q_GeV, beta_GeV))
    S0 = _calibrated_S0(convention, S_at_B, A, qbar_B)
    return StrongDecayModel(A, S0, beta_GeV)
end

_calibrated_S0(::LeadingS0, S_at_B, A, qbar_B) = S_at_B
_calibrated_S0(::TableIVPolynomial, S_at_B, A, qbar_B) = S_at_B + 0.5 * A * qbar_B^2

# =============================================================================
# Row-oriented API: DecayChannel -> StrongDecayAmplitude
# =============================================================================

"""
    DecayChannel(parent, daughter1, daughter2, coefficient, class, qbar_power;
                 label="", section="")

A named strong-decay channel. `parent`/`daughter1`/`daughter2` are meson names
resolved to masses by a `MesonMasses`; `coefficient` is the signed
flavor-spin factor `c`, `class` the reduced-amplitude class, and `qbar_power`
the orbital power `L`. Paper-specific row loaders live in the comparison layer;
channels can always be constructed directly.

`heavy_fraction` is `r = m_Q/(m_Q + m_q)`, the constituent-mass ratio driving the
form factor (see `spatial_overlap`). It defaults to `0.5` (equal masses,
the light rows) but is **required** for the unequal-mass `:A_c`/`:S_c` classes —
constructing one without it throws rather than silently applying the light
Gaussian. Unlike the other fields it is *derived*, not digitized: the canonical
CSV has no quark-content column, so the caller resolves it from the parent's
flavor content. The A_c P-wave recoil multiplier is still keyed on `class`.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
@assert channel.qbar_power == 1
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
struct DecayChannel
    parent::String
    daughter1::String
    daughter2::String
    coefficient::Float64
    class::Symbol
    qbar_power::Int
    label::String
    section::String
    heavy_fraction::Float64
end

function DecayChannel(
    parent, daughter1, daughter2, coefficient, class, qbar_power;
    label::AbstractString = "", section::AbstractString = "",
    heavy_fraction::Union{Nothing,Real} = nothing,
)
    cls = Symbol(class)
    hf = if heavy_fraction === nothing
        # Loud, not silent: an unequal-mass row that falls back to the equal-mass
        # default would quietly get the light form factor (0.5 vs ~0.881).
        is_charm_class(cls) && throw(ArgumentError(
            "channel `$parent -> $daughter1 $daughter2` has unequal-mass class " *
            "`$cls` but no `heavy_fraction`; pass r = m_Q/(m_Q+m_q) explicitly",
        ))
        0.5
    else
        Float64(heavy_fraction)
    end
    0 < hf < 1 || throw(ArgumentError(
        "heavy_fraction is m_Q/(m_Q+m_q) and must lie in (0,1), got $hf",
    ))
    return DecayChannel(String(parent), String(daughter1), String(daughter2),
        Float64(coefficient), cls, Int(qbar_power),
        String(label), String(section), hf)
end

is_charm_class(class::Symbol) = class in (:A_c, :S_c)
# Footnote d: only the A_c P-wave rows (qbar^L, L >= 2) carry the recoil factor.
channel_has_recoil(ch::DecayChannel) = ch.class === :A_c && ch.qbar_power >= 2

"""
    StrongDecayAmplitude

The factorized decomposition of a Table V amplitude:
- `coefficient` — dimensionless flavor-spin factor `c` [PAPER];
- `reduced` — dimensionless Table IV reduced amplitude `X(qbar)`;
- `spatial_overlap` — SHO momentum overlap in `MeV^(1/2)` [DERIVED];
- `total` — their product in `MeV^(1/2)` (the tabulated amplitude);
- `q_GeV` — the breakup momentum used.

See `reduced_matrix_element` (`= coefficient*reduced`) and
[`decay_width`](@ref) (`= total^2`).

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
amplitude = QMT.decay_amplitude(model, channel, 0.36)
@assert amplitude isa QMT.StrongDecayAmplitude
```

## Related

- `TableVReference` — evaluate a frozen paper row.
- `decay_amplitude` — evaluate a reference row.
"""
struct StrongDecayAmplitude
    coefficient::Float64
    reduced::Float64
    spatial_overlap::Float64
    total::Float64
    q_GeV::Float64
end

"""
    reduced_matrix_element(a::StrongDecayAmplitude)

The dimensionless flavor-spin/reduced interaction factor
`coefficient * reduced`, before the spatial overlap.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
amplitude = QMT.decay_amplitude(model, channel, 0.36)
@assert QMT.reduced_matrix_element(amplitude) ≈ sqrt(4/3)
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
reduced_matrix_element(a::StrongDecayAmplitude) = a.coefficient * a.reduced


decay_width(a::StrongDecayAmplitude) = a.total^2

"""
    decay_amplitude(model, ch::DecayChannel, q_GeV::Real;
                    convention=LeadingS0())
    decay_amplitude(model, ch::DecayChannel, masses::MesonMasses;
                    convention=LeadingS0())

Row-oriented Table V amplitude, returning the full `StrongDecayAmplitude`
decomposition. With a `MesonMasses` the breakup momentum is resolved internally
from the channel's parent/daughter names, so nothing is passed positionally.
Charmed rows (`:A_c`/`:S_c`) automatically use the footnote-d form factor and
(for A_c P-waves) the recoil multiplier. Defaults to `LeadingS0()` (the
paper's numeric column); pass `convention=TableIVPolynomial()` for the printed
Table IV formula.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
amplitude = QMT.decay_amplitude(model, channel, 0.36)
@assert decay_width(amplitude) ≈ abs2(amplitude.total)
```

## Related

- `TableVReference` — evaluate a frozen paper row.
- [`decay_width`](@ref) — convert a transition result to MeV.
- `DecayChannel` — one explicitly supplied reference row.
"""
function decay_amplitude(
    model::StrongDecayModel, ch::DecayChannel, q_GeV::Real;
    convention::ReducedAmplitudeConvention = LeadingS0(),
)
    hf = ch.heavy_fraction
    qbar = q_GeV / model.beta_GeV
    reduced = reduced_decay_amplitude(model, ch.class, qbar; convention, heavy_fraction = hf)
    overlap = spatial_overlap(q_GeV, ch.qbar_power, model.beta_GeV;
        heavy_fraction = hf, recoil = channel_has_recoil(ch))
    total = ch.coefficient * reduced * overlap
    return StrongDecayAmplitude(ch.coefficient, reduced, overlap, total, Float64(q_GeV))
end

"""
    TableVReference(model, channel, partial_wave; convention=LeadingS0())

Frozen Table IV/V reference backend for the typed transition API. The explicit
`partial_wave` records information already selected by the paper row; it is not
used by solver-native operators, which derive all allowed waves together.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
model = QMT.StrongDecayModel(1.0, 3.27, 0.4)
channel = QMT.DecayChannel("rho", "pi+", "pi-", sqrt(4/3), :A, 1)
initial = QMT.ReferenceState("rho", 0.77; J=1, parity=-1)
final = TwoMesonChannel(QMT.ReferenceState("pi+", 0.14; J=0, parity=-1), QMT.ReferenceState("pi-", 0.14; J=0, parity=-1))
operator = QMT.TableVReference(model, channel, PartialWave(1, 0))
amplitude = matrix_element(final, operator, initial)
@assert decay_width(amplitude) > 0
```

## Related

- `DecayChannel` — one explicitly supplied reference row.
- `GITableVNormalization` — squared reference amplitude convention.
- `LeadingS0` — leading constant convention.
- `MesonMasses` — reference meson mass lookup.
- `ReferenceState` — wave-free paper-reference state.
- `STRONG_DECAY_DEFAULTS` — reference calibration inputs.
- `StrongDecayAmplitude` — reference-row factorization.
- `StrongDecayModel` — frozen Table IV/V parameters.
- `TableIVPolynomial` — full printed polynomial convention.
- `calibrate_strong_decay_model` — fit the two reference strengths.
- `decay_amplitude` — evaluate a reference row.
- `decay_momentum` — two-body breakup momentum in GeV.
- [`matrix_element`](@ref) — evaluate an operator between states.
- `meson_mass` — retrieve a reference mass in GeV.
- `reduced_decay_amplitude` — Table IV reduced factor.
- `reduced_matrix_element` — inspect the nonspatial reference factor.
- `spatial_overlap` — reference SHO overlap factor.
- `strong_decay_amplitude` — scalar reference-row compatibility helper.
"""
struct TableVReference{C<:ReducedAmplitudeConvention} <: StrongDecayOperator
    model::StrongDecayModel
    channel::DecayChannel
    partial_wave::PartialWave
    convention::C
end

function TableVReference(
    model::StrongDecayModel,
    channel::DecayChannel,
    partial_wave::PartialWave;
    convention::ReducedAmplitudeConvention = LeadingS0(),
)
    partial_wave.relative_L == channel.qbar_power || throw(ArgumentError(
        "legacy partial-wave L=$(partial_wave.relative_L) disagrees with " *
        "DecayChannel qbar_power=$(channel.qbar_power)",
    ))
    return TableVReference{typeof(convention)}(model, channel, partial_wave, convention)
end

function _validate_transition(
    final::TwoMesonChannel,
    operator::TableVReference,
    initial::TransitionState,
)
    initial.label == operator.channel.parent || throw(ArgumentError(
        "reference parent `$(initial.label)` does not match row parent " *
        "`$(operator.channel.parent)`",
    ))
    actual = sort([final.first.label, final.second.label])
    expected = sort([operator.channel.daughter1, operator.channel.daughter2])
    actual == expected || throw(ArgumentError(
        "reference daughters $(join(actual, ", ")) do not match row daughters " *
        "$(join(expected, ", "))",
    ))
    return nothing
end

function matrix_element(
    final::TwoMesonChannel,
    operator::TableVReference,
    initial::TransitionState,
)
    _validate_transition(final, operator, initial)
    resolved = _resolve_kinematics(final, initial, OnShell())
    q = resolved.momentum_GeV
    q isa Real || throw(ArgumentError(
        "TableVReference has no complex-momentum continuation; use a native operator",
    ))
    q > 0 || throw(ArgumentError(
        "TableVReference requires positive real momentum; zero is reserved for closed widths",
    ))
    legacy = decay_amplitude(operator.model, operator.channel, q;
        convention = operator.convention)
    term = TransitionTerm(
        "Table V row",
        legacy.coefficient,
        legacy.reduced,
        legacy.spatial_overlap,
        legacy.total,
        (source = :paper, section = operator.channel.section),
    )
    partial_wave_amplitudes = (operator.partial_wave => legacy.total,)
    provenance = (
        backend = :table_v_reference,
        convention = nameof(typeof(operator.convention)),
        helicity_available = false,
        row_label = operator.channel.label,
    )
    return TransitionAmplitude(
        operator,
        initial,
        final,
        resolved,
        GITableVNormalization(),
        (),
        partial_wave_amplitudes,
        (term,),
        provenance,
    )
end

partial_width(::GITableVNormalization, amplitude::TransitionAmplitude) =
    sum(abs2(last(item)) for item in amplitude.partial_wave_amplitudes)

# =============================================================================
# MesonMasses: name -> mass resolver
# =============================================================================

"""
    MesonMasses(lookup::Dict{String,Float64})

A simple meson-name -> mass (GeV) resolver shared by every `DecayChannel`. Built
by the reproduction harness from experimental values and/or model-predicted
masses (`compute_spectrum`). Access with `meson_mass` or indexing.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
masses = QMT.MesonMasses(Dict("rho" => 0.77, "pi" => 0.14))
@assert masses["rho"] == 0.77
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
struct MesonMasses
    lookup::Dict{String,Float64}
end
MesonMasses() = MesonMasses(Dict{String,Float64}())

"""
    meson_mass(m::MesonMasses, name) -> Float64

Mass (GeV) of a meson by name; throws with a clear message if the name is not
registered.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
masses = QMT.MesonMasses(Dict("rho" => 0.77))
@assert QMT.meson_mass(masses, "rho") == 0.77
```

## Related

- `TableVReference` — evaluate a frozen paper row.
"""
function meson_mass(m::MesonMasses, name::AbstractString)
    haskey(m.lookup, name) || throw(KeyError("no mass registered for meson `$name`"))
    return m.lookup[name]
end
Base.getindex(m::MesonMasses, name::AbstractString) = meson_mass(m, name)
Base.setindex!(m::MesonMasses, v::Real, name::AbstractString) = (m.lookup[name] = Float64(v))
Base.haskey(m::MesonMasses, name::AbstractString) = haskey(m.lookup, name)

function decay_amplitude(
    model::StrongDecayModel, ch::DecayChannel, masses::MesonMasses;
    convention::ReducedAmplitudeConvention = LeadingS0(),
)
    q = decay_momentum(meson_mass(masses, ch.parent),
        meson_mass(masses, ch.daughter1), meson_mass(masses, ch.daughter2))
    return decay_amplitude(model, ch, q; convention)
end

# --- Display -----------------------------------------------------------------
# The point of `StrongDecayAmplitude` is that the three factors stay separate;
# print them that way, with the product shown as the product it is.

function Base.show(io::IO, ::MIME"text/plain", a::StrongDecayAmplitude)
    println(io, "StrongDecayAmplitude  (q = ", @sprintf("%.4f", a.q_GeV), " GeV)")
    println(io, "  coefficient      ", @sprintf("%10.4f", a.coefficient), "   [PAPER]  App. B flavor-spin")
    println(io, "  reduced          ", @sprintf("%10.4f", a.reduced), "   [PAPER]  Table IV class")
    println(io, "  spatial_overlap  ", @sprintf("%10.4f", a.spatial_overlap), "   [DERIVED] SHO integral")
    println(io, "  " * "-"^58)
    println(io, "  total            ", @sprintf("%10.4f", a.total), "   MeV^(1/2)  = the product")
    print(io, "  width            ", @sprintf("%10.4f", decay_width(a)), "   MeV       = total^2")
    return nothing
end

Base.show(io::IO, a::StrongDecayAmplitude) =
    print(io, "StrongDecayAmplitude(total = ", @sprintf("%.4f", a.total), " MeV^(1/2))")

function Base.show(io::IO, ::MIME"text/plain", ch::DecayChannel)
    reaction = string(ch.parent, " -> ", ch.daughter1, " + ", ch.daughter2)
    println(io, "DecayChannel: ", reaction)
    isempty(ch.label) || ch.label == reaction || println(io, "  label          ", ch.label)
    println(io, "  coefficient    ", ch.coefficient)
    println(io, "  class          :", ch.class, "   (Table IV)")
    println(io, "  qbar_power L   ", ch.qbar_power)
    print(io, "  heavy_fraction ", ch.heavy_fraction,
        ch.heavy_fraction == 0.5 ? "   (equal-mass: the light Gaussian)" : "   r = m_Q/(m_Q+m_q)")
    isempty(ch.section) || print(io, "\n  section        ", ch.section)
    return nothing
end

Base.show(io::IO, ch::DecayChannel) = print(
    io, "DecayChannel(", ch.parent, " -> ", ch.daughter1, " + ", ch.daughter2, ", :", ch.class, ")",
)
