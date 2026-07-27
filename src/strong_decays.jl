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
# STRUCTURE-DEPENDENT CONVENTION (:table_iv vs :leading)
#   Table IV writes the structure-dependent classes as S0 - k A qbar^2 with
#   k = 1/2, 3/10, 3/4 for S, D, P. That is `:table_iv`. But the paper's NUMERIC
#   column was computed with only the LEADING constant S0 = 3 h beta (the
#   -k A qbar^2 polynomial dropped); that is `:leading`, and it is what
#   reproduces both the reported S0 ~ 3.27 fit and the tabulated D/P numbers to
#   ~1%. `:table_iv` is the default (faithful to the table's printed formula);
#   the reproduction harness uses `:leading`.
#
# Public API (exported from GIModel.jl):
#   StrongDecayModel, decay_momentum, reduced_decay_amplitude, spatial_overlap,
#   strong_decay_amplitude, charm_decay_amplitude, calibrate_strong_decay_model,
#   DecayChannel, StrongDecayAmplitude, decay_amplitude, matrix_element,
#   decay_width, MesonMasses, mass, load_table_v

"""
    StrongDecayModel(A, S0, beta_GeV)

Two-parameter Table IV/V decay model: `A` is the structure-independent reduced
amplitude (fit to `rho -> pi pi`), `S0 = 3 h beta` the structure-dependent
strength (fit to `B -> [omega pi]_S`), with oscillator scale `beta` (0.40 GeV
in the paper).
"""
struct StrongDecayModel
    A::Float64
    S0::Float64
    beta_GeV::Float64
end

"""
    decay_momentum(M, m1, m2)

[DERIVED] Two-body breakup momentum of `M -> m1 + m2` (GeV); zero below
threshold. Pure Kallen kinematics, no paper input.
"""
function decay_momentum(M::Real, m1::Real, m2::Real)
    M <= m1 + m2 && return 0.0
    return sqrt((M^2 - (m1 + m2)^2) * (M^2 - (m1 - m2)^2)) / (2M)
end

"""
    reduced_decay_amplitude(model, class, qbar; convention=:table_iv)

[PAPER form, DERIVED strengths] Table IV reduced partial-wave amplitude for
`class` at `qbar = q/beta`. Structure-independent classes (`:A`, `:Aprime`,
`:Adoubleprime`, `:A0`, `:A_c`) return the fitted `A` (the paper's
`A' ~ A'' ~ A0 ~ A` convention).

Structure-dependent classes (`:S`, `:S_c`, `:D`, `:P`) depend on `convention`:
- `:table_iv` (default) returns `S0 - k A qbar^2` with `k = 1/2, 3/10, 3/4` — the
  formula printed in Table IV.
- `:leading` returns the constant `S0` — the convention the paper actually used
  for the numeric column (drops the `-k A qbar^2` polynomial), reproducing the
  tabulated D/P amplitudes to ~1%.
"""
function reduced_decay_amplitude(
    model::StrongDecayModel, class::Symbol, qbar::Real; convention::Symbol = :table_iv,
)
    class in (:A, :Aprime, :Adoubleprime, :A0, :A_c) && return model.A
    if convention === :leading
        class in (:S, :S_c, :D, :P) && return model.S0
    elseif convention === :table_iv
        class in (:S, :S_c) && return model.S0 - 0.5 * model.A * qbar^2
        class === :D && return model.S0 - 0.3 * model.A * qbar^2
        class === :P && return model.S0 - 0.75 * model.A * qbar^2
    else
        throw(ArgumentError("unknown convention `$convention` (:table_iv | :leading)"))
    end
    throw(ArgumentError("unknown reduced-amplitude class `$class`"))
end

# Default charm-meson quark masses (GeV): m_c and the light spectator m_d, from
# data/parameters.provisional.toml (Table II). Used by the charm form factor.
const CHARM_M_C_GEV = 1.628
const CHARM_M_D_GEV = 0.220

"""
    spatial_overlap(q_GeV, L, beta_GeV; heavy_fraction=0.5, recoil=false,
                    m_c_GeV=CHARM_M_C_GEV, m_d_GeV=CHARM_M_D_GEV, beta_c_GeV=beta_GeV)

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
"""
function spatial_overlap(
    q_GeV::Real, L::Integer, beta_GeV::Real;
    heavy_fraction::Real = 0.5, recoil::Bool = false,
    m_c_GeV::Real = CHARM_M_C_GEV, m_d_GeV::Real = CHARM_M_D_GEV,
    beta_c_GeV::Real = beta_GeV,
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / beta_c_GeV
    form_factor = exp(-0.25 * heavy_fraction^2 * q_GeV^2 / beta_c_GeV^2)
    recoil_mult = recoil ? m_c_GeV * beta_GeV / ((m_c_GeV + m_d_GeV) * beta_c_GeV) : 1.0
    return qbar^L * sqrt(1000q_GeV / (2π)) * form_factor * recoil_mult
end

# The charm heavy fraction while the class still selects it (removed in S4,
# when DecayChannel carries the fraction directly).
_charm_heavy_fraction() = CHARM_M_C_GEV / (CHARM_M_C_GEV + CHARM_M_D_GEV)

# Back-compat private aliases (used by the scalar amplitude functions below).
_suppressed_factor(q_GeV::Real, beta_GeV::Real) = spatial_overlap(q_GeV, 0, beta_GeV)

"""
    strong_decay_amplitude(model, coefficient, class, qbar_power, q_GeV; convention=:table_iv)

Scalar (compat) Table V amplitude in `MeV^(1/2)`: `coefficient` is the signed
flavor/spin factor, `class` the reduced-amplitude class, `qbar_power` the
explicit `qbar^L` power, `q_GeV` the breakup momentum. Prefer the row-oriented
[`decay_amplitude`](@ref) for new code; this returns only the product.
"""
function strong_decay_amplitude(
    model::StrongDecayModel,
    coefficient::Real,
    class::Symbol,
    qbar_power::Integer,
    q_GeV::Real;
    convention::Symbol = :table_iv,
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / model.beta_GeV
    return coefficient *
           reduced_decay_amplitude(model, class, qbar; convention) *
           spatial_overlap(q_GeV, qbar_power, model.beta_GeV)
end

"""
    charm_decay_amplitude(model, coefficient, class, qbar_power, q_GeV;
                          recoil=false, m_c_GeV=CHARM_M_C_GEV, m_d_GeV=CHARM_M_D_GEV,
                          beta_c_GeV=model.beta_GeV, convention=:table_iv)

Scalar (compat) Table V charmed-meson amplitude in `MeV^(1/2)` (footnote d):
the charmed Gaussian form factor with `beta_c = beta` numerically, the
`:A_c` / `:S_c` classes, and (with `recoil=true`) the unequal-mass multiplier
`m_c beta / ((m_c+m_d) beta_c)` printed on the A_c P-wave rows.
"""
function charm_decay_amplitude(
    model::StrongDecayModel,
    coefficient::Real,
    class::Symbol,
    qbar_power::Integer,
    q_GeV::Real;
    recoil::Bool = false,
    m_c_GeV::Real = CHARM_M_C_GEV,
    m_d_GeV::Real = CHARM_M_D_GEV,
    beta_c_GeV::Real = model.beta_GeV,
    convention::Symbol = :table_iv,
)
    q_GeV <= 0 && return 0.0
    qbar = q_GeV / beta_c_GeV
    return coefficient *
           reduced_decay_amplitude(model, class, qbar; convention) *
           spatial_overlap(q_GeV, qbar_power, model.beta_GeV;
               heavy_fraction = m_c_GeV / (m_c_GeV + m_d_GeV), recoil = recoil,
               m_c_GeV = m_c_GeV, m_d_GeV = m_d_GeV, beta_c_GeV = beta_c_GeV)
end

"""
    calibrate_strong_decay_model(rho_q_GeV, B_q_GeV; rho_amplitude, B_amplitude,
                                 beta_GeV, convention=:table_iv)

Fix `A` from `rho -> pi pi` (`+(4/3)^(1/2) A qbar`, paper `+12.4 MeV^(1/2)`)
and then `S0` from `B -> [omega pi]_S` (`-(2/9)^(1/2) S(qbar)`, paper `-11`),
given the breakup momenta of the two fit decays. `A` is convention-independent.
`S0` differs: `:table_iv` back-solves the full `S0 - (1/2) A qbar_B^2`, while
`:leading` sets `S0` directly to the value at `q_B` (the paper's numeric
convention, giving `S0 ~ 3.27`).
"""
function calibrate_strong_decay_model(
    rho_q_GeV::Real,
    B_q_GeV::Real;
    rho_amplitude::Real = 12.4,
    B_amplitude::Real = -11.0,
    beta_GeV::Real = 0.40,
    convention::Symbol = :table_iv,
)
    qbar_rho = rho_q_GeV / beta_GeV
    A = rho_amplitude /
        (sqrt(4 / 3) * qbar_rho * _suppressed_factor(rho_q_GeV, beta_GeV))
    qbar_B = B_q_GeV / beta_GeV
    S_at_B = B_amplitude / (-sqrt(2 / 9) * _suppressed_factor(B_q_GeV, beta_GeV))
    S0 = convention === :leading ? S_at_B : S_at_B + 0.5 * A * qbar_B^2
    return StrongDecayModel(A, S0, beta_GeV)
end

# =============================================================================
# Row-oriented API: DecayChannel -> StrongDecayAmplitude
# =============================================================================

"""
    DecayChannel(parent, daughter1, daughter2, coefficient, class, qbar_power;
                 label="", section="")

A named Table V decay row. `parent`/`daughter1`/`daughter2` are meson names
resolved to masses by a [`MesonMasses`](@ref); `coefficient` is the signed
flavor-spin factor `c` [PAPER], `class` the reduced-amplitude class, `qbar_power`
the orbital power `L`. Built by [`load_table_v`](@ref) from the canonical CSV or
by hand. Charmed-row behavior (footnote d form factor + A_c P-wave recoil) is
derived from `class`.
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
end

function DecayChannel(
    parent, daughter1, daughter2, coefficient, class, qbar_power;
    label::AbstractString = "", section::AbstractString = "",
)
    return DecayChannel(String(parent), String(daughter1), String(daughter2),
        Float64(coefficient), Symbol(class), Int(qbar_power),
        String(label), String(section))
end

is_charm_class(class::Symbol) = class in (:A_c, :S_c)
# Footnote d: only the A_c P-wave rows (qbar^L, L >= 2) carry the recoil factor.
channel_has_recoil(ch::DecayChannel) = ch.class === :A_c && ch.qbar_power >= 2

"""
    StrongDecayAmplitude

The 3-way decomposition of a Table V amplitude, all in `MeV^(1/2)`:
- `coefficient` — the flavor-spin factor `c` [PAPER];
- `reduced` — the Table IV reduced amplitude `X(qbar)`;
- `spatial_overlap` — the SHO momentum overlap [DERIVED];
- `total` — their product (the tabulated amplitude);
- `q_GeV` — the breakup momentum used.

See [`matrix_element`](@ref) (`= coefficient*reduced`) and [`decay_width`](@ref)
(`= total^2`).
"""
struct StrongDecayAmplitude
    coefficient::Float64
    reduced::Float64
    spatial_overlap::Float64
    total::Float64
    q_GeV::Float64
end

"""
    matrix_element(a::StrongDecayAmplitude)

The flavor-spin matrix element `coefficient * reduced` (MeV^(1/2)) — the
interaction part of the amplitude, before the spatial overlap.
"""
matrix_element(a::StrongDecayAmplitude) = a.coefficient * a.reduced

"""
    decay_width(a::StrongDecayAmplitude)
    decay_width(model, ch, q_or_masses; convention=:leading)

Partial width `|amplitude|^2` in MeV (GI normalization: the tabulated MeV^(1/2)
amplitude squared is the partial width).
"""
decay_width(a::StrongDecayAmplitude) = a.total^2

"""
    decay_amplitude(model, ch::DecayChannel, q_GeV::Real; convention=:leading)
    decay_amplitude(model, ch::DecayChannel, masses::MesonMasses; convention=:leading)

Row-oriented Table V amplitude, returning the full [`StrongDecayAmplitude`]
decomposition. With a `MesonMasses` the breakup momentum is resolved internally
from the channel's parent/daughter names, so nothing is passed positionally.
Charmed rows (`:A_c`/`:S_c`) automatically use the footnote-d form factor and
(for A_c P-waves) the recoil multiplier. Defaults to the `:leading` convention
(the paper's numeric column); pass `convention=:table_iv` for the printed
Table IV formula.
"""
function decay_amplitude(
    model::StrongDecayModel, ch::DecayChannel, q_GeV::Real;
    convention::Symbol = :leading,
)
    hf = is_charm_class(ch.class) ? _charm_heavy_fraction() : 0.5
    qbar = q_GeV / model.beta_GeV
    reduced = reduced_decay_amplitude(model, ch.class, qbar; convention)
    overlap = spatial_overlap(q_GeV, ch.qbar_power, model.beta_GeV;
        heavy_fraction = hf, recoil = channel_has_recoil(ch))
    total = ch.coefficient * reduced * overlap
    return StrongDecayAmplitude(ch.coefficient, reduced, overlap, total, Float64(q_GeV))
end

# =============================================================================
# MesonMasses: name -> mass resolver
# =============================================================================

"""
    MesonMasses(lookup::Dict{String,Float64})

A simple meson-name -> mass (GeV) resolver shared by every `DecayChannel`. Built
by the reproduction harness from experimental values and/or model-predicted
masses (`compute_spectrum`). Access with [`mass`](@ref) or indexing.
"""
struct MesonMasses
    lookup::Dict{String,Float64}
end
MesonMasses() = MesonMasses(Dict{String,Float64}())

"""
    mass(m::MesonMasses, name) -> Float64

Mass (GeV) of a meson by name; throws with a clear message if the name is not
registered.
"""
function mass(m::MesonMasses, name::AbstractString)
    haskey(m.lookup, name) || throw(KeyError("no mass registered for meson `$name`"))
    return m.lookup[name]
end
Base.getindex(m::MesonMasses, name::AbstractString) = mass(m, name)
Base.setindex!(m::MesonMasses, v::Real, name::AbstractString) = (m.lookup[name] = Float64(v))
Base.haskey(m::MesonMasses, name::AbstractString) = haskey(m.lookup, name)

function decay_amplitude(
    model::StrongDecayModel, ch::DecayChannel, masses::MesonMasses;
    convention::Symbol = :leading,
)
    q = decay_momentum(mass(masses, ch.parent),
        mass(masses, ch.daughter1), mass(masses, ch.daughter2))
    return decay_amplitude(model, ch, q; convention)
end

# =============================================================================
# Canonical Table V loader
# =============================================================================

# Map the canonical CSV `amp_class` strings to reduced-amplitude class symbols.
const _CLASS_SYMBOL = Dict(
    "A" => :A, "Aprime" => :Aprime, "Adoubleprime" => :Adoubleprime, "A0" => :A0,
    "S" => :S, "D" => :D, "P" => :P, "Ac" => :A_c, "Sc" => :S_c,
)

"""
    load_table_v(rows) -> Vector{DecayChannel}

Turn parsed canonical `table_v_strong_decays.csv` rows into `DecayChannel`s.
Accepts any row iterator (e.g. `CSV.File(path)`) so `src/` takes no CSV
dependency. Rows with a non-amplitude class (`mixing_only`, `unlisted`) are
skipped. Each row needs the columns `parent`, `daughter1`, `daughter2`,
`coefficient`, `amp_class`, `qbar_power`, `decay`, `section`.
"""
function load_table_v(rows)
    channels = DecayChannel[]
    for row in rows
        class_str = String(row.amp_class)
        haskey(_CLASS_SYMBOL, class_str) || continue  # skip mixing_only / unlisted
        push!(channels, DecayChannel(
            String(row.parent), String(row.daughter1), String(row.daughter2),
            Float64(row.coefficient), _CLASS_SYMBOL[class_str], Int(row.qbar_power);
            label = String(row.decay), section = String(row.section),
        ))
    end
    return channels
end
