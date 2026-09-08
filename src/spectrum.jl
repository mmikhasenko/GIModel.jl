# Pure model spectrum for one meson. The production path solves each complete
# fixed-(L,S,J) Hamiltonian once, then applies cross-sector spectroscopic
# mixing. A central-only spectrum remains an independent diagnostic and is not
# a prerequisite for the physical calculation.
#
# Public API (exported from GIModel.jl):
#   spectrum_levels, StateMixing, CentralState, CorrectedState, MixedState,
#   Spectrum, CentralSpectrum, CorrectedSpectrum, MixedSpectrum,
#   central_spectrum, fixed_spectrum, add_intra_meson_mixing,
#   compute_spectrum, spectrum_state, parameters

"""
    spectrum_levels(nmax; L_labels=("S", "P", "D")) -> Vector{BasisState}

Enumerate the ``n\\,^{2S+1}L_J`` multiplets for `n in 1:nmax` and each orbital
letter: the spin singlet (`J = L`) and the spin triplets (`J = 1` for `S` waves,
`J in L-1:L+1` otherwise). Use as the `levels` argument of [`compute_spectrum`](@ref);
hand-built [`BasisState`](@ref) vectors work the same way for custom selections.
Entries have `flavors = nothing`: the same selection can be used for any
[`Meson`](@ref). Solved states carry explicit channel flavors.

## Example

```julia
using GIModel
levels = spectrum_levels(2; L_labels=("S", "P"))
levels[1].multiplicity, levels[1].flavors
```
"""
function spectrum_levels(nmax::Integer; L_labels = ("S", "P", "D"))
    nmax >= 1 || throw(ArgumentError("spectrum_levels: nmax must be ≥ 1, got $nmax"))
    levels = BasisState[]
    for L_label in L_labels
        haskey(L_SYMBOLS, String(L_label)) ||
            throw(ArgumentError("spectrum_levels: unknown orbital label `$L_label`"))
        L = L_SYMBOLS[String(L_label)]
        for n in 1:nmax
            push!(levels, BasisState(n, L_label, 1, L))
            triplet_J = L == 0 ? (1:1) : ((L-1):(L+1))
            for J in triplet_J
                push!(levels, BasisState(n, L_label, 3, J))
            end
        end
    end
    return levels
end

"""
    StateMixing

Provenance of one mixing transformation applied to a [`MixedState`](@ref):
the shared [`MixingResult`](@ref), the eigenstate column selected for this
member, and its pre-mixing mass. Derived properties expose `mechanism`,
`block_label`, `partner_labels`, `components`, `partner_masses_GeV`, and the
largest `offdiag_GeV` in the shared block without storing copies of the block
eigensystem or undefined pairwise angles for a multi-level block.
"""
struct StateMixing
    result::MixingResult
    eigenstate::Int
    unmixed_GeV::Float64
    function StateMixing(
        result::MixingResult,
        eigenstate::Integer,
        unmixed_GeV::Real,
    )
        column = Int(eigenstate)
        1 <= column <= length(result.masses) || throw(ArgumentError(
            "StateMixing: eigenstate column $column is outside the mixing result",
        ))
        mass = float(unmixed_GeV)
        isfinite(mass) || throw(ArgumentError("StateMixing: unmixed mass must be finite"))
        return new(result, column, mass)
    end
end

function _maximum_offdiag(matrix::AbstractMatrix)
    n = size(matrix, 1)
    n <= 1 && return 0.0
    return maximum(abs(matrix[i, j]) for i in 1:n for j in (i + 1):n)
end

function Base.getproperty(m::StateMixing, name::Symbol)
    name === :mechanism && return getfield(m, :result).block.mechanism
    name === :block_label && return getfield(m, :result).block.name
    name === :partner_labels && return [s.label for s in getfield(m, :result).block.basis]
    name === :components && return @view getfield(m, :result).vectors[:, getfield(m, :eigenstate)]
    name === :partner_masses_GeV && return getfield(m, :result).masses
    name === :offdiag_GeV && return _maximum_offdiag(getfield(m, :result).block.matrix)
    return getfield(m, name)
end

Base.propertynames(::StateMixing) = (
    :result, :eigenstate, :unmixed_GeV, :offdiag_GeV,
    :mechanism, :block_label, :partner_labels, :components, :partner_masses_GeV,
)

"""
    CentralState

One ``n\\,^{2S+1}L_J`` level after the spin-independent central solve only:
the quantum labels and the central eigenvalue `central_GeV`. Element type of
[`CentralSpectrum`](@ref); no spin-dependent data exists at this stage.

## Next steps

Read `state.central_GeV` for the central eigenvalue and `state.basis` for its
quantum numbers. With its parent spectrum `spec`, use `radial_wave(spec, state)`
to obtain the wavefunction. For masses including spin interactions and mixing,
run [`compute_spectrum`](@ref) with the same parameters and meson.

## Related

- [`central_spectrum`](@ref): central levels.
- [`spectrum_state`](@ref): select a level.
- [`radial_wave`](@ref): retrieve its wave.
- [`wave_norm`](@ref): check normalization.
"""
struct CentralState
    basis::BasisState
    central_GeV::Float64
end

function Base.getproperty(s::CentralState, name::Symbol)
    name === :n && return getfield(s, :basis).n
    name === :L && return getfield(s, :basis).L_label
    name === :multiplicity && return getfield(s, :basis).multiplicity
    name === :J && return getfield(s, :basis).J
    name === :label && return getfield(s, :basis).label
    return getfield(s, name)
end
Base.propertynames(::CentralState) =
    (:basis, :n, :L, :multiplicity, :J, :label, :central_GeV)

"""
    CorrectedState

One eigenstate of the complete fixed-sector Hamiltonian. `central_GeV` is the
central-Hamiltonian contribution in this spin-distorted eigenstate, not the
eigenvalue of a separate central-only solve. Together with the contact,
vector/Thomas spin-orbit, and tensor contributions it sums to `mass_GeV`.
The spectroscopic identity reuses [`BasisState`](@ref); no central-state holder
is embedded in the production state.

## Next steps

Read `state.mass_GeV` for the fixed-sector mass and `state.basis` for the quantum
numbers. `propertynames(state)` lists the energy contributions, including
`state.contact_shift_GeV` and `state.tensor_shift_GeV`.
Use `radial_wave(spec, state)` with the parent spectrum to obtain its radial wave.
Apply [`add_intra_meson_mixing`](@ref) to the whole [`CorrectedSpectrum`](@ref)
to obtain a [`MixedSpectrum`](@ref).

## Related

- [`fixed_spectrum`](@ref): fixed-sector levels.
- [`spectrum_state`](@ref): select a level.
- [`radial_wave`](@ref): retrieve its wave.
- [`compute_spectrum`](@ref): solve and mix.
"""
struct CorrectedState
    basis::BasisState
    central_GeV::Float64
    contact_shift_GeV::Float64
    spin_orbit_vector_shift_GeV::Float64
    spin_orbit_thomas_shift_GeV::Float64
    spin_orbit_shift_GeV::Float64
    tensor_shift_GeV::Float64
    fine_structure_shift_GeV::Float64
    fine_structure_mass_convention::String
    mass_GeV::Float64
end

function Base.getproperty(s::CorrectedState, name::Symbol)
    name === :n && return getfield(s, :basis).n
    name === :L && return getfield(s, :basis).L_label
    name === :multiplicity && return getfield(s, :basis).multiplicity
    name === :J && return getfield(s, :basis).J
    name === :label && return getfield(s, :basis).label
    return getfield(s, name)
end
Base.propertynames(::CorrectedState) = (
    fieldnames(CorrectedState)..., :n, :L, :multiplicity, :J, :label,
)

"""
    MixedState

One level of a [`MixedSpectrum`](@ref), obtained with [`spectrum_state`](@ref)
or `spec.states[i]`.

- `state.mass_GeV`: final mass.
- `state.corrected`: pre-mixing [`CorrectedState`](@ref).
- `state.mixings`: [`StateMixing`](@ref) records; empty for unmixed levels.
- `state.basis`: assigned quantum numbers and flavors.

Quantum labels and energy contributions are also accessible directly, e.g.
`state.n`, `state.L`, `state.J`, and `state.contact_shift_GeV`.
Use `propertynames(state)` for the full list. `fine_structure_mass_convention`
reflects mixing; `state.corrected` retains the pre-mixing convention.

## Example

Given `spec` from [`compute_spectrum`](@ref), including P levels:

```julia
state = spectrum_state(spec, "1^1P_1")
state.mass_GeV
state.corrected.mass_GeV  # before mixing
state.mixings
```

For its wavefunction, use `physical_components(spec, state)`: each component
has a signed coefficient and a radial wave. `radial_wave(spec, state.corrected)`
returns only the assigned precursor wave.

## Related

- [`spectrum_state`](@ref): select a level.
- [`physical_components`](@ref): signed components and waves.
- [`radial_wave`](@ref): unmixed or precursor wave.
"""
struct MixedState
    corrected::CorrectedState
    mixings::Vector{StateMixing}
    fine_structure_mass_convention::String
    mass_GeV::Float64
end

Base.getproperty(s::MixedState, name::Symbol) =
    hasfield(MixedState, name) ? getfield(s, name) :
    getproperty(getfield(s, :corrected), name)
Base.propertynames(::MixedState) = Tuple(union(
    fieldnames(MixedState),
    (
        fieldnames(CorrectedState)...,
        :n, :L, :multiplicity, :J, :label,
    ),
))

# Immutable update helper: append one mixing and set the mixed mass (and
# optionally the mixing-stage fine-structure convention note). The wrapped
# corrected state is provenance and never changes.
function _with_mixing(
    state::MixedState,
    mixing::StateMixing,
    mass_GeV::Real;
    fine_structure_mass_convention::AbstractString = state.fine_structure_mass_convention,
)
    return MixedState(
        getfield(state, :corrected),
        vcat(state.mixings, [mixing]),
        String(fine_structure_mass_convention),
        float(mass_GeV),
    )
end

"""
    Spectrum{S}

Staged model spectrum over one or more flavor `channels`, the `states` (in the
request order), and the underlying [`SectorComputation`](@ref). Ordinary
spectra contain one channel; the final isoscalar spectrum contains the
nonstrange and strange channels whose native solves participate in flavor
mixing. The stage is the element type `S` of `states`:

  - [`CentralSpectrum`](@ref)` = Spectrum{CentralState}` — [`central_spectrum`](@ref)
  - [`CorrectedSpectrum`](@ref)` = Spectrum{CorrectedState}` — [`fixed_spectrum`](@ref)
  - [`MixedSpectrum`](@ref)` = Spectrum{MixedState}` — [`add_intra_meson_mixing`](@ref)

[`compute_spectrum`](@ref) composes the fixed-sector solve and spectroscopic
mixing. The [`GIParameters`](@ref) that produced the spectrum live in
`computation.params`; use [`parameters`](@ref) to retrieve them.
"""
struct Spectrum{S,C<:SectorComputation}
    channels::Vector{Meson}
    states::Vector{S}
    computation::C
    function Spectrum(
        channels::AbstractVector{<:Meson},
        states::AbstractVector{S},
        computation::C,
    ) where {S,C<:SectorComputation}
        isempty(channels) && throw(ArgumentError("Spectrum needs at least one flavor channel"))
        unique_channels = collect(Meson, channels)
        length(unique(unique_channels)) == length(unique_channels) || throw(ArgumentError(
            "Spectrum flavor channels must be unique",
        ))
        owned_states = collect(S, states)
        channel_flavors = Set(
            (meson.flavor1, meson.flavor2) for meson in unique_channels
        )
        length(channel_flavors) == length(unique_channels) || throw(ArgumentError(
            "Spectrum channels must have distinct flavor identities",
        ))
        for state in owned_states
            isnothing(state.basis.flavors) && throw(ArgumentError(
                "Spectrum state $(state.label) must carry explicit flavor identity",
            ))
            state.basis.flavors in channel_flavors || throw(ArgumentError(
                "Spectrum state $(state.label) belongs to $(state.basis.flavors), " *
                "which is absent from the spectrum channels",
            ))
        end
        return new{S,C}(unique_channels, owned_states, computation)
    end
end

Spectrum(meson::Meson, states::AbstractVector, computation::SectorComputation) =
    Spectrum([meson], states, computation)

_single_channel(spec::Spectrum) = length(spec.channels) == 1 ? only(spec.channels) :
    throw(ArgumentError("operation requires a single-channel spectrum; got $(length(spec.channels)) channels"))

function _channel_for(spec::Spectrum, basis::BasisState)
    if isnothing(basis.flavors)
        return _single_channel(spec)
    end
    idx = findfirst(
        meson -> (meson.flavor1, meson.flavor2) == basis.flavors,
        spec.channels,
    )
    isnothing(idx) && throw(ArgumentError(
        "spectrum has no flavor channel $(basis.flavors) for $(basis.label)",
    ))
    return spec.channels[idx]
end

"""[`Spectrum`](@ref) after the central solve: `Spectrum{CentralState}`."""
const CentralSpectrum = Spectrum{CentralState}

"""[`Spectrum`](@ref) with spin-dependent shifts attached: `Spectrum{CorrectedState}`."""
const CorrectedSpectrum = Spectrum{CorrectedState}

"""
    MixedSpectrum

[`Spectrum`](@ref) with spectroscopic and/or flavor mixing applied:
`Spectrum{MixedState}`. Obtain one with [`compute_spectrum`](@ref) or
[`compute_isoscalar_spectrum`](@ref).

## Properties

- `spec.states`: [`MixedState`](@ref) entries in request order.
- `spec.channels`: the [`Meson`](@ref) channels represented in the spectrum.
- `spec.computation`: parameters and cached radial solutions; use [`parameters`](@ref).

Mixing belongs to individual states: use `spec.states[i].mixings`, not
`spec.mixings` or `spec.mixing`. Each entry is a [`StateMixing`](@ref).

## Example

Given `spec` from [`compute_spectrum`](@ref):

```julia
[(s.label, s.mass_GeV) for s in spec.states]
state = spectrum_state(spec, first(spec.states).basis)
state.mixings
physical_components(spec, state)
```

## Related

- [`spectrum_state`](@ref): select a level.
- [`physical_components`](@ref): signed components and waves.
- [`radial_wave`](@ref): unmixed wave.
- [`parameters`](@ref): retrieve model parameters.
"""
const MixedSpectrum = Spectrum{MixedState}

# Stages that carry a spin-resolved `mass_GeV` per state (annihilation-block inputs).
const SpinResolvedSpectrum = Union{CorrectedSpectrum,MixedSpectrum}

"""
    parameters(spec::Spectrum) -> GIParameters

The parameters used to build `spec` (single copy, stored with the cached radial
solves in [`SectorComputation`](@ref)).
"""
parameters(spec::Spectrum) = spec.computation.params

"""
    central_spectrum(params, meson; levels=spectrum_levels(2), ...) -> CentralSpectrum

Stage 1: solve the spin-independent radial problem once per distinct orbital in
`levels` and return the central eigenvalues as [`CentralState`](@ref)s together
with the filled [`SectorComputation`](@ref).

`nlevels_per_channel` bounds the radial levels kept per channel; a level with
`n` beyond it throws `ArgumentError`.

Every wave in the returned spectrum comes from `solver`. There used to be a
second, oscillator-basis wave cache here for the Table III annihilation matrix
elements; see the note above `fix_annihilation_phase!` in the mixing code for
why it is gone.
"""
function central_spectrum(
    params::GIParameters,
    meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    solver::RadialSolver = FiniteDifferenceSolver(),
)
    nlevels_per_channel = solver.nlevels_per_channel
    isempty(levels) && throw(ArgumentError("central_spectrum: empty `levels`"))
    level_keys = [(level.n, level.L_label, level.multiplicity, level.J) for level in levels]
    length(unique(level_keys)) == length(level_keys) || throw(ArgumentError(
        "central_spectrum: duplicate requested spectroscopic levels",
    ))
    masses = meson.constituent_masses
    for level in levels
        (!isnothing(level.flavors) && level.flavors != (meson.flavor1, meson.flavor2)) &&
            throw(ArgumentError(
                "central_spectrum: level $(level.label) belongs to $(level.flavors), not $(flavor_label(meson))",
            ))
        haskey(L_SYMBOLS, level.L_label) ||
            throw(ArgumentError("central_spectrum: unknown orbital label `$(level.L_label)`"))
        1 <= level.n <= nlevels_per_channel || throw(ArgumentError(
            "central_spectrum: level $(level.label) has n=$(level.n) outside 1:$nlevels_per_channel (raise nlevels_per_channel)",
        ))
    end

    # One solve per distinct orbital, with the requested solver. That is all of
    # them: nothing downstream gets a wave from anywhere else.
    observed_L = sort(unique(level.L_label for level in levels))
    channel_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    for L_label in observed_L
        requested_n = maximum(level.n for level in levels if level.L_label == L_label)
        key = RadialChannelKey(masses, L_label)
        sol = channel_solution(
            params,
            masses,
            L_SYMBOLS[L_label];
            nlevels = requested_n,
            solver = solver,
        )
        channel_cache[key] = sol
    end
    computation = SectorComputation(params, solver, channel_cache)

    states = map(collect(levels)) do level
        sol = channel_cache[RadialChannelKey(masses, level.L_label)]
        level.n <= length(sol.eigenvalues_GeV) || throw(ArgumentError(
            "central_spectrum: channel `$(level.L_label)` returned only $(length(sol.eigenvalues_GeV)) levels; requested n=$(level.n)",
        ))
        CentralState(_with_flavors(level, meson), sol.eigenvalues_GeV[level.n])
    end
    return Spectrum(meson, states, computation)
end

"""
    fixed_spectrum(params, meson; levels=spectrum_levels(2), solver, terms)

Assemble and diagonalize the complete fixed `(L,S,J)` Hamiltonian for every
requested sector and cache its native waves. This is the production entry point
before cross-sector mixing. It does not call [`central_spectrum`](@ref).

The reported central/contact/spin-orbit/tensor values are contributions in the
same solved eigenstate and therefore sum to its fixed-sector eigenvalue. With a
switch off, that operator family is absent from both the Hamiltonian and the
decomposition.
"""
function fixed_spectrum(
    params::GIParameters,
    meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    solver::RadialSolver = FiniteDifferenceSolver(),
    terms::SpinTerms = SpinTerms(),
)
    isempty(levels) && throw(ArgumentError("fixed_spectrum: empty `levels`"))
    level_keys = [(level.n, level.L_label, level.multiplicity, level.J) for level in levels]
    length(unique(level_keys)) == length(level_keys) || throw(ArgumentError(
        "fixed_spectrum: duplicate requested spectroscopic levels",
    ))
    masses = meson.constituent_masses
    for level in levels
        (!isnothing(level.flavors) && level.flavors != (meson.flavor1, meson.flavor2)) &&
            throw(ArgumentError(
                "fixed_spectrum: level $(level.label) belongs to $(level.flavors), not $(flavor_label(meson))",
            ))
        haskey(L_SYMBOLS, level.L_label) || throw(ArgumentError(
            "fixed_spectrum: unknown orbital label `$(level.L_label)`",
        ))
        1 <= level.n <= solver.nlevels_per_channel || throw(ArgumentError(
            "fixed_spectrum: level $(level.label) has n=$(level.n) outside 1:$(solver.nlevels_per_channel) (raise nlevels_per_channel)",
        ))
    end
    channel_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    required_levels = Dict{Tuple{String,Int,Int},Int}()
    for level in levels
        sector = (level.L_label, level.multiplicity, level.J)
        required_levels[sector] = max(get(required_levels, sector, 0), level.n)
    end
    states = map(collect(levels)) do level
        multiplet = FineStructureMultiplet(
            level.L_label, level.multiplicity, level.J,
        )
        key = RadialChannelKey(
            masses, level.L_label, level.multiplicity, level.J,
        )
        sol = get!(channel_cache, key) do
            fixed_channel_solution(
                params,
                masses,
                multiplet;
                solver = solver,
                terms = terms,
                nlevels = required_levels[(
                    level.L_label, level.multiplicity, level.J,
                )],
            )
        end
        wave = radial_wave(sol, level.n)
        mass = sol.eigenvalues_GeV[level.n]
        contact_shift = terms.contact_hyperfine ?
            contact_hyperfine_shift_active(params, masses, multiplet, wave) : 0.0
        so_vector = 0.0
        so_thomas = 0.0
        so_total = 0.0
        tensor = 0.0
        fs_total = 0.0
        fs_convention = "disabled"
        fine_structure_requested = terms.fine_structure && params.fine_structure.enabled &&
                                   level.L_label != "S"
        fine_structure_active = fine_structure_requested &&
                                level.L_label != "S" && level.multiplicity == 3
        if fine_structure_active
            comp = fine_structure_components(
                params,
                masses,
                multiplet,
                wave;
                enabled = true,
            )
            so_vector = comp.spin_orbit_vector
            so_thomas = comp.spin_orbit_thomas
            so_total = so_vector + so_thomas
            tensor = comp.tensor
            fs_total = comp.total
            fs_convention =
                isapprox(masses.m1_GeV, masses.m2_GeV; rtol = 0.0, atol = 0.0) ?
                "equal_mass" : "unequal_mass_equal_share_LdotS"
        elseif fine_structure_requested
            fs_convention =
                isapprox(masses.m1_GeV, masses.m2_GeV; rtol = 0.0, atol = 0.0) ?
                "equal_mass" : "unequal_mass_equal_share_LdotS"
        end
        central_contribution = mass - contact_shift - fs_total
        CorrectedState(
            _with_flavors(level, meson),
            central_contribution,
            contact_shift,
            so_vector,
            so_thomas,
            so_total,
            tensor,
            fs_total,
            fs_convention,
            mass,
        )
    end
    computation = SectorComputation(params, solver, channel_cache)
    return Spectrum(meson, states, computation)
end

"""
    add_intra_meson_mixing(spec::CorrectedSpectrum; terms=SpinTerms()) -> MixedSpectrum

Stage 3: apply intra-meson same-`J` antisymmetric spin-orbit mixing (unequal
flavor only) and triplet tensor `L = J∓1` mixing on the corrected masses.
Mixing needs the fine-structure operators, so it only acts when stage 2
actually applied them (any state with a convention other than `"disabled"`);
otherwise every state passes through with empty `mixings`.

Mixed eigenvalues are assigned order-preservingly: ascending mixed mass to
ascending unmixed diagonal. Each state's [`StateMixing`](@ref) records the full
block so other orderings can be reconstructed downstream.
"""
function add_intra_meson_mixing(
    spec::CorrectedSpectrum;
    terms::SpinTerms = SpinTerms(),
)
    same_j_spin_orbit_mixing, tensor_mixing = terms.same_j_spin_orbit, terms.tensor
    params = parameters(spec)
    meson = _single_channel(spec)
    masses = meson.constituent_masses
    channel_cache = spec.computation.channel_cache
    states = [
        MixedState(s, StateMixing[], s.fine_structure_mass_convention, s.mass_GeV) for
        s in spec.states
    ]
    fine_structure_applied =
        any(s -> s.fine_structure_mass_convention != "disabled", spec.states)
    if same_j_spin_orbit_mixing && fine_structure_applied && !is_equal_flavor(meson)
        _apply_same_j_spin_orbit_mixing!(states, params, masses, channel_cache)
    end
    if tensor_mixing && fine_structure_applied
        _apply_tensor_mixing!(states, params, masses, channel_cache)
    end
    return Spectrum(spec.channels, states, spec.computation)
end

"""
    compute_spectrum(params, meson; levels=spectrum_levels(2), ...) -> MixedSpectrum

Run the physical production stages: [`fixed_spectrum`](@ref) →
[`add_intra_meson_mixing`](@ref). The independent [`central_spectrum`](@ref)
diagnostic is not evaluated. Returns a [`MixedSpectrum`](@ref) containing
[`MixedState`](@ref) entries, including levels that did not mix.

## Example

```julia
using GIModel
path = joinpath(pkgdir(GIModel), "data", "parameters.provisional.toml")
params, mq = load_parameters_and_quark_masses(path)
meson = Meson(mq, :c, :b)
```

Compute the spectrum, then inspect one level:

```julia
spec = compute_spectrum(params, meson; levels=spectrum_levels(1))
state = spectrum_state(spec, "1^1P_1")
state.mass_GeV
physical_components(spec, state)
```

## Related

[`load_parameters_and_quark_masses`](@ref), [`Meson`](@ref), [`spectrum_levels`](@ref),
[`RadialSolver`](@ref), [`SpinTerms`](@ref), [`spectrum_state`](@ref),
[`physical_components`](@ref).
"""
function compute_spectrum(
    params::GIParameters,
    meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    solver::RadialSolver = FiniteDifferenceSolver(),
    terms::SpinTerms = SpinTerms(),
)
    corrected = fixed_spectrum(
        params, meson; levels = levels, solver = solver, terms = terms,
    )
    return add_intra_meson_mixing(corrected; terms = terms)
end

# Assign ascending block eigenvalues to block members ordered by ascending
# unmixed mass, recording full-block provenance on each member.
function _assign_block_members!(
    states::Vector{MixedState},
    member_indices::Vector{Int},
    result::MixingResult,
    ;
    fine_structure_mass_convention = nothing,
)
    length(member_indices) == length(result.block.basis) || throw(ArgumentError(
        "mixing member count does not match the result basis",
    ))
    ordered = sort(member_indices; by = i -> states[i].mass_GeV)
    ascending = result.masses
    for (rank, i) in enumerate(ordered)
        # `vectors` columns are ordered by ascending eigenvalue (see
        # diagonalize_mixing_block), matching `ascending`.
        mixing = StateMixing(
            result,
            rank,
            states[i].mass_GeV,
        )
        states[i] = if isnothing(fine_structure_mass_convention)
            _with_mixing(states[i], mixing, ascending[rank])
        else
            _with_mixing(
                states[i],
                mixing,
                ascending[rank];
                fine_structure_mass_convention = fine_structure_mass_convention,
            )
        end
    end
    return states
end

function _cached_state_wave(channel_cache, masses, state)
    fixed = RadialChannelKey(masses, state.L, state.multiplicity, state.J)
    key = haskey(channel_cache, fixed) ? fixed : RadialChannelKey(masses, state.L)
    return radial_wave(channel_cache[key], state.n)
end

function _cached_state_wave(spec::Spectrum, state)
    meson = _channel_for(spec, state.basis)
    return _cached_state_wave(
        spec.computation.channel_cache, meson.constituent_masses, state,
    )
end

function _apply_same_j_spin_orbit_mixing!(
    states::Vector{MixedState},
    params::GIParameters,
    masses::ConstituentMasses,
    channel_cache::Dict{RadialChannelKey,ChannelRadialSolution},
)
    groups = Dict{Tuple{String,Int},Vector{Int}}()
    for (i, s) in pairs(states)
        haskey(L_SYMBOLS, s.L) || continue
        Lval = L_SYMBOLS[s.L]
        (Lval > 0 && s.J == Lval && s.multiplicity in (1, 3)) || continue
        push!(get!(groups, (s.L, s.J), Int[]), i)
    end
    for indices in values(groups)
        singlets = [i for i in indices if states[i].multiplicity == 1]
        triplets = [i for i in indices if states[i].multiplicity == 3]
        (isempty(singlets) || isempty(triplets)) && continue
        ordered = vcat(sort(singlets; by = i -> states[i].n),
                       sort(triplets; by = i -> states[i].n))
        basis = [s.basis for s in states[ordered]]
        matrix = Matrix(Diagonal([states[i].mass_GeV for i in ordered]))
        for ia in eachindex(singlets), ib in eachindex(triplets)
            i, j = singlets[ia], triplets[ib]
            left = _cached_state_wave(channel_cache, masses, states[i])
            right = _cached_state_wave(channel_cache, masses, states[j])
            element = spin_orbit_mixing_components(
                params, masses, states[i].L, left, right;
                enabled = true,
            ).total
            row = findfirst(==(i), ordered)
            col = findfirst(==(j), ordered)
            matrix[row, col] = matrix[col, row] = element
        end
        block = MixingBlock(
            "same-J antisymmetric spin-orbit",
            basis,
            matrix;
            mechanism = "antisymmetric_spin_orbit",
            source = "Godfrey-Isgur Eqs. (6)-(7)",
        )
        result = diagonalize_mixing_block(block)
        _assign_block_members!(
            states,
            ordered,
            result,
            fine_structure_mass_convention = "unequal_mass_same_j_mixed",
        )
    end
    return states
end

function _apply_tensor_mixing!(
    states::Vector{MixedState},
    params::GIParameters,
    masses::ConstituentMasses,
    channel_cache::Dict{RadialChannelKey,ChannelRadialSolution},
)
    groups = Dict{Int,Vector{Int}}()
    for (i, s) in pairs(states)
        s.multiplicity == 3 || continue
        haskey(L_SYMBOLS, s.L) || continue
        s.J > 0 || continue
        L = L_SYMBOLS[s.L]
        (L == s.J - 1 || L == s.J + 1) || continue
        push!(get!(groups, s.J, Int[]), i)
    end
    for indices in values(groups)
        low_indices = [i for i in indices if L_SYMBOLS[states[i].L] == states[i].J - 1]
        high_indices = [i for i in indices if L_SYMBOLS[states[i].L] == states[i].J + 1]
        (isempty(low_indices) || isempty(high_indices)) && continue
        ordered = vcat(sort(low_indices; by = i -> states[i].n),
                       sort(high_indices; by = i -> states[i].n))
        basis = [s.basis for s in states[ordered]]
        matrix = Matrix(Diagonal([states[i].mass_GeV for i in ordered]))
        for i in low_indices, j in high_indices
            low_radial = _cached_state_wave(channel_cache, masses, states[i])
            high_radial = _cached_state_wave(channel_cache, masses, states[j])
            element = tensor_mixing_components(
                params, masses, low_radial, high_radial, states[i].J;
                enabled = true,
            ).total
            row = findfirst(==(i), ordered)
            col = findfirst(==(j), ordered)
            matrix[row, col] = matrix[col, row] = element
        end
        block = MixingBlock(
            "same-J tensor triplet L/L'",
            basis,
            matrix;
            mechanism = "tensor_mixing",
        )
        mix = diagonalize_mixing_block(block)
        _assign_block_members!(
            states,
            ordered,
            mix,
        )
    end
    return states
end

"""
    spectrum_state(spectrum, n, L_label, multiplicity, J)
    spectrum_state(spectrum, level::BasisState)
    spectrum_state(spectrum, label::AbstractString)

Look up one state by quantum numbers at any stage; throws `ArgumentError` when
absent. The return type is the spectrum's state type.

The string form takes the state's own `label` (`"1^3S_1"`), which is the handle
every state carries and every report prints. If several flavor channels share
the quantum numbers or label, use a [`BasisState`](@ref) with explicit `flavors`;
an ambiguous lookup throws `ArgumentError`.

## Example

Given `spec` from [`compute_spectrum`](@ref):

```julia
state = spectrum_state(spec, "1^1P_1")
state.mass_GeV
physical_components(spec, state)
```

## Related

[`MixedState`](@ref), [`MixedSpectrum`](@ref), [`radial_wave`](@ref),
[`physical_components`](@ref).
"""
function spectrum_state(
    spec::Spectrum,
    n::Integer,
    L_label::AbstractString,
    multiplicity::Integer,
    J::Integer,
)
    indices = findall(
        s ->
            s.n == n &&
            s.L == String(L_label) &&
            s.multiplicity == multiplicity &&
            s.J == J,
        spec.states,
    )
    isempty(indices) && throw(ArgumentError(
        "spectrum has no state $(n)^$(multiplicity)$(L_label)_$(J)",
    ))
    length(indices) == 1 || throw(ArgumentError(
        "state $(n)^$(multiplicity)$(L_label)_$(J) is flavor-ambiguous; pass a BasisState with explicit `flavors`",
    ))
    return spec.states[only(indices)]
end

function spectrum_state(spec::Spectrum, level::BasisState)
    indices = findall(
        s ->
            s.n == level.n &&
            s.L == level.L_label &&
            s.multiplicity == level.multiplicity &&
            s.J == level.J &&
            (isnothing(level.flavors) || s.basis.flavors == level.flavors),
        spec.states,
    )
    isempty(indices) && throw(ArgumentError(
        "spectrum has no state matching $(level.label) with flavors=$(level.flavors)",
    ))
    length(indices) == 1 || throw(ArgumentError(
        "state $(level.label) is flavor-ambiguous; set `flavors` on BasisState",
    ))
    return spec.states[only(indices)]
end

function spectrum_state(spec::Spectrum, label::AbstractString)
    indices = findall(s -> s.label == label, spec.states)
    isempty(indices) && throw(ArgumentError(
        "spectrum has no state `$label`; " *
        "available: $(join((s.label for s in spec.states), ", "))",
    ))
    length(indices) == 1 || throw(ArgumentError(
        "state label `$label` is flavor-ambiguous; use spectrum_state(spec, BasisState(...; flavors=(...)))",
    ))
    return spec.states[only(indices)]
end

"""
    radial_wave(spec, state) -> RadialWave
    radial_wave(spec, label) -> RadialWave

The native radial wavefunction behind an unmixed level of `spec`. For a
spin-resolved spectrum this is the eigenvector of the complete fixed-sector
Hamiltonian, not the central precursor.

This is the hand-off from the spectrum to every wavefunction-level observable:
annihilation and leptonic widths, two-photon amplitudes, charge radii and the
radiative transition moments all take a `RadialWave`. The solve is
already cached on the spectrum, so this is a lookup, not a recomputation.

For a [`MixedState`](@ref) with nonempty `mixings`, this method throws
`ArgumentError`: use [`physical_components`](@ref) to obtain all signed
components. `radial_wave(spec, state.corrected)` retrieves only the assigned
pre-mixing basis wave, not the physical superposition.

## Example

Given `spec` from [`compute_spectrum`](@ref):

```julia
state = spectrum_state(spec, "1^1S_0")
wave = radial_wave(spec, state)  # this level is unmixed in this example
wave_norm(wave)
```

## Related

[`RadialWave`](@ref), [`sample_wave`](@ref), [`wave_norm`](@ref),
[`physical_components`](@ref), [`spectrum_state`](@ref).

Throws `ArgumentError` when the orbital was not in `levels` or when `n` exceeds
the levels kept per channel.
"""
function radial_wave(spec::Spectrum, state::Union{CentralState,CorrectedState,MixedState})
    if state isa MixedState && !isempty(state.mixings)
        throw(ArgumentError(
            "$(state.label) is a mixed physical state; use physical_components " *
            "instead of selecting one precursor radial wave",
        ))
    end
    return _cached_state_wave(spec, state)
end

"""
    physical_components(spec, state_or_label)

Resolve a physical state's fully flattened spectroscopic/flavor composition.
Each returned named tuple contains the pure `basis` identity, its signed
`coefficient`, and that component's native radial `wave`. Unmixed states return
one unit component. Sequential spin and flavor transformations are composed
from their shared [`StateMixing`](@ref) eigensystems, and paths ending at the
same native basis state are combined.

## Example

Given `spec` from [`compute_spectrum`](@ref):

```julia
components = physical_components(spec, "1^1P_1")
[(c.basis.flavors, c.basis.label, c.coefficient) for c in components]
[wave_norm(c.wave) for c in components]
sum(abs2(c.coefficient) for c in components)  # approximately 1
```

Coefficients are signed amplitudes; squaring them gives component weights.
Keep their signs when computing interference-sensitive observables.

## Related

[`MixedState`](@ref), [`spectrum_state`](@ref), [`radial_wave`](@ref),
[`RadialWave`](@ref), [`sample_wave`](@ref), [`physical_state_amplitude`](@ref).
"""
function physical_components(
    spec::Spectrum,
    state::Union{CentralState,CorrectedState,MixedState},
)
    if !(state isa MixedState) || isempty(state.mixings)
        basis = state.basis
        return [(basis = basis, coefficient = 1.0,
                 wave = _cached_state_wave(spec, state))]
    end
    mixing = last(state.mixings)
    result = mixing.result
    coefficients = @view result.vectors[:, mixing.eigenstate]
    components = NamedTuple[]
    for (i, basis) in pairs(result.block.basis)
        component_state = spectrum_state(spec, basis)
        source_indices = findall(
            source -> source.result === result,
            component_state.mixings,
        )
        length(source_indices) == 1 || throw(ArgumentError(
            "$(basis.label) does not carry the shared $(result.block.name) result",
        ))
        source_index = only(source_indices)
        source_mixing = component_state.mixings[source_index]
        prior_mixings = component_state.mixings[1:(source_index - 1)]
        prior_state = MixedState(
            component_state.corrected,
            prior_mixings,
            component_state.fine_structure_mass_convention,
            source_mixing.unmixed_GeV,
        )
        for component in physical_components(spec, prior_state)
            push!(components, (
                basis = component.basis,
                coefficient = coefficients[i] * component.coefficient,
                wave = component.wave,
            ))
        end
    end
    merged = NamedTuple[]
    positions = Dict{BasisState,Int}()
    for component in components
        position = get(positions, component.basis, 0)
        if position == 0
            push!(merged, component)
            positions[component.basis] = length(merged)
        else
            prior = merged[position]
            merged[position] = (
                basis = prior.basis,
                coefficient = prior.coefficient + component.coefficient,
                wave = prior.wave,
            )
        end
    end
    return merged
end

physical_components(spec::Spectrum, label::AbstractString) =
    physical_components(spec, spectrum_state(spec, label))

"""
    physical_state_amplitude(kernel, spec, state_or_label)

Coherently apply a linear amplitude `kernel(component)` to every signed native
component of one physical state. This is the observable-side counterpart of
[`physical_components`](@ref): consumers do not need to read a mixing matrix or
reconstruct a parallel flavor/radial coefficient vector.

`kernel` receives one named tuple with `basis`, `coefficient`, and `wave`
fields. Its return value must exclude the mixing coefficient; this function
multiplies by that coefficient exactly once.
"""
function physical_state_amplitude(
    kernel,
    spec::Spectrum,
    state::Union{CentralState,CorrectedState,MixedState},
)
    return sum(
        component.coefficient * kernel(component) for
        component in physical_components(spec, state)
    )
end

physical_state_amplitude(kernel, spec::Spectrum, label::AbstractString) =
    physical_state_amplitude(kernel, spec, spectrum_state(spec, label))

"""
    physical_transition_amplitude(kernel, spec, left, right)

Coherently apply a bilinear transition `kernel(left_component,
right_component)` between two physical states. The signed coefficients of both
states are composed automatically. The kernel owns the observable's selection
rules: it should return zero for flavor, spin, or angular components that the
operator does not connect.
"""
function physical_transition_amplitude(
    kernel,
    spec::Spectrum,
    left::Union{CentralState,CorrectedState,MixedState},
    right::Union{CentralState,CorrectedState,MixedState},
)
    left_components = physical_components(spec, left)
    right_components = physical_components(spec, right)
    return sum(
        a.coefficient * b.coefficient * kernel(a, b) for
        a in left_components for b in right_components
    )
end

physical_transition_amplitude(
    kernel,
    spec::Spectrum,
    left::AbstractString,
    right::AbstractString,
) = physical_transition_amplitude(
    kernel, spec, spectrum_state(spec, left), spectrum_state(spec, right),
)

function radial_expect(
    spec::Spectrum,
    state::Union{CentralState,CorrectedState,MixedState},
    f,
)
    components = physical_components(spec, state)
    # A scalar radial operator is diagonal in L,S,J, but not in radial n.
    # Retain interference between radial levels in the same spectroscopic
    # sector; angular/spin-orthogonal sectors contribute incoherently.
    return sum(
        left.coefficient * right.coefficient *
        radial_overlap(left.wave, right.wave, f) for
        left in components for right in components if
        left.basis.L_label == right.basis.L_label &&
        left.basis.multiplicity == right.basis.multiplicity &&
        left.basis.J == right.basis.J &&
        left.basis.flavors == right.basis.flavors
    )
end

function radial_wave(spec::Spectrum, label::AbstractString)
    state = spectrum_state(spec, label)
    return radial_wave(spec, state)
end

# --- Display -----------------------------------------------------------------
# A spectrum carries its radial eigenvectors, so the default struct dump is
# ~200 kB of nested constructors. These methods print the physics instead: one
# row per level, one column per contribution, so the mass decomposition is
# readable at the REPL without knowing any field names.

_stage_name(::Type{CentralState}) = "CentralSpectrum"
_stage_name(::Type{CorrectedState}) = "CorrectedSpectrum"
_stage_name(::Type{MixedState}) = "MixedSpectrum"

_stage_columns(::Type{CentralState}) = ("central",)
_stage_columns(::Type{CorrectedState}) = ("central", "contact", "fine str", "mass")
_stage_columns(::Type{MixedState}) = ("central", "contact", "fine str", "mixing", "mass")

_stage_values(s::CentralState) = (s.central_GeV,)
_stage_values(s::CorrectedState) =
    (s.central_GeV, s.contact_shift_GeV, s.fine_structure_shift_GeV, s.mass_GeV)
_stage_values(s::MixedState) = (
    s.central_GeV,
    s.contact_shift_GeV,
    s.fine_structure_shift_GeV,
    s.mass_GeV - s.corrected.mass_GeV,
    s.mass_GeV,
)

function Base.show(io::IO, ::MIME"text/plain", spec::Spectrum{S}) where {S}
    cols = _stage_columns(S)
    channel_text = join((flavor_label(m) for m in spec.channels), " + ")
    println(
        io,
        _stage_name(S), ": ", channel_text, ", ",
        length(spec.states), " levels — all values in GeV",
    )
    width = isempty(spec.states) ? 8 : maximum(length(s.label) for s in spec.states)
    print(io, "  ", rpad("level", width))
    for c in cols
        print(io, lpad(c, 11))
    end
    for s in spec.states
        print(io, "\n  ", rpad(s.label, width))
        for v in _stage_values(s)
            print(io, lpad(@sprintf("%.4f", v), 11))
        end
    end
    if S !== CentralState
        nmix = count(s -> !isempty(s.mixings), spec.states)
        nmix > 0 && print(io, "\n  (", nmix, " levels carry mixing; see `spec.states[i].mixings`)")
    end
    return nothing
end

Base.show(io::IO, spec::Spectrum{S}) where {S} = print(
    io, _stage_name(S), "(", join((flavor_label(m) for m in spec.channels), "+"), ", ",
    length(spec.states), " levels)",
)
