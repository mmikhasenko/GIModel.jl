# Pure model spectrum for one meson, staged: (1) spin-independent central
# solve, (2) complete fixed-sector spin Hamiltonians, (3) intra-meson
# (same-J spin-orbit, tensor) mixing on the model's own eigenvalues — no
# reference data involved. Each stage is a first-class value; the stage is
# carried by the element type of `Spectrum.states`, not by a tag.
#
# Public API (exported from GIModel.jl):
#   spectrum_levels, StateMixing, CentralState, CorrectedState, MixedState,
#   Spectrum, CentralSpectrum, CorrectedSpectrum, MixedSpectrum,
#   central_spectrum, add_spin_corrections, add_intra_meson_mixing,
#   compute_spectrum, spectrum_state, parameters

"""
    spectrum_levels(nmax; L_labels=("S", "P", "D")) -> Vector{BasisState}

Enumerate the ``n\\,^{2S+1}L_J`` multiplets for `n in 1:nmax` and each orbital
letter: the spin singlet (`J = L`) and the spin triplets (`J = 1` for `S` waves,
`J in L-1:L+1` otherwise). Use as the `levels` argument of [`compute_spectrum`](@ref);
hand-built [`BasisState`](@ref) vectors work the same way for custom selections.
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

Provenance of one intra-meson mixing applied to a [`MixedState`](@ref):
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
"""
struct CentralState
    n::Int
    L::String
    multiplicity::Int
    J::Int
    label::String
    central_GeV::Float64
end

"""
    CorrectedState

A [`CentralState`](@ref) plus a reporting decomposition of the eigenvalue from
the complete fixed-sector Hamiltonian (contact, vector/Thomas spin-orbit, and
tensor) and its `mass_GeV`. These fields are not perturbative inputs. Element type of
[`CorrectedSpectrum`](@ref). Properties of the wrapped `central` state
(`n`, `L`, `central_GeV`, …) forward transparently.
"""
struct CorrectedState
    central::CentralState
    contact_shift_GeV::Float64
    spin_orbit_vector_shift_GeV::Float64
    spin_orbit_thomas_shift_GeV::Float64
    spin_orbit_shift_GeV::Float64
    tensor_shift_GeV::Float64
    fine_structure_shift_GeV::Float64
    fine_structure_mass_convention::String
    mass_GeV::Float64
end

Base.getproperty(s::CorrectedState, name::Symbol) =
    hasfield(CorrectedState, name) ? getfield(s, name) :
    getproperty(getfield(s, :central), name)
Base.propertynames(::CorrectedState) =
    Tuple(union(fieldnames(CorrectedState), fieldnames(CentralState)))

"""
    MixedState

A [`CorrectedState`](@ref) plus the intra-meson `mixings` applied to it (empty
for unmixed states) and the final `mass_GeV`. `fine_structure_mass_convention`
is the mixing-stage view: same-J mixing overrides it to
`"unequal_mass_same_j_mixed"` while the wrapped `corrected` state keeps the
pre-mixing convention. Element type of [`MixedSpectrum`](@ref). Properties of
the wrapped stages forward transparently.
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
Base.propertynames(::MixedState) = Tuple(
    union(fieldnames(MixedState), fieldnames(CorrectedState), fieldnames(CentralState)),
)

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

Staged model spectrum for one meson: the `meson`, the `states` (in the request
`levels` order), and the underlying [`SectorComputation`](@ref) (kept so later
stages and two-meson flavor mixing reuse the cached radial solves; see
`flavor_mixing.jl`). The stage is the element type `S` of `states`:

  - [`CentralSpectrum`](@ref)` = Spectrum{CentralState}` — [`central_spectrum`](@ref)
  - [`CorrectedSpectrum`](@ref)` = Spectrum{CorrectedState}` — [`add_spin_corrections`](@ref)
  - [`MixedSpectrum`](@ref)` = Spectrum{MixedState}` — [`add_intra_meson_mixing`](@ref)

[`compute_spectrum`](@ref) composes the three stages. The [`GIParameters`](@ref)
that produced the spectrum live in `computation.params`; use [`parameters`](@ref)
to retrieve them.
"""
struct Spectrum{S,M<:Meson,C<:SectorComputation}
    meson::M
    states::Vector{S}
    computation::C
end

"""[`Spectrum`](@ref) after the central solve: `Spectrum{CentralState}`."""
const CentralSpectrum = Spectrum{CentralState}

"""[`Spectrum`](@ref) with spin-dependent shifts attached: `Spectrum{CorrectedState}`."""
const CorrectedSpectrum = Spectrum{CorrectedState}

"""[`Spectrum`](@ref) with intra-meson mixing applied: `Spectrum{MixedState}`."""
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
    masses = meson.constituent_masses
    for level in levels
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
        CentralState(
            level.n,
            level.L_label,
            level.multiplicity,
            level.J,
            level.label,
            sol.eigenvalues_GeV[level.n],
        )
    end
    return Spectrum(meson, states, computation)
end

"""
    add_spin_corrections(spec::CentralSpectrum; terms=SpinTerms()) -> CorrectedSpectrum

Stage 2: assemble and diagonalize the complete fixed `(L,S,J)` Hamiltonian for
every requested sector, cache its native waves, and report the displacement
from the central eigenvalue by operator family. With a switch off, that family
is absent from the Hamiltonian.
"""
function add_spin_corrections(
    spec::CentralSpectrum;
    terms::SpinTerms = SpinTerms(),
)
    params = parameters(spec)
    masses = spec.meson.constituent_masses
    # A spin-resolved solve is physics-dependent. Never write it into the
    # central spectrum's cache: the same central value may legitimately feed a
    # second correction stage with different SpinTerms.
    channel_cache = copy(spec.computation.channel_cache)
    required_levels = Dict{Tuple{String,Int,Int},Int}()
    for state in spec.states
        sector = (state.L, state.multiplicity, state.J)
        required_levels[sector] = max(get(required_levels, sector, 0), state.n)
    end
    states = map(spec.states) do state
        multiplet = FineStructureMultiplet(state.L, state.multiplicity, state.J)
        key = RadialChannelKey(
            masses, state.L, state.multiplicity, state.J,
        )
        sol = get!(channel_cache, key) do
            fixed_channel_solution(
                params,
                masses,
                multiplet;
                solver = spec.computation.solver,
                terms = terms,
                nlevels = required_levels[(state.L, state.multiplicity, state.J)],
            )
        end
        wave = radial_wave(sol, state.n)
        mass = sol.eigenvalues_GeV[state.n]
        total_shift = mass - state.central_GeV
        is_contact_sector = terms.contact_hyperfine && state.L == "S" &&
                            state.multiplicity in (1, 3)
        contact_shift = is_contact_sector ? total_shift : 0.0
        so_vector = 0.0
        so_thomas = 0.0
        so_total = 0.0
        tensor = 0.0
        fs_total = 0.0
        fs_convention = "disabled"
        fine_structure_requested = terms.fine_structure && params.fine_structure.enabled &&
                                   state.L != "S"
        fine_structure_active = fine_structure_requested &&
                                state.L != "S" && state.multiplicity == 3
        if fine_structure_active
            comp = fine_structure_components(
                params,
                masses,
                multiplet,
                wave;
                enabled = true,
                k_spin_orbit = params.fine_structure.k_spin_orbit,
                k_tensor = params.fine_structure.k_tensor,
            )
            # Full diagonalization changes the radial state, so the exact mass
            # displacement is not merely the first-order expectation. Preserve
            # the operator decomposition by distributing the exact displacement
            # according to the component expectations on the solved wave.
            scale = iszero(comp.total) ? 0.0 : total_shift / comp.total
            so_vector = scale * comp.spin_orbit_vector
            so_thomas = scale * comp.spin_orbit_thomas
            so_total = so_vector + so_thomas
            tensor = scale * comp.tensor
            fs_total = total_shift
            fs_convention =
                isapprox(masses.m1_GeV, masses.m2_GeV; rtol = 0.0, atol = 0.0) ?
                "equal_mass" : "unequal_mass_equal_share_LdotS"
        elseif fine_structure_requested
            fs_convention =
                isapprox(masses.m1_GeV, masses.m2_GeV; rtol = 0.0, atol = 0.0) ?
                "equal_mass" : "unequal_mass_equal_share_LdotS"
        end
        CorrectedState(
            state,
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
    computation = SectorComputation(params, spec.computation.solver, channel_cache)
    return Spectrum(spec.meson, states, computation)
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
    masses = spec.meson.constituent_masses
    channel_cache = spec.computation.channel_cache
    states = [
        MixedState(s, StateMixing[], s.fine_structure_mass_convention, s.mass_GeV) for
        s in spec.states
    ]
    fine_structure_applied =
        any(s -> s.fine_structure_mass_convention != "disabled", spec.states)
    if same_j_spin_orbit_mixing && fine_structure_applied && !is_equal_flavor(spec.meson)
        _apply_same_j_spin_orbit_mixing!(states, params, masses, channel_cache)
    end
    if tensor_mixing && fine_structure_applied
        _apply_tensor_mixing!(states, params, masses, channel_cache)
    end
    return Spectrum(spec.meson, states, spec.computation)
end

"""
    compute_spectrum(params, meson; levels=spectrum_levels(2), ...) -> MixedSpectrum

Run all three stages: [`central_spectrum`](@ref) →
[`add_spin_corrections`](@ref) → [`add_intra_meson_mixing`](@ref). Keyword
arguments are forwarded to the stage they belong to; call the stages directly
to inspect the intermediate [`CentralSpectrum`](@ref) / [`CorrectedSpectrum`](@ref).
"""
function compute_spectrum(
    params::GIParameters,
    meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    solver::RadialSolver = FiniteDifferenceSolver(),
    terms::SpinTerms = SpinTerms(),
)
    central = central_spectrum(params, meson; levels = levels, solver = solver)
    corrected = add_spin_corrections(central; terms = terms)
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
        basis = [
            BasisState(s.n, s.L, s.multiplicity, s.J; label = s.label) for
            s in states[ordered]
        ]
        matrix = Matrix(Diagonal([states[i].mass_GeV for i in ordered]))
        for ia in eachindex(singlets), ib in eachindex(triplets)
            i, j = singlets[ia], triplets[ib]
            left = _cached_state_wave(channel_cache, masses, states[i])
            right = _cached_state_wave(channel_cache, masses, states[j])
            element = spin_orbit_mixing_components(
                params, masses, states[i].L, left, right;
                enabled = true,
                k_spin_orbit = params.fine_structure.k_spin_orbit,
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
        basis = [
            BasisState(s.n, s.L, 3, s.J; label = s.label) for s in states[ordered]
        ]
        matrix = Matrix(Diagonal([states[i].mass_GeV for i in ordered]))
        for i in low_indices, j in high_indices
            low_radial = _cached_state_wave(channel_cache, masses, states[i])
            high_radial = _cached_state_wave(channel_cache, masses, states[j])
            element = tensor_mixing_components(
                params, masses, low_radial, high_radial, states[i].J;
                enabled = true,
                k_tensor = params.fine_structure.k_tensor,
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
every state carries and every report prints.
"""
function spectrum_state(
    spec::Spectrum,
    n::Integer,
    L_label::AbstractString,
    multiplicity::Integer,
    J::Integer,
)
    idx = findfirst(
        s ->
            s.n == n &&
            s.L == String(L_label) &&
            s.multiplicity == multiplicity &&
            s.J == J,
        spec.states,
    )
    isnothing(idx) && throw(ArgumentError(
        "spectrum has no state $(n)^$(multiplicity)$(L_label)_$(J) for meson $(flavor_label(spec.meson))",
    ))
    return spec.states[idx]
end

spectrum_state(spec::Spectrum, level::BasisState) =
    spectrum_state(spec, level.n, level.L_label, level.multiplicity, level.J)

function spectrum_state(spec::Spectrum, label::AbstractString)
    idx = findfirst(s -> s.label == label, spec.states)
    isnothing(idx) && throw(ArgumentError(
        "spectrum has no state `$label` for meson $(flavor_label(spec.meson)); " *
        "available: $(join((s.label for s in spec.states), ", "))",
    ))
    return spec.states[idx]
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

There is one wave per level, produced by the spectrum's own solver — the same
wavefunction the eigenvalues came from. A `wave_basis` keyword used to select a
separately-cached oscillator wave here; that is gone, and the reason it existed
is recorded at `annihilation_basis_input`.

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
    masses = spec.meson.constituent_masses
    cache = spec.computation.channel_cache
    return _cached_state_wave(cache, masses, state)
end

"""
    physical_components(spec, state_or_label)

Resolve a physical state's spectroscopic composition. Each returned named tuple
contains the pure `basis` label, its signed `coefficient`, and that component's
native radial `wave`. Unmixed states return one unit component; mixed states use
the exact eigenvector stored in their existing [`StateMixing`](@ref).
"""
function physical_components(
    spec::Spectrum,
    state::Union{CentralState,CorrectedState,MixedState},
)
    if !(state isa MixedState) || isempty(state.mixings)
        basis = BasisState(
            state.n, state.L, state.multiplicity, state.J; label = state.label,
        )
        return [(basis = basis, coefficient = 1.0,
                 wave = _cached_state_wave(
                     spec.computation.channel_cache,
                     spec.meson.constituent_masses,
                     state,
                 ))]
    end
    length(state.mixings) == 1 || throw(ArgumentError(
        "physical composition for sequential mixing is not implemented; " *
        "assemble the mechanisms in one block",
    ))
    mixing = only(state.mixings)
    result = mixing.result
    coefficients = @view result.vectors[:, mixing.eigenstate]
    return [
        begin
            component_state = spectrum_state(spec, basis.label)
            (
                basis = basis,
                coefficient = coefficients[i],
                wave = _cached_state_wave(
                    spec.computation.channel_cache,
                    spec.meson.constituent_masses,
                    component_state,
                ),
            )
        end for (i, basis) in pairs(result.block.basis)
    ]
end

physical_components(spec::Spectrum, label::AbstractString) =
    physical_components(spec, spectrum_state(spec, label))

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
        left.basis.J == right.basis.J
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
    m = spec.meson
    cols = _stage_columns(S)
    println(
        io,
        _stage_name(S), ": ", m.flavor1, " ", m.flavor2, "bar  (m = ",
        m.constituent_masses.m1_GeV, ", ", m.constituent_masses.m2_GeV, " GeV), ",
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
        nmix > 0 && print(io, "\n  (", nmix, " levels carry intra-meson mixing; see `.mixings`)")
    end
    return nothing
end

Base.show(io::IO, spec::Spectrum{S}) where {S} = print(
    io, _stage_name(S), "(", flavor_label(spec.meson), ", ",
    length(spec.states), " levels)",
)
