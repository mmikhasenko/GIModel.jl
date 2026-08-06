# Pure model spectrum for one meson, staged: (1) spin-independent central
# solve, (2) contact-hyperfine + fine-structure shifts, (3) intra-meson
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
which mechanism, the block basis labels, this state's eigenvector `components`
(in block-basis order), and `partner_masses_GeV` — **all** block eigenvalues in
ascending order, so downstream consumers can reassign eigenvalues under a
different ordering convention without re-diagonalizing.
"""
struct StateMixing
    mechanism::String
    block_label::String
    partner_labels::Vector{String}
    components::Vector{Float64}
    offdiag_GeV::Float64
    mixing_angle_deg::Float64
    unmixed_GeV::Float64
    partner_masses_GeV::Vector{Float64}
end

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

A [`CentralState`](@ref) plus the spin-dependent first-order shifts (contact
hyperfine, spin-orbit vector/Thomas, tensor) and the corrected
`mass_GeV = central + contact + fine structure`. Element type of
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
struct Spectrum{S}
    meson::Meson
    states::Vector{S}
    computation::SectorComputation
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
`n` beyond it throws `ArgumentError`. The HO wave cache (`annihilation_wave_basis
= :ho`, orbitals in `ho_wave_L`) stores phase-fixed harmonic-oscillator waves for
annihilation matrix elements, exactly as the paper path requires.
"""
function central_spectrum(
    params::GIParameters,
    meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    solver::RadialSolver = FiniteDifferenceSolver(),
    ngrid = nothing,
    rmax = nothing,
    kinetic = nothing,
    eigensolver = nothing,
    nlevels_per_channel = nothing,
    annihilation_wave_basis::Symbol = :ho,
    ho_wave_L::Tuple{Vararg{String}} = ("S",),
)
    solver = _solver_with_legacy(solver, "central_spectrum"; ngrid = ngrid, rmax = rmax,
        kinetic = kinetic, eigensolver = eigensolver,
        nlevels_per_channel = nlevels_per_channel)
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

    # One solve per distinct orbital with the requested solver; HO waves for the
    # annihilation phase convention (see fix_annihilation_phase!) where requested.
    observed_L = sort(unique(level.L_label for level in levels))
    channel_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    for L_label in observed_L
        key = RadialChannelKey(masses, L_label)
        ev, vecs, r = channel_solution(
            params,
            masses,
            L_SYMBOLS[L_label];
            nlevels = nlevels_per_channel,
            solver = solver,
        )
        channel_cache[key] = ChannelRadialSolution(ev, vecs, r)
    end
    ho_wave_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    if annihilation_wave_basis == :ho
        # Deliberately the oscillator path regardless of `solver`: the Table III
        # annihilation amplitudes are defined on the paper's own basis. It borrows
        # only the reporting mesh, so these waves land on the same `r` as the
        # central ones and can be compared with them.
        ho_solver = solver isa OscillatorSolver ? solver :
                    OscillatorSolver(ngrid = solver.ngrid, rmax = solver.rmax)
        for L_label in observed_L
            L_label in ho_wave_L || continue
            key = RadialChannelKey(masses, L_label)
            ev, vecs, r = channel_solution(
                params,
                masses,
                L_SYMBOLS[L_label];
                nlevels = 2,
                solver = ho_solver,
            )
            # The eigensolver returns arbitrary-sign columns; without a fixed
            # phase the S_L smearing factor can flip sign between quark masses
            # and radial levels, randomizing off-diagonal annihilation matrix
            # elements. The GI Table III convention is Φ(0) > 0 in momentum
            # space (see fix_annihilation_phase!).
            fix_annihilation_phase!(vecs, r)
            ho_wave_cache[key] = ChannelRadialSolution(ev, vecs, r)
        end
    end
    computation = SectorComputation(params, solver, channel_cache, ho_wave_cache)

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
    add_spin_corrections(spec::CentralSpectrum; contact_hyperfine=true,
                         use_fine_structure=parameters(spec).fine_structure.enabled)
        -> CorrectedSpectrum

Stage 2: attach the contact-hyperfine and fine-structure shifts per state on
the cached radial waves and set each state's corrected mass. With a switch
off, the corresponding shifts are zero (and the fine-structure convention is
`"disabled"`), so the corrected mass falls back to the central eigenvalue.
"""
function add_spin_corrections(
    spec::CentralSpectrum;
    terms::SpinTerms = SpinTerms(),
    contact_hyperfine = nothing,
    use_fine_structure = nothing,
)
    terms = _terms_with_legacy(terms, "add_spin_corrections";
        contact_hyperfine = contact_hyperfine, use_fine_structure = use_fine_structure)
    contact_hyperfine, use_fine_structure = terms.contact_hyperfine, terms.fine_structure
    params = parameters(spec)
    masses = spec.meson.constituent_masses
    channel_cache = spec.computation.channel_cache
    fine_structure_active = use_fine_structure && params.fine_structure.enabled
    contact_level_cache = Dict{Tuple{RadialChannelKey,Int},Vector{Float64}}()
    states = map(spec.states) do state
        key = RadialChannelKey(masses, state.L)
        sol = channel_cache[key]
        multiplet = FineStructureMultiplet(state.L, state.multiplicity, state.J)
        wave = RadialWaveOnUniformMesh(sol, state.n)
        contact_shift = 0.0
        if contact_hyperfine
            nonperturbative = get!(contact_level_cache, (key, state.multiplicity)) do
                contact_hyperfine_nonperturbative_levels(
                    params,
                    masses,
                    state.L,
                    state.multiplicity,
                    sol.r,
                    length(sol.eigenvalues_GeV);
                    solver = spec.computation.solver,
                )
            end
            if !isempty(nonperturbative) && state.n <= length(nonperturbative)
                contact_shift = nonperturbative[state.n] - state.central_GeV
            else
                contact_shift = contact_hyperfine_shift_active(params, masses, multiplet, wave)
            end
        end
        so_vector = 0.0
        so_thomas = 0.0
        so_total = 0.0
        tensor = 0.0
        fs_total = 0.0
        fs_convention = "disabled"
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
            so_vector = comp.spin_orbit_vector
            so_thomas = comp.spin_orbit_thomas
            so_total = comp.spin_orbit
            tensor = comp.tensor
            fs_total = comp.total
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
            state.central_GeV + contact_shift + fs_total,
        )
    end
    return Spectrum(spec.meson, states, spec.computation)
end

"""
    add_intra_meson_mixing(spec::CorrectedSpectrum;
                           same_j_spin_orbit_mixing=true, tensor_mixing=true)
        -> MixedSpectrum

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
    same_j_spin_orbit_mixing = nothing,
    tensor_mixing = nothing,
)
    terms = _terms_with_legacy(terms, "add_intra_meson_mixing";
        same_j_spin_orbit_mixing = same_j_spin_orbit_mixing, tensor_mixing = tensor_mixing)
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
    ngrid = nothing,
    rmax = nothing,
    kinetic = nothing,
    eigensolver = nothing,
    nlevels_per_channel = nothing,
    contact_hyperfine = nothing,
    use_fine_structure = nothing,
    same_j_spin_orbit_mixing = nothing,
    tensor_mixing = nothing,
    annihilation_wave_basis::Symbol = :ho,
    ho_wave_L::Tuple{Vararg{String}} = ("S",),
)
    # Fold any deprecated keywords in once, here, then hand the stages the
    # objects — so a plain `compute_spectrum(params, meson)` never trips the
    # deprecation path on its way down.
    solver = _solver_with_legacy(solver, "compute_spectrum"; ngrid = ngrid, rmax = rmax,
        kinetic = kinetic, eigensolver = eigensolver,
        nlevels_per_channel = nlevels_per_channel)
    terms = _terms_with_legacy(terms, "compute_spectrum";
        contact_hyperfine = contact_hyperfine, use_fine_structure = use_fine_structure,
        same_j_spin_orbit_mixing = same_j_spin_orbit_mixing, tensor_mixing = tensor_mixing)
    central = central_spectrum(
        params,
        meson;
        levels = levels,
        solver = solver,
        annihilation_wave_basis = annihilation_wave_basis,
        ho_wave_L = ho_wave_L,
    )
    corrected = add_spin_corrections(central; terms = terms)
    return add_intra_meson_mixing(corrected; terms = terms)
end

# Assign ascending block eigenvalues to block members ordered by ascending
# unmixed mass, recording full-block provenance on each member.
function _assign_block_members!(
    states::Vector{MixedState},
    member_indices::Vector{Int},
    mechanism::AbstractString,
    block_label::AbstractString,
    basis_labels::Vector{String},
    mixed_masses::Vector{Float64},
    vectors::Matrix{Float64},
    offdiag_GeV::Real,
    mixing_angle_deg::Real;
    fine_structure_mass_convention = nothing,
)
    ordered = sort(member_indices; by = i -> states[i].mass_GeV)
    ascending = sort(mixed_masses)
    for (rank, i) in enumerate(ordered)
        # `vectors` columns are ordered by ascending eigenvalue (see
        # diagonalize_mixing_block), matching `ascending`.
        mixing = StateMixing(
            String(mechanism),
            String(block_label),
            basis_labels,
            collect(vectors[:, rank]),
            float(offdiag_GeV),
            float(mixing_angle_deg),
            states[i].mass_GeV,
            ascending,
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

function _apply_same_j_spin_orbit_mixing!(
    states::Vector{MixedState},
    params::GIParameters,
    masses::ConstituentMasses,
    channel_cache::Dict{RadialChannelKey,ChannelRadialSolution},
)
    groups = Dict{Tuple{Int,String,Int},Vector{Int}}()
    for (i, s) in pairs(states)
        haskey(L_SYMBOLS, s.L) || continue
        Lval = L_SYMBOLS[s.L]
        (Lval > 0 && s.J == Lval && s.multiplicity in (1, 3)) || continue
        push!(get!(groups, (s.n, s.L, s.J), Int[]), i)
    end
    for indices in values(groups)
        length(indices) == 2 || continue
        isinglet = findfirst(i -> states[i].multiplicity == 1, indices)
        itriplet = findfirst(i -> states[i].multiplicity == 3, indices)
        (isnothing(isinglet) || isnothing(itriplet)) && continue
        singlet = states[indices[isinglet]]
        triplet = states[indices[itriplet]]
        key = RadialChannelKey(masses, singlet.L)
        sol = channel_cache[key]
        radial = RadialWaveOnUniformMesh(sol, singlet.n)
        offdiag = spin_orbit_mixing_components(
            params,
            masses,
            singlet.L,
            radial;
            enabled = true,
            k_spin_orbit = params.fine_structure.k_spin_orbit,
        )
        mix = same_j_mixing(singlet.mass_GeV, triplet.mass_GeV, offdiag.total)
        _assign_block_members!(
            states,
            [indices[isinglet], indices[itriplet]],
            "antisymmetric_spin_orbit",
            mix.block.name,
            [singlet.label, triplet.label],
            collect(mix.masses),
            Matrix(mix.vectors),
            offdiag.total,
            mix.theta_deg;
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
    # Pair n ^3(J-1)_J with (n-1) ^3(J+1)_J, e.g. 2^3S_1 with 1^3D_1.
    groups = Dict{Tuple{Int,Int},Vector{Int}}()
    for (i, s) in pairs(states)
        s.multiplicity == 3 || continue
        haskey(L_SYMBOLS, s.L) || continue
        s.J > 0 || continue
        L = L_SYMBOLS[s.L]
        (L == s.J - 1 || L == s.J + 1) || continue
        partner_level = L == s.J - 1 ? s.n - 1 : s.n
        partner_level >= 1 || continue
        push!(get!(groups, (s.J, partner_level), Int[]), i)
    end
    for indices in values(groups)
        length(indices) == 2 || continue
        ilow = findfirst(i -> L_SYMBOLS[states[i].L] == states[i].J - 1, indices)
        ihigh = findfirst(i -> L_SYMBOLS[states[i].L] == states[i].J + 1, indices)
        (isnothing(ilow) || isnothing(ihigh)) && continue
        low = states[indices[ilow]]
        high = states[indices[ihigh]]
        low_sol = channel_cache[RadialChannelKey(masses, low.L)]
        high_sol = channel_cache[RadialChannelKey(masses, high.L)]
        low_radial = RadialWaveOnUniformMesh(low_sol, low.n)
        high_radial = RadialWaveOnUniformMesh(high_sol, high.n)
        offdiag = tensor_mixing_components(
            params,
            masses,
            low_radial,
            high_radial,
            low.J;
            enabled = true,
            k_tensor = params.fine_structure.k_tensor,
        )
        basis = [
            BasisState(low.n, low.L, 3, low.J),
            BasisState(high.n, high.L, 3, high.J),
        ]
        block = MixingBlock(
            "same-J tensor triplet L/L'",
            basis,
            [low.mass_GeV offdiag.total; offdiag.total high.mass_GeV];
            mechanism = "tensor_mixing",
        )
        mix = diagonalize_mixing_block(block)
        _assign_block_members!(
            states,
            [indices[ilow], indices[ihigh]],
            "tensor_mixing",
            block.name,
            [low.label, high.label],
            collect(mix.masses),
            Matrix(mix.vectors),
            offdiag.total,
            NaN,
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
    radial_wave(spec, label; wave_basis=:fd) -> RadialWaveOnUniformMesh
    radial_wave(spec, L_label, n; wave_basis=:fd) -> RadialWaveOnUniformMesh

The radial wavefunction behind a level of `spec` — `radial_wave(spec, "1^3S_1")`
is the wave whose eigenvalue is `spectrum_state(spec, "1^3S_1").central_GeV`.

This is the hand-off from the spectrum to every wavefunction-level observable:
annihilation and leptonic widths, two-photon amplitudes, charge radii and the
radiative transition moments all take a `RadialWaveOnUniformMesh`. The solve is
already cached on the spectrum, so this is a lookup, not a recomputation.

`wave_basis` selects which cached solve to read: `:fd` is the model's own
finite-difference wave (the default — this is the wavefunction the eigenvalues
came from), while `:ho` reads the phase-fixed harmonic-oscillator wave that the
paper's annihilation prescription requires, falling back to `:fd` when no HO
wave was cached for that channel (see `central_spectrum`'s `ho_wave_L`).

Throws `ArgumentError` when the orbital was not in `levels` or when `n` exceeds
the levels kept per channel.
"""
function radial_wave(
    spec::Spectrum,
    L_label::AbstractString,
    n::Integer;
    wave_basis::Symbol = :fd,
)
    wave_basis in (:fd, :ho) || throw(ArgumentError(
        "radial_wave: `wave_basis` must be :fd or :ho, got `$wave_basis`",
    ))
    key = RadialChannelKey(spec.meson.constituent_masses, String(L_label))
    fd = spec.computation.channel_cache
    haskey(fd, key) || throw(ArgumentError(
        "spectrum for $(flavor_label(spec.meson)) has no `$L_label` channel; include it in `levels`",
    ))
    sol = wave_basis === :ho ?
          get(spec.computation.ho_wave_cache, key, fd[key]) : fd[key]
    n <= size(sol.eigenvectors, 2) || throw(ArgumentError(
        "cached `$L_label` waves for $(flavor_label(spec.meson)) hold $(size(sol.eigenvectors, 2)) levels; requested n=$n",
    ))
    return RadialWaveOnUniformMesh(sol, n)
end

function radial_wave(spec::Spectrum, label::AbstractString; wave_basis::Symbol = :fd)
    state = spectrum_state(spec, label)
    return radial_wave(spec, state.L, state.n; wave_basis = wave_basis)
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
