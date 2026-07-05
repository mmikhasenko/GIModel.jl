# Pure model spectrum for one meson: explicit level specification, cached radial
# solves, contribution breakdown, and intra-meson (same-J spin-orbit, tensor)
# mixing applied on the model's own eigenvalues — no reference data involved.
#
# Public API (exported from GIModel.jl):
#   spectrum_levels, StateMixing, SpectrumState, Spectrum, compute_spectrum,
#   spectrum_state, parameters

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

Provenance of one intra-meson mixing applied to a [`SpectrumState`](@ref):
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
    SpectrumState

One ``n\\,^{2S+1}L_J`` model state: contribution breakdown relative to the
central eigenvalue (contact hyperfine, spin-orbit vector/Thomas, tensor), the
mixings applied (empty for unmixed states), and the final `mass_GeV`.
"""
struct SpectrumState
    n::Int
    L::String
    multiplicity::Int
    J::Int
    label::String
    central_GeV::Float64
    contact_shift_GeV::Float64
    spin_orbit_vector_shift_GeV::Float64
    spin_orbit_thomas_shift_GeV::Float64
    spin_orbit_shift_GeV::Float64
    tensor_shift_GeV::Float64
    fine_structure_shift_GeV::Float64
    fine_structure_mass_convention::String
    mixings::Vector{StateMixing}
    mass_GeV::Float64
end

# Immutable update helper: append one mixing and set the mixed mass (and
# optionally the fine-structure convention note).
function _with_mixing(
    state::SpectrumState,
    mixing::StateMixing,
    mass_GeV::Real;
    fine_structure_mass_convention::AbstractString = state.fine_structure_mass_convention,
)
    return SpectrumState(
        state.n,
        state.L,
        state.multiplicity,
        state.J,
        state.label,
        state.central_GeV,
        state.contact_shift_GeV,
        state.spin_orbit_vector_shift_GeV,
        state.spin_orbit_thomas_shift_GeV,
        state.spin_orbit_shift_GeV,
        state.tensor_shift_GeV,
        state.fine_structure_shift_GeV,
        String(fine_structure_mass_convention),
        vcat(state.mixings, [mixing]),
        float(mass_GeV),
    )
end

"""
    Spectrum

Output of [`compute_spectrum`](@ref): the input `meson`, the computed `states`
(in `levels` order), and the underlying [`SectorComputation`](@ref) (kept so
two-meson flavor mixing can reuse the cached radial solves; see
`flavor_mixing.jl`). The [`GIParameters`](@ref) that produced the spectrum live
in `computation.params`; use [`parameters`](@ref) to retrieve them.
"""
struct Spectrum
    meson::Meson
    states::Vector{SpectrumState}
    computation::SectorComputation
end

"""
    parameters(spec::Spectrum) -> GIParameters

The parameters used to build `spec` (single copy, stored with the cached radial
solves in [`SectorComputation`](@ref)).
"""
parameters(spec::Spectrum) = spec.computation.params

"""
    compute_spectrum(params, meson; levels=spectrum_levels(2), ...) -> Spectrum

Solve the spin-independent radial problem once per distinct orbital in `levels`,
attach contact-hyperfine and fine-structure shifts per state, apply intra-meson
same-`J` antisymmetric spin-orbit mixing (unequal flavor only) and triplet
tensor `L = J∓1` mixing on the model's own eigenvalues, and return a
[`Spectrum`](@ref).

Mixed eigenvalues are assigned order-preservingly: ascending mixed mass to
ascending unmixed diagonal. Each state's [`StateMixing`](@ref) records the full
block so other orderings can be reconstructed downstream.

`nlevels_per_channel` bounds the radial levels kept per channel; a level with
`n` beyond it throws `ArgumentError`. The HO wave cache (`annihilation_wave_basis
= :ho`, orbitals in `ho_wave_L`) stores phase-fixed harmonic-oscillator waves for
annihilation matrix elements, exactly as the paper path requires.
"""
function compute_spectrum(
    params::GIParameters,
    meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
    nlevels_per_channel::Integer = 6,
    contact_hyperfine::Bool = true,
    use_fine_structure::Bool = params.fine_structure.enabled,
    same_j_spin_orbit_mixing::Bool = true,
    tensor_mixing::Bool = true,
    annihilation_wave_basis::Symbol = :ho,
    ho_wave_L::Tuple{Vararg{String}} = ("S",),
)
    isempty(levels) && throw(ArgumentError("compute_spectrum: empty `levels`"))
    masses = meson.constituent_masses
    for level in levels
        haskey(L_SYMBOLS, level.L_label) ||
            throw(ArgumentError("compute_spectrum: unknown orbital label `$(level.L_label)`"))
        1 <= level.n <= nlevels_per_channel || throw(ArgumentError(
            "compute_spectrum: level $(level.label) has n=$(level.n) outside 1:$nlevels_per_channel (raise nlevels_per_channel)",
        ))
    end

    # One FD solve per distinct orbital; HO waves for the annihilation phase
    # convention (see fix_annihilation_phase!) where requested.
    observed_L = sort(unique(level.L_label for level in levels))
    channel_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    for L_label in observed_L
        key = RadialChannelKey(masses, L_label)
        ev, vecs, r = channel_solution(
            params,
            masses,
            L_SYMBOLS[L_label];
            nlevels = nlevels_per_channel,
            ngrid = ngrid,
            rmax = rmax,
            kinetic = kinetic,
            eigensolver = eigensolver,
        )
        channel_cache[key] = ChannelRadialSolution(ev, vecs, r)
    end
    ho_wave_cache = Dict{RadialChannelKey,ChannelRadialSolution}()
    if annihilation_wave_basis == :ho
        ho_params = with_basis(params, HarmonicOscillatorBasis)
        for L_label in observed_L
            L_label in ho_wave_L || continue
            key = RadialChannelKey(masses, L_label)
            ev, vecs, r = channel_solution(
                ho_params,
                masses,
                L_SYMBOLS[L_label];
                nlevels = 2,
                ngrid = ngrid,
                rmax = rmax,
                kinetic = kinetic,
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
    computation = SectorComputation(params, channel_cache, ho_wave_cache)

    # Fixed-sector breakdown per requested level.
    contact_level_cache = Dict{Tuple{RadialChannelKey,Int},Vector{Float64}}()
    states = SpectrumState[]
    for level in levels
        key = RadialChannelKey(masses, level.L_label)
        sol = channel_cache[key]
        level.n <= length(sol.eigenvalues_GeV) || throw(ArgumentError(
            "compute_spectrum: channel `$(level.L_label)` returned only $(length(sol.eigenvalues_GeV)) levels; requested n=$(level.n)",
        ))
        central = sol.eigenvalues_GeV[level.n]
        multiplet = FineStructureMultiplet(level.L_label, level.multiplicity, level.J)
        wave = RadialWaveOnUniformMesh(sol, level.n)
        contact_shift = 0.0
        if contact_hyperfine
            nonperturbative = get!(contact_level_cache, (key, level.multiplicity)) do
                contact_hyperfine_nonperturbative_levels(
                    params,
                    masses,
                    level.L_label,
                    level.multiplicity,
                    sol.r,
                    length(sol.eigenvalues_GeV),
                )
            end
            if !isempty(nonperturbative) && level.n <= length(nonperturbative)
                contact_shift = nonperturbative[level.n] - central
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
        if use_fine_structure && params.fine_structure.enabled
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
        push!(
            states,
            SpectrumState(
                level.n,
                level.L_label,
                level.multiplicity,
                level.J,
                level.label,
                central,
                contact_shift,
                so_vector,
                so_thomas,
                so_total,
                tensor,
                fs_total,
                fs_convention,
                StateMixing[],
                central + contact_shift + fs_total,
            ),
        )
    end

    fine_structure_active = use_fine_structure && params.fine_structure.enabled
    if same_j_spin_orbit_mixing && fine_structure_active && !is_equal_flavor(meson)
        _apply_same_j_spin_orbit_mixing!(states, params, masses, channel_cache)
    end
    if tensor_mixing && fine_structure_active
        _apply_tensor_mixing!(states, params, masses, channel_cache)
    end

    return Spectrum(meson, states, computation)
end

# Assign ascending block eigenvalues to block members ordered by ascending
# unmixed mass, recording full-block provenance on each member.
function _assign_block_members!(
    states::Vector{SpectrumState},
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
    states::Vector{SpectrumState},
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
    states::Vector{SpectrumState},
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
    spectrum_state(spectrum, n, L_label, multiplicity, J) -> SpectrumState
    spectrum_state(spectrum, level::BasisState) -> SpectrumState

Look up one state by quantum numbers; throws `ArgumentError` when absent.
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
