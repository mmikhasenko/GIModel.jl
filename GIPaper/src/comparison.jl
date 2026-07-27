# Model-vs-paper comparison: group reference rows by Meson, compute one
# GIModel.Spectrum per meson, match states back to rows, and emit the
# comparison NamedTuple rows consumed by write_residual_report.
#
# Public API (exported from GIPaper.jl): compare_reference

const L_SYMBOLS = GIModel.L_SYMBOLS

# The base comparison row wants the pre-mixing values, which are exactly the
# wrapped corrected-stage state (the same-J assignment below re-applies the
# mixed convention note under reference ordering).
_pre_mixing_GeV(s::MixedState) = s.corrected.mass_GeV

_pre_mixing_convention(s::MixedState) = s.corrected.fine_structure_mass_convention

function _mixing_of(s::MixedState, mechanism::AbstractString)
    idx = findfirst(m -> m.mechanism == mechanism, s.mixings)
    return isnothing(idx) ? nothing : s.mixings[idx]
end

function _base_row(row::ReferenceState, meson::Meson, s::MixedState)
    predicted = _pre_mixing_GeV(s)
    return (
        sector = row.sector,
        state = row.composition,
        L = row.L,
        n = row.n,
        J = row.J,
        multiplicity = row.multiplicity,
        m1_GeV = meson.constituent_masses.m1_GeV,
        m2_GeV = meson.constituent_masses.m2_GeV,
        fine_structure_mass_convention = _pre_mixing_convention(s),
        reference_GeV = row.mass_GeV,
        central_GeV = s.central_GeV,
        contact_shift_GeV = s.contact_shift_GeV,
        spin_orbit_vector_shift_GeV = s.spin_orbit_vector_shift_GeV,
        spin_orbit_thomas_shift_GeV = s.spin_orbit_thomas_shift_GeV,
        spin_orbit_shift_GeV = s.spin_orbit_shift_GeV,
        tensor_shift_GeV = s.tensor_shift_GeV,
        fine_structure_shift_GeV = s.fine_structure_shift_GeV,
        annihilation_shift_GeV = 0.0,
        annihilation_scheme = "none",
        same_j_mixing_scheme = "none",
        same_j_unmixed_GeV = predicted,
        same_j_offdiag_GeV = 0.0,
        same_j_mixing_angle_deg = 0.0,
        same_j_component_singlet = NaN,
        same_j_component_triplet = NaN,
        tensor_mixing_scheme = "none",
        tensor_unmixed_GeV = predicted,
        tensor_offdiag_GeV = 0.0,
        tensor_component_lowL = NaN,
        tensor_component_highL = NaN,
        isoscalar_annihilation_scheme = "none",
        isoscalar_annihilation_unmixed_GeV = predicted,
        predicted_GeV = predicted,
        residual_MeV = 1000 * (predicted - row.mass_GeV),
        confidence = row.confidence,
    )
end

# Reconstruct the 2x2 eigenvector matrix of a mixing block from its two member
# states: the member with the lower final mass carries the lower-eigenvalue
# column (compute_spectrum assigns ascending mixed mass to ascending unmixed
# diagonal, and each StateMixing stores its own column).
function _block_columns(s1::MixedState, s2::MixedState, mechanism::AbstractString)
    m1 = _mixing_of(s1, mechanism)
    m2 = _mixing_of(s2, mechanism)
    (isnothing(m1) || isnothing(m2)) && return nothing
    by_model = s1.mass_GeV <= s2.mass_GeV ? (m1, m2) : (m2, m1)
    vectors = hcat(by_model[1].components, by_model[2].components)
    return (masses = m1.partner_masses_GeV, vectors = vectors, mixing = m1)
end

# Two rows sharing a 2x2 mixing block: assign ascending mixed masses either by
# reference-mass ordering (paper convention, byte-compatible with the old
# comparison) or by ascending unmixed model prediction (:model_order).
function _pair_order(rows, i1::Int, i2::Int, mixed_assignment::Symbol)
    if mixed_assignment == :reference_order
        return sort([i1, i2]; by = i -> rows[i].reference_GeV)
    elseif mixed_assignment == :model_order
        return sort([i1, i2]; by = i -> rows[i].predicted_GeV)
    end
    throw(ArgumentError("unsupported mixed_assignment `$mixed_assignment`"))
end

function _assign_same_j_rows!(
    rows::Vector{NamedTuple},
    group_rows::Vector{Int},
    states::Vector{MixedState},
    mixed_assignment::Symbol,
)
    pairs_by_key = Dict{Tuple{Int,String,Int},Vector{Int}}()
    for i in group_rows
        isnothing(_mixing_of(states[i], "antisymmetric_spin_orbit")) && continue
        push!(get!(pairs_by_key, (rows[i].n, rows[i].L, rows[i].J), Int[]), i)
    end
    for indices in values(pairs_by_key)
        length(indices) == 2 || continue
        isinglet = findfirst(i -> rows[i].multiplicity == 1, indices)
        itriplet = findfirst(i -> rows[i].multiplicity == 3, indices)
        (isnothing(isinglet) || isnothing(itriplet)) && continue
        i1 = indices[isinglet]
        i2 = indices[itriplet]
        block = _block_columns(states[i1], states[i2], "antisymmetric_spin_orbit")
        isnothing(block) && continue
        ordered = _pair_order(rows, i1, i2, mixed_assignment)
        for (level, irow) in enumerate(ordered)
            row = rows[irow]
            predicted = block.masses[level]
            vec = block.vectors[:, level]
            rows[irow] = merge(
                row,
                (
                    same_j_mixing_scheme = "antisymmetric_spin_orbit",
                    same_j_unmixed_GeV = row.predicted_GeV,
                    same_j_offdiag_GeV = block.mixing.offdiag_GeV,
                    same_j_mixing_angle_deg = block.mixing.mixing_angle_deg,
                    same_j_component_singlet = vec[1],
                    same_j_component_triplet = vec[2],
                    predicted_GeV = predicted,
                    residual_MeV = 1000 * (predicted - row.reference_GeV),
                    fine_structure_mass_convention = "unequal_mass_same_j_mixed",
                ),
            )
        end
    end
    return rows
end

function _assign_tensor_rows!(
    rows::Vector{NamedTuple},
    group_rows::Vector{Int},
    states::Vector{MixedState},
    mixed_assignment::Symbol,
)
    pairs_by_key = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in group_rows
        isnothing(_mixing_of(states[i], "tensor_mixing")) && continue
        L = L_SYMBOLS[rows[i].L]
        partner_level = L == rows[i].J - 1 ? rows[i].n - 1 : rows[i].n
        push!(get!(pairs_by_key, (rows[i].J, partner_level), Int[]), i)
    end
    for indices in values(pairs_by_key)
        length(indices) == 2 || continue
        ilow = findfirst(i -> L_SYMBOLS[rows[i].L] == rows[i].J - 1, indices)
        ihigh = findfirst(i -> L_SYMBOLS[rows[i].L] == rows[i].J + 1, indices)
        (isnothing(ilow) || isnothing(ihigh)) && continue
        i1 = indices[ilow]
        i2 = indices[ihigh]
        block = _block_columns(states[i1], states[i2], "tensor_mixing")
        isnothing(block) && continue
        ordered = _pair_order(rows, i1, i2, mixed_assignment)
        for (level, irow) in enumerate(ordered)
            row = rows[irow]
            predicted = block.masses[level]
            vec = block.vectors[:, level]
            rows[irow] = merge(
                row,
                (
                    tensor_mixing_scheme = "tensor_mixing",
                    tensor_unmixed_GeV = row.predicted_GeV,
                    tensor_offdiag_GeV = block.mixing.offdiag_GeV,
                    tensor_component_lowL = vec[1],
                    tensor_component_highL = vec[2],
                    predicted_GeV = predicted,
                    residual_MeV = 1000 * (predicted - row.reference_GeV),
                ),
            )
        end
    end
    return rows
end

# --- isoscalar annihilation schemes -----------------------------------------
#
# Flavor-mixed rows share quantum numbers (eta vs eta-prime), so eigenvalues
# are always assigned by reference-mass ordering here regardless of
# `mixed_assignment` — the reference identity is what distinguishes the rows.

function _assign_pseudoscalar_rows!(
    rows::Vector{NamedTuple},
    scheme::Symbol,
    params::GIParameters,
    nn_spec::Spectrum,
    ss_spec::Spectrum,
)
    matches = [
        (i, row) for (i, row) in pairs(rows) if
        row.sector == "isoscalar" &&
        row.L == "S" &&
        row.multiplicity == 1 &&
        row.J == 0 &&
        row.n in (1, 2)
    ]
    length(matches) == 4 || return rows
    seen_n = Set(row.n for (_, row) in matches)
    (1 in seen_n && 2 in seen_n) || return rows

    solution = if scheme == :calibrated_p1
        targets = sort([row.reference_GeV for (_, row) in matches])
        pseudoscalar_annihilation_block(
            CalibratedP1Annihilation(),
            params,
            nn_spec,
            ss_spec;
            targets = targets,
        )
    elseif scheme in (:p1, :paper_p1)
        pseudoscalar_annihilation_block(PaperP1Annihilation(), params, nn_spec, ss_spec)
    else
        pseudoscalar_annihilation_block(PaperP2Annihilation(), params, nn_spec, ss_spec)
    end

    ordered = sort(matches; by = item -> item[2].reference_GeV)
    for (level, (i, row)) in enumerate(ordered)
        predicted = solution.masses[level]
        rows[i] = merge(
            row,
            (
                annihilation_shift_GeV = predicted - row.predicted_GeV,
                annihilation_scheme = String(scheme),
                isoscalar_annihilation_scheme = String(scheme),
                isoscalar_annihilation_unmixed_GeV = row.predicted_GeV,
                predicted_GeV = predicted,
                residual_MeV = 1000 * (predicted - row.reference_GeV),
            ),
        )
    end
    return rows
end

function _assign_general_s1_rows!(
    rows::Vector{NamedTuple},
    params::GIParameters,
    nn_spec::Spectrum,
    ss_spec::Spectrum,
)
    matches = [
        (i, row) for (i, row) in pairs(rows) if
        row.sector == "isoscalar" &&
        row.L == "S" &&
        row.multiplicity == 3 &&
        row.J == 1 &&
        row.n == 1
    ]
    length(matches) == 2 || return rows
    solution = isoscalar_annihilation_block(
        params,
        nn_spec,
        ss_spec,
        BasisState(1, "S", 3, 1);
        amplitude_A = params.annihilation.s1_A,
    )
    ordered = sort(matches; by = item -> item[2].reference_GeV)
    for (level, (i, row)) in enumerate(ordered)
        predicted = solution.masses[level]
        rows[i] = merge(
            row,
            (
                annihilation_shift_GeV = predicted - row.predicted_GeV,
                annihilation_scheme = "general_s1",
                isoscalar_annihilation_scheme = "general_s1",
                isoscalar_annihilation_unmixed_GeV = row.predicted_GeV,
                predicted_GeV = predicted,
                residual_MeV = 1000 * (predicted - row.reference_GeV),
            ),
        )
    end
    return rows
end

function _assign_table_iii_rows!(
    rows::Vector{NamedTuple},
    params::GIParameters,
    nn_spec::Spectrum,
    ss_spec::Spectrum,
)
    groups = Dict{Tuple{Int,Int,String,Int},Vector{Int}}()
    for (i, row) in pairs(rows)
        row.sector == "isoscalar" || continue
        # Pseudoscalar rows were already assigned by the P1 block.
        row.isoscalar_annihilation_scheme == "none" || continue
        push!(get!(groups, (row.n, row.multiplicity, row.L, row.J), Int[]), i)
    end
    for ((n, multiplicity, L_label, J), indices) in groups
        length(indices) == 2 || continue
        ordered = sort(indices; by = i -> rows[i].reference_GeV)
        ns_row = rows[ordered[1]]
        ns_pred = ns_row.predicted_GeV
        level = BasisState(n, L_label, multiplicity, J)
        ss = spectrum_state(ss_spec, level)
        ss_pred = ss.mass_GeV

        amplitude = table_iii_amplitude(params, L_label, multiplicity, J)
        if isnothing(amplitude)
            predicted = [ns_pred, ss_pred]
            scheme = "ideal"
            components = [(1.0, 0.0), (0.0, 1.0)]
        else
            solution = isoscalar_annihilation_block(
                params,
                nn_spec,
                ss_spec,
                level;
                amplitude_A = amplitude,
            )
            predicted = solution.masses
            scheme = "general_eq16"
            components = [
                (solution.vectors[1, 1], solution.vectors[2, 1]),
                (solution.vectors[1, 2], solution.vectors[2, 2]),
            ]
        end

        for (lvl, irow) in enumerate(ordered)
            row = rows[irow]
            unmixed = lvl == 1 ? ns_pred : ss_pred
            # The heavier row is dominantly `s sbar`: rewrite its fixed-sector
            # breakdown columns from the strange channel so reports stay truthful.
            diag_update = lvl == 2 ?
                (
                    m1_GeV = ss_spec.meson.constituent_masses.m1_GeV,
                    m2_GeV = ss_spec.meson.constituent_masses.m2_GeV,
                    central_GeV = ss.central_GeV,
                    contact_shift_GeV = ss.contact_shift_GeV,
                    spin_orbit_vector_shift_GeV = ss.spin_orbit_vector_shift_GeV,
                    spin_orbit_thomas_shift_GeV = ss.spin_orbit_thomas_shift_GeV,
                    spin_orbit_shift_GeV = ss.spin_orbit_shift_GeV,
                    tensor_shift_GeV = ss.tensor_shift_GeV,
                    fine_structure_shift_GeV = ss.fine_structure_shift_GeV,
                ) : NamedTuple()
            rows[irow] = merge(
                row,
                diag_update,
                (
                    annihilation_shift_GeV = predicted[lvl] - unmixed,
                    annihilation_scheme = scheme,
                    isoscalar_annihilation_scheme = scheme,
                    isoscalar_annihilation_unmixed_GeV = unmixed,
                    isoscalar_component_ns = components[lvl][1],
                    isoscalar_component_ss = components[lvl][2],
                    predicted_GeV = predicted[lvl],
                    residual_MeV = 1000 * (predicted[lvl] - row.reference_GeV),
                ),
            )
        end
    end
    return rows
end

"""
    compare_reference(params, quark_masses, reference; ...) -> Vector{NamedTuple}

Compare the model to a vector of [`ReferenceState`](@ref) rows:

 1. group rows by [`reference_meson`](@ref) and run one
    [`GIModel.compute_spectrum`](@ref) per meson (the reference rows define
    exactly which levels are computed);
 2. match each row to its model state by `(n, multiplicity, L, J)`;
 3. reassign 2×2 mixing-block eigenvalues per `mixed_assignment`:
    `:reference_order` (default) reproduces the paper convention of assigning
    ascending mixed masses by ascending reference mass; `:model_order` keeps
    the model's own assignment (ascending unmixed prediction);
 4. optionally apply the isoscalar annihilation `scheme`
    (`:calibrated_p1`, `:p1`/`:paper_p1`, `:p2`/`:paper_p2`, `:general_s1`,
    `:p1_and_s1`, `:table_iii`) — requires `strange_mass_GeV`.

Rows with `n > 6` are skipped with a warning. The output rows feed
[`write_residual_report`](@ref) / [`nonmixing_deviation_summary`](@ref).
"""
function compare_reference(
    params::GIParameters,
    quark_masses::QuarkMassTable,
    reference::AbstractVector{ReferenceState};
    contact_hyperfine::Bool = true,
    use_fine_structure::Bool = true,
    antisymmetric_spin_orbit_mixing::Bool = true,
    tensor_mixing::Bool = true,
    isoscalar_pseudoscalar_annihilation::Symbol = :none,
    strange_mass_GeV = nothing,
    mixed_assignment::Symbol = :reference_order,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    eigensolver::Symbol = :full,
    annihilation_wave_basis::Symbol = :ho,
    ho_wave_L::Tuple{Vararg{String}} = ("S",),
)
    scheme = isoscalar_pseudoscalar_annihilation
    scheme in (:none, :calibrated_p1, :p1, :paper_p1, :p2, :paper_p2, :general_s1, :p1_and_s1, :table_iii) ||
        throw(ArgumentError("unsupported isoscalar annihilation scheme `$scheme`"))
    mixed_assignment in (:reference_order, :model_order) ||
        throw(ArgumentError("unsupported mixed_assignment `$mixed_assignment`"))
    # The Table III `^3P_2` block needs HO-basis P waves.
    if scheme == :table_iii && !("P" in ho_wave_L)
        ho_wave_L = (ho_wave_L..., "P")
    end

    kept = ReferenceState[]
    for row in reference
        if row.n > 6
            @warn "skipping reference row with n > 6" sector = row.sector state = row.composition n = row.n
            continue
        end
        push!(kept, row)
    end
    isempty(kept) && return NamedTuple[]

    # Group rows per meson and compute one spectrum each; keep the per-row
    # state association (duplicate quantum numbers, e.g. eta/eta-prime flavor
    # pairs, map to identical states).
    mesons = Meson[]
    group_of = Dict{Meson,Int}()
    group_rows = Vector{Vector{Int}}()
    meson_of = Vector{Int}(undef, length(kept))
    for (i, row) in pairs(kept)
        meson = reference_meson(quark_masses, row)
        g = get!(group_of, meson) do
            push!(mesons, meson)
            push!(group_rows, Int[])
            length(mesons)
        end
        push!(group_rows[g], i)
        meson_of[i] = g
    end
    specs = Vector{MixedSpectrum}(undef, length(mesons))
    states = Vector{MixedState}(undef, length(kept))
    for g in eachindex(mesons)
        levels = [
            BasisState(kept[i].n, kept[i].L, kept[i].multiplicity, kept[i].J) for
            i in group_rows[g]
        ]
        specs[g] = compute_spectrum(
            params,
            mesons[g];
            levels = levels,
            ngrid = ngrid,
            rmax = rmax,
            kinetic = kinetic,
            eigensolver = eigensolver,
            contact_hyperfine = contact_hyperfine,
            use_fine_structure = use_fine_structure,
            same_j_spin_orbit_mixing = antisymmetric_spin_orbit_mixing,
            tensor_mixing = tensor_mixing,
            annihilation_wave_basis = annihilation_wave_basis,
            ho_wave_L = ho_wave_L,
        )
        for (k, i) in enumerate(group_rows[g])
            states[i] = specs[g].states[k]
        end
    end

    rows = NamedTuple[
        _base_row(kept[i], mesons[meson_of[i]], states[i]) for i in eachindex(kept)
    ]

    for g in eachindex(mesons)
        _assign_same_j_rows!(rows, group_rows[g], states, mixed_assignment)
        _assign_tensor_rows!(rows, group_rows[g], states, mixed_assignment)
    end

    if scheme != :none
        isnothing(strange_mass_GeV) && throw(ArgumentError(
            "strange_mass_GeV is required for isoscalar pseudoscalar annihilation",
        ))
        iso = findfirst(m -> any(kept[i].sector == "isoscalar" for i in group_rows[group_of[m]]), mesons)
        if !isnothing(iso)
            nn_spec = specs[iso]
            strange_meson = Meson(StrangeQuark(strange_mass_GeV), StrangeQuark(strange_mass_GeV))
            seen = Set{Tuple{Int,String,Int,Int}}()
            ss_levels = BasisState[]
            for i in group_rows[iso]
                key = (kept[i].n, kept[i].L, kept[i].multiplicity, kept[i].J)
                key in seen && continue
                push!(seen, key)
                push!(ss_levels, BasisState(key...))
            end
            # The strange partner diagonals mirror the old comparison: unmixed
            # fixed-sector predictions (no intra-meson mixing blocks).
            ss_spec = compute_spectrum(
                params,
                strange_meson;
                levels = ss_levels,
                ngrid = ngrid,
                rmax = rmax,
                kinetic = kinetic,
                eigensolver = eigensolver,
                contact_hyperfine = contact_hyperfine,
                use_fine_structure = use_fine_structure,
                same_j_spin_orbit_mixing = false,
                tensor_mixing = false,
                annihilation_wave_basis = annihilation_wave_basis,
                ho_wave_L = ho_wave_L,
            )
            if scheme == :general_s1
                _assign_general_s1_rows!(rows, params, nn_spec, ss_spec)
            elseif scheme == :p1_and_s1
                _assign_pseudoscalar_rows!(rows, :calibrated_p1, params, nn_spec, ss_spec)
                _assign_general_s1_rows!(rows, params, nn_spec, ss_spec)
            elseif scheme == :table_iii
                _assign_pseudoscalar_rows!(rows, :calibrated_p1, params, nn_spec, ss_spec)
                _assign_table_iii_rows!(rows, params, nn_spec, ss_spec)
            else
                _assign_pseudoscalar_rows!(rows, scheme, params, nn_spec, ss_spec)
            end
        end
    end

    return rows
end
