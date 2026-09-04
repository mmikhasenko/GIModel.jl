# Two-meson flavor mixing bridge: build annihilation blocks from Spectrum
# objects, reusing their cached radial solves. The Eq. (16)-(18) numeric kernels
# live in pseudoscalar_annihilation.jl; the paper's channel-amplitude
# prescription (which explicit model/control and amplitudes a report requests)
# is a comparison-layer concern and does not live here.
#
# Public API (exported from GIModel.jl):
#   annihilation_basis_input, isoscalar_annihilation_block,
#   pseudoscalar_annihilation_block, add_isoscalar_annihilation,
#   compute_isoscalar_spectrum

# The sqrt(2) flavor-coherence factor for the isoscalar
# n nbar = (u ubar + d dbar)/sqrt(2) combination is explicit input metadata;
# labels only reproduce the "1 ns" / "1 ss" convention used by reports.
function _annihilation_flavor_tag(m::Meson)
    is_equal_flavor(m) || throw(ArgumentError(
        "annihilation basis inputs need self-conjugate flavor content, got $(flavor_label(m))",
    ))
    m.flavor1 in (:u, :d, :q) && return "ns"
    return string(m.flavor1, m.flavor1)
end

"""
    annihilation_basis_input(spec::Spectrum, level::BasisState)

One flavor-channel entry for an annihilation mixing block: the state's
pre-annihilation model mass as the diagonal, and the signed projection of its
physical radial components into the requested `(L,S,J)` annihilation channel.
Each wave is phase-fixed to the GI convention `Φ(0) > 0` (see
[`fix_annihilation_phase!`](@ref)), with the compensating sign retained in its
coefficient. The
eigensolver hands back arbitrary column signs, and without the convention the
`S_L` factor flips between quark masses and radial levels, randomizing the
off-diagonal block elements.

The wave comes from the spectrum's own solver, like every other observable.

**Historical note, because the alternative looks principled and is not.** This
used to take a `wave_basis` keyword defaulting to `:ho`, so annihilation read a
separately-cached oscillator wave no matter which solver produced the spectrum.
The stated reason was that the oscillator basis gave a ~3x larger
wavefunction-at-origin for light quarks and Table III was built on that basis.
The real cause was a normalization bug: finite differences then returned
Euclidean eigenvectors (`sum u^2 = 1`) against the oscillator path's physical
ones (`int u^2 dr = 1`), a ratio of `1/sqrt(h)` = 3.17 on the Table III audit's
own 220-point mesh. That bug is fixed (see `physically_normalized_waves`), and
the two solvers now agree on the smeared origin factor to 0.1%, so there is
nothing left for a second basis to correct.
"""
# `fix_annihilation_phase!` works on eigenvector columns; a cached wave must not
# be mutated in place, so flip a copy.
function fix_annihilation_phase(w::MeshWave)
    u = reshape(copy(w.u), :, 1)
    fix_annihilation_phase!(u, w.r)
    return MeshWave(vec(u), w.r, w.h)
end

function fix_annihilation_phase(w::OscillatorWave)
    phase, _ = quadgk(r -> r * _oscillator_radial_value(w, r), 0.0, Inf; rtol = 1e-10)
    phase >= 0 && return w
    return OscillatorWave(w.L, w.beta, -w.coefficients)
end

function _annihilation_phase_sign(w::MeshWave)
    return sum(w.r .* w.u) < 0 ? -1.0 : 1.0
end

function _annihilation_phase_sign(w::OscillatorWave)
    phase, _ = quadgk(r -> r * _oscillator_radial_value(w, r), 0.0, Inf; rtol = 1e-10)
    return phase < 0 ? -1.0 : 1.0
end

function annihilation_basis_input(spec::SpinResolvedSpectrum, level::BasisState)
    meson = _single_channel(spec)
    state = spectrum_state(spec, level)
    components = [
        component for component in physical_components(spec, state) if
        component.basis.L_label == level.L_label &&
        component.basis.multiplicity == level.multiplicity &&
        component.basis.J == level.J
    ]
    isempty(components) && throw(ArgumentError(
        "$(state.label) has no $(level.L_label), multiplicity=$(level.multiplicity), J=$(level.J) component for annihilation",
    ))
    # Express every component in the annihilation phase convention without
    # changing relative physical phases, then fix the remaining global state
    # phase by making the first component coefficient positive.
    radial_components = [
        begin
            sign = _annihilation_phase_sign(component.wave)
            (component.coefficient * sign, fix_annihilation_phase(component.wave))
        end for component in components
    ]
    global_sign = first(radial_components)[1] < 0 ? -1.0 : 1.0
    radial_components = [
        (global_sign * coefficient, wave) for (coefficient, wave) in radial_components
    ]
    label = "$(level.n) $(_annihilation_flavor_tag(meson))"
    return pseudoscalar_annihilation_basis_input(
        BasisState(
            level.n, level.L_label, level.multiplicity, level.J;
            label = label,
            flavors = (meson.flavor1, meson.flavor2),
        ),
        meson.constituent_masses.m1_GeV,
        state.mass_GeV,
        radial_components;
        # Coherent exactly for the nonstrange (u ubar + d dbar)/sqrt(2) channel —
        # the same predicate `_annihilation_flavor_tag` uses to emit "ns", so this
        # reproduces the old label-substring behavior identically.
        isoscalar_coherent = meson.flavor1 in (:u, :d, :q),
    )
end

"""
    isoscalar_annihilation_block(params, nn, ss, level; amplitude_A)

General Eq. (16) two-flavor annihilation block for one ``n\\,^{2S+1}L_J``
channel across the `nn` and `ss` [`Spectrum`](@ref)s (basis order `[nn, ss]`).
Wraps [`isoscalar_general_annihilation_solution`](@ref); the caller supplies the
channel amplitude `A(^{2S+1}L_J)`.
"""
function isoscalar_annihilation_block(
    params::GIParameters,
    nn::SpinResolvedSpectrum,
    ss::SpinResolvedSpectrum,
    level::BasisState;
    amplitude_A::Real,
)
    basis = [
        annihilation_basis_input(nn, level),
        annihilation_basis_input(ss, level),
    ]
    return isoscalar_general_annihilation_solution(
        params,
        basis;
        amplitude_A = amplitude_A,
        L = L_SYMBOLS[level.L_label],
        multiplicity = level.multiplicity,
        J = level.J,
    )
end

function _pseudoscalar_block_basis(nn::SpinResolvedSpectrum, ss::SpinResolvedSpectrum)
    return [
        annihilation_basis_input(nn, BasisState(1, "S", 1, 0)),
        annihilation_basis_input(ss, BasisState(1, "S", 1, 0)),
        annihilation_basis_input(nn, BasisState(2, "S", 1, 0)),
        annihilation_basis_input(ss, BasisState(2, "S", 1, 0)),
    ]
end

"""
    pseudoscalar_annihilation_block(model, params, nn, ss; targets=nothing)

`^1S_0` isoscalar annihilation block over the `[1 nn̄, 1 ss̄, 2 nn̄, 2 ss̄]`
basis built from two [`Spectrum`](@ref)s. Dispatches on the
[`PseudoscalarAnnihilationModel`](@ref): [`CalibratedP1Annihilation`](@ref)
requires explicit `targets` (four masses to calibrate the rank-one block to —
digitized paper values live in the comparison layer, not here);
[`PaperP1Annihilation`](@ref)/[`PaperP2Annihilation`](@ref) use the Eq. (18a,b)
formulas with the cached waves.
"""
function pseudoscalar_annihilation_block(
    ::CalibratedP1Annihilation,
    params::GIParameters,
    nn::SpinResolvedSpectrum,
    ss::SpinResolvedSpectrum;
    targets = nothing,
)
    isnothing(targets) && throw(ArgumentError(
        "CalibratedP1Annihilation requires explicit `targets` (four pseudoscalar masses in GeV)",
    ))
    basis = _pseudoscalar_block_basis(nn, ss)
    diagonal = [input.diagonal_GeV for input in basis]
    return isoscalar_pseudoscalar_annihilation_solution(
        CalibratedP1Annihilation(),
        diagonal;
        targets = targets,
        basis = [input.basis for input in basis],
    )
end

function pseudoscalar_annihilation_block(
    model::PaperP1Annihilation,
    params::GIParameters,
    nn::SpinResolvedSpectrum,
    ss::SpinResolvedSpectrum;
    targets = nothing,
)
    basis = _pseudoscalar_block_basis(nn, ss)
    return isoscalar_pseudoscalar_annihilation_solution(model, params, basis)
end

function pseudoscalar_annihilation_block(
    model::PaperP2Annihilation,
    params::GIParameters,
    nn::SpinResolvedSpectrum,
    ss::SpinResolvedSpectrum;
    targets = nothing,
)
    basis = _pseudoscalar_block_basis(nn, ss)
    return isoscalar_pseudoscalar_annihilation_solution(model, params, basis)
end

function _validate_isoscalar_channels(nn::MixedSpectrum, ss::MixedSpectrum)
    nn_meson = _single_channel(nn)
    ss_meson = _single_channel(ss)
    is_equal_flavor(nn_meson) && nn_meson.flavor1 in (:u, :d, :q) ||
        throw(ArgumentError(
            "first spectrum must be the self-conjugate nonstrange isoscalar channel",
        ))
    is_equal_flavor(ss_meson) && ss_meson.flavor1 == :s || throw(ArgumentError(
        "second spectrum must be the self-conjugate strange channel",
    ))
    parameters(nn) == parameters(ss) || throw(ArgumentError(
        "isoscalar flavor channels must use identical GI parameters",
    ))
    nn.computation.solver == ss.computation.solver || throw(ArgumentError(
        "isoscalar flavor channels must use the same radial solver settings",
    ))
    return nn_meson, ss_meson
end

function _combined_isoscalar_spectrum(nn::MixedSpectrum, ss::MixedSpectrum)
    nn_meson, ss_meson = _validate_isoscalar_channels(nn, ss)
    cache = copy(nn.computation.channel_cache)
    for (key, solution) in ss.computation.channel_cache
        haskey(cache, key) && cache[key] !== solution && throw(ArgumentError(
            "isoscalar channel caches contain two solutions for the same radial key",
        ))
        cache[key] = solution
    end
    computation = SectorComputation(parameters(nn), nn.computation.solver, cache)
    return Spectrum([nn_meson, ss_meson], vcat(nn.states, ss.states), computation)
end

function _apply_annihilation_result!(spec::MixedSpectrum, result::MixingResult)
    occursin("annihilation", result.block.mechanism) || throw(ArgumentError(
        "$(result.block.name) is not an annihilation mixing result",
    ))
    all(!isnothing(basis.flavors) for basis in result.block.basis) ||
        throw(ArgumentError(
            "annihilation mixing basis must carry explicit flavor identity",
        ))
    indices = [
        begin
            state = spectrum_state(spec, basis)
            idx = findfirst(candidate -> candidate === state, spec.states)
            isnothing(idx) && error("internal state lookup failure")
            idx
        end for basis in result.block.basis
    ]
    length(unique(indices)) == length(indices) || throw(ArgumentError(
        "annihilation mixing basis resolves more than once to the same spectrum state",
    ))
    any(
        mixing -> occursin("annihilation", mixing.mechanism),
        Iterators.flatten(spec.states[i].mixings for i in indices),
    ) && throw(ArgumentError(
        "a state cannot participate in more than one annihilation block",
    ))
    _assign_block_members!(spec.states, indices, result)
    return spec
end

function _annihilation_channel_key(key)
    length(key) == 3 || throw(ArgumentError(
        "annihilation amplitude keys must be (L_label, multiplicity, J)",
    ))
    L_label, multiplicity, J = String(key[1]), Int(key[2]), Int(key[3])
    haskey(L_SYMBOLS, L_label) || throw(ArgumentError(
        "unknown annihilation orbital label `$L_label`",
    ))
    return (L_label, multiplicity, J)
end

function _common_radial_levels(
    nn::MixedSpectrum,
    ss::MixedSpectrum,
    key::Tuple{String,Int,Int},
)
    L_label, multiplicity, J = key
    matching(spec) = Set(
        state.n for state in spec.states if
        state.L == L_label && state.multiplicity == multiplicity && state.J == J
    )
    return sort!(collect(intersect(matching(nn), matching(ss))))
end

"""
    add_isoscalar_annihilation(params, nn, ss;
        pseudoscalar=nothing, pseudoscalar_targets=nothing, amplitudes=Dict())

Combine solved nonstrange and strange self-conjugate spectra and apply the
requested Eq. (16)-Eq. (18) flavor-annihilation blocks. The returned
[`MixedSpectrum`](@ref) owns both flavor channels, final masses, and the shared
[`MixingResult`](@ref)s used by [`physical_components`](@ref).

Nothing is selected from reference data: eigenstates are assigned by ascending
model mass to ascending pre-annihilation model mass. `pseudoscalar` must be an
explicit `PaperP1Annihilation`, `PaperP2Annihilation`, or calibrated control;
the latter also requires explicit `pseudoscalar_targets`. `amplitudes` maps
`(L_label, multiplicity, J)` to the general Eq. (16) amplitude and is applied
to every radial level present in both flavor channels.
"""
function add_isoscalar_annihilation(
    params::GIParameters,
    nn::MixedSpectrum,
    ss::MixedSpectrum;
    pseudoscalar::Union{Nothing,PseudoscalarAnnihilationModel} = nothing,
    pseudoscalar_targets = nothing,
    amplitudes = Dict{Tuple{String,Int,Int},Float64}(),
)
    params == parameters(nn) == parameters(ss) || throw(ArgumentError(
        "explicit parameters must match both solved flavor channels",
    ))
    isnothing(pseudoscalar) && !isnothing(pseudoscalar_targets) &&
        throw(ArgumentError("pseudoscalar_targets require an explicit pseudoscalar model"))
    spec = _combined_isoscalar_spectrum(nn, ss)
    if !isnothing(pseudoscalar)
        result = if pseudoscalar isa CalibratedP1Annihilation
            isnothing(pseudoscalar_targets) && throw(ArgumentError(
                "CalibratedP1Annihilation requires explicit pseudoscalar_targets",
            ))
            pseudoscalar_annihilation_block(
                pseudoscalar, params, nn, ss; targets = pseudoscalar_targets,
            )
        else
            isnothing(pseudoscalar_targets) || throw(ArgumentError(
                "pseudoscalar_targets are only valid for CalibratedP1Annihilation",
            ))
            pseudoscalar_annihilation_block(pseudoscalar, params, nn, ss)
        end
        _apply_annihilation_result!(spec, result)
    end
    for raw_key in sort!(collect(keys(amplitudes)); by = string)
        key = _annihilation_channel_key(raw_key)
        radial_levels = _common_radial_levels(nn, ss, key)
        isempty(radial_levels) && throw(ArgumentError(
            "annihilation channel $key has no common nonstrange/strange radial levels",
        ))
        for n in radial_levels
            level = BasisState(n, key...)
            result = isoscalar_annihilation_block(
                params, nn, ss, level; amplitude_A = amplitudes[raw_key],
            )
            _apply_annihilation_result!(spec, result)
        end
    end
    return spec
end

"""
    compute_isoscalar_spectrum(params, nn_meson, ss_meson; ...)

End-to-end model-level isoscalar calculation: solve both native flavor
channels through fixed-sector and spectroscopic mixing, then call
[`add_isoscalar_annihilation`](@ref). No paper reference rows are accepted.
The default is native HO because this is the paper-algorithm entry point;
callers may pass another `RadialSolver` explicitly for an independent comparator.
"""
function compute_isoscalar_spectrum(
    params::GIParameters,
    nn_meson::Meson,
    ss_meson::Meson;
    levels::AbstractVector{BasisState} = spectrum_levels(2),
    solver::RadialSolver = OscillatorSolver(),
    terms::SpinTerms = SpinTerms(),
    pseudoscalar::Union{Nothing,PseudoscalarAnnihilationModel} = nothing,
    pseudoscalar_targets = nothing,
    amplitudes = Dict{Tuple{String,Int,Int},Float64}(),
)
    nn = compute_spectrum(params, nn_meson; levels = levels, solver = solver, terms = terms)
    ss = compute_spectrum(params, ss_meson; levels = levels, solver = solver, terms = terms)
    return add_isoscalar_annihilation(
        params, nn, ss;
        pseudoscalar = pseudoscalar,
        pseudoscalar_targets = pseudoscalar_targets,
        amplitudes = amplitudes,
    )
end
