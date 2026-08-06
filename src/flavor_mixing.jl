# Two-meson flavor mixing bridge: build annihilation blocks from Spectrum
# objects, reusing their cached radial solves. The Eq. (16)-(18) numeric kernels
# live in pseudoscalar_annihilation.jl; the paper's channel-amplitude
# prescription (which channel gets which A, ideal-mixing defaults) is a
# comparison-layer concern and does not live here.
#
# Public API (exported from GIModel.jl):
#   annihilation_basis_input, isoscalar_annihilation_block,
#   pseudoscalar_annihilation_block

# The sqrt(2) flavor-coherence factor for the isoscalar n nbar = (u ubar + d dbar)/sqrt(2)
# combination is keyed off the basis label (see _flavor_coherence_factor); these
# labels reproduce the "1 ns" / "1 ss" convention of the Table III audits.
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
pre-annihilation model mass as the diagonal, and its radial wave, phase-fixed to
the GI convention `Φ(0) > 0` (see [`fix_annihilation_phase!`](@ref)). The
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
function _phase_fixed(w::RadialWaveOnUniformMesh)
    u = reshape(copy(w.u), :, 1)
    fix_annihilation_phase!(u, w.r)
    return RadialWaveOnUniformMesh(vec(u), w.r, w.h)
end

function annihilation_basis_input(spec::SpinResolvedSpectrum, level::BasisState)
    state = spectrum_state(spec, level)
    wave = _phase_fixed(radial_wave(spec, level.L_label, level.n))
    label = "$(level.n) $(_annihilation_flavor_tag(spec.meson))"
    return pseudoscalar_annihilation_basis_input(
        label,
        spec.meson.constituent_masses.m1_GeV,
        state.mass_GeV,
        wave;
        # Coherent exactly for the nonstrange (u ubar + d dbar)/sqrt(2) channel —
        # the same predicate `_annihilation_flavor_tag` uses to emit "ns", so this
        # reproduces the old label-substring behavior identically.
        isoscalar_coherent = spec.meson.flavor1 in (:u, :d, :q),
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
