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

function _annihilation_wave(
    spec::SpinResolvedSpectrum,
    L_label::AbstractString,
    n::Integer;
    wave_basis::Symbol = :ho,
)
    key = RadialChannelKey(spec.meson.constituent_masses, L_label)
    fd = spec.computation.channel_cache
    haskey(fd, key) || throw(ArgumentError(
        "spectrum for $(flavor_label(spec.meson)) has no `$L_label` channel; include it in `levels`",
    ))
    sol = wave_basis == :ho ?
          get(spec.computation.ho_wave_cache, key, fd[key]) : fd[key]
    n <= size(sol.eigenvectors, 2) || throw(ArgumentError(
        "cached `$L_label` waves for $(flavor_label(spec.meson)) hold $(size(sol.eigenvectors, 2)) levels; requested n=$n",
    ))
    return RadialWaveOnUniformMesh(sol, n)
end

"""
    annihilation_basis_input(spec::Spectrum, level::BasisState; wave_basis=:ho)

One flavor-channel entry for an annihilation mixing block: the state's
pre-annihilation model mass as the diagonal, and the cached radial wave
(harmonic-oscillator basis when available and `wave_basis = :ho`, finite
difference otherwise — the same fallback rule the paper audits use).
"""
function annihilation_basis_input(
    spec::SpinResolvedSpectrum,
    level::BasisState;
    wave_basis::Symbol = :ho,
)
    state = spectrum_state(spec, level)
    wave = _annihilation_wave(spec, level.L_label, level.n; wave_basis = wave_basis)
    label = "$(level.n) $(_annihilation_flavor_tag(spec.meson))"
    return pseudoscalar_annihilation_basis_input(
        label,
        spec.meson.constituent_masses.m1_GeV,
        state.mass_GeV,
        wave,
    )
end

"""
    isoscalar_annihilation_block(params, nn, ss, level; amplitude_A, wave_basis=:ho)

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
    wave_basis::Symbol = :ho,
)
    basis = [
        annihilation_basis_input(nn, level; wave_basis = wave_basis),
        annihilation_basis_input(ss, level; wave_basis = wave_basis),
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

function _pseudoscalar_block_basis(nn::SpinResolvedSpectrum, ss::SpinResolvedSpectrum; wave_basis::Symbol)
    return [
        annihilation_basis_input(nn, BasisState(1, "S", 1, 0); wave_basis = wave_basis),
        annihilation_basis_input(ss, BasisState(1, "S", 1, 0); wave_basis = wave_basis),
        annihilation_basis_input(nn, BasisState(2, "S", 1, 0); wave_basis = wave_basis),
        annihilation_basis_input(ss, BasisState(2, "S", 1, 0); wave_basis = wave_basis),
    ]
end

"""
    pseudoscalar_annihilation_block(model, params, nn, ss; targets=nothing, wave_basis=:ho)

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
    wave_basis::Symbol = :ho,
)
    isnothing(targets) && throw(ArgumentError(
        "CalibratedP1Annihilation requires explicit `targets` (four pseudoscalar masses in GeV)",
    ))
    basis = _pseudoscalar_block_basis(nn, ss; wave_basis = wave_basis)
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
    wave_basis::Symbol = :ho,
)
    basis = _pseudoscalar_block_basis(nn, ss; wave_basis = wave_basis)
    return isoscalar_pseudoscalar_annihilation_solution(model, params, basis)
end

function pseudoscalar_annihilation_block(
    model::PaperP2Annihilation,
    params::GIParameters,
    nn::SpinResolvedSpectrum,
    ss::SpinResolvedSpectrum;
    targets = nothing,
    wave_basis::Symbol = :ho,
)
    basis = _pseudoscalar_block_basis(nn, ss; wave_basis = wave_basis)
    return isoscalar_pseudoscalar_annihilation_solution(model, params, basis)
end
