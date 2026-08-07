# Core layout types: constituent masses, spin multiplet labels, radial FD samples.
# Experimental/catalog row types live in the GIPaper comparison package.

"""
    ConstituentMasses(m1_GeV, m2_GeV)

Constituent quark masses in GeV for a meson **sector** / spin-independent radial channel.

Components promote together and preserve the resulting real values exactly.
Cache normalization belongs to [`RadialChannelKey`](@ref), not to physics input.

This is plain physics input (not a subset of [`GIParameters`](@ref)); masses appear in
spin-dependent operators, smearing widths, and kinetic factors throughout the package.
"""
struct ConstituentMasses{T<:Real}
    m1_GeV::T
    m2_GeV::T
    function ConstituentMasses(m1::Real, m2::Real)
        promoted = promote(m1, m2)
        new{typeof(first(promoted))}(promoted...)
    end
end

"""
    reduced_mass(masses::ConstituentMasses)

Reduced mass ``μ = m_1 m_2 / (m_1 + m_2)`` in GeV (same units as the constituents).
"""
function reduced_mass(m::ConstituentMasses)
    m.m1_GeV * m.m2_GeV / (m.m1_GeV + m.m2_GeV)
end

"""
    FineStructureMultiplet(L_label, multiplicity, J)

Spectroscopic spin/orbital labels shared by spin-dependent corrections:

  - [`fine_structure_components`](@ref) uses `L_label`, `multiplicity`, and `J`.
  - [`contact_hyperfine_shift`](@ref) uses only `L_label` and `multiplicity` (`J` is ignored).

The comparison layer (GIPaper) adds an overload for its reference-catalog rows.
"""
struct FineStructureMultiplet
    L_label::String
    multiplicity::Int
    J::Int
    function FineStructureMultiplet(
        L_label::AbstractString,
        multiplicity::Integer,
        J::Integer,
    )
        new(String(L_label), Int(multiplicity), Int(J))
    end
end

"""
    abstract type RadialWave

One radial eigenlevel, however it was computed. Consumers ask a `RadialWave`
for operations and never touch its representation:

| operation | meaning |
|---|---|
| `radial_expect(w, f)` | integral of u^2 f(r) dr |
| `radial_overlap(wx, wy, f)` | integral of u_x u_y f(r) dr |
| `momentum_wave(w, L)` | the momentum-space wave Phi(p) |
| `momentum_expect(mw, g)` | integral of p^2 |Phi|^2 g(p) dp |
| `wave_norm(w)` | integral of u^2 dr, guaranteed 1 |

Two implementations, in separate files, that never refer to each other:

  - [`MeshWave`](@ref) — samples on a uniform mesh. This is the finite-difference
    solver's native form.
  - [`OscillatorWave`](@ref) — the analytic oscillator representation (β + expansion
    coefficients). It implements an operation **only** when that operation has a
    closed form; anything not yet derived has no method and fails loudly rather
    than quietly discretizing.
"""
abstract type RadialWave end

"""Return the same normalized radial state with its outermost lobe positive."""
function fix_outer_phase end

"""
    MeshWave(u, r, h)
    MeshWave(u, r)

Reduced radial wavefunction ``u(r)`` on a **uniform** interior grid: samples `uᵢ` and radii `rᵢ`
with spacing `h` (for `length(r) ≥ 2`, the two-argument form sets `h = r[2] - r[1]`).

This bundles the data [`fine_structure_components`](@ref), [`contact_hyperfine_shift`](@ref),
and [`physical_u_norm`](@ref) rely on. It is **one radial eigenlevel** on the mesh — not the
full multi-level output of [`channel_solution`](@ref).

For a cached [`ChannelRadialSolution`](@ref), use [`radial_wave`](@ref) to
retrieve the stored native representation. Sample an `OscillatorWave` only for
an explicit plot/export grid with `sample_wave(wave, r)`.

The explicit `h` argument must agree with the uniform spacing implied by `r` (guardrail).
"""
struct MeshWave <: RadialWave
    u::Vector{Float64}
    r::Vector{Float64}
    h::Float64
    function MeshWave(
        u::AbstractVector{<:Real},
        r::AbstractVector{<:Real},
        h::Real,
    )
        length(u) == length(r) ||
            throw(ArgumentError("MeshWave: length(u) != length(r)"))
        length(r) >= 1 || throw(ArgumentError("MeshWave: empty r"))
        hf = float(h)
        isfinite(hf) && hf > 0 ||
            throw(ArgumentError("MeshWave: invalid mesh spacing h=$h"))
        if length(r) >= 2
            hinfer = float(r[2]) - float(r[1])
            isapprox(hinfer, hf; rtol = 1e-10, atol = 1e-12) ||
                throw(
                    ArgumentError(
                        "MeshWave: h=$hf inconsistent with r spacing $hinfer",
                    ),
                )
        end
        return new(collect(Float64, u), collect(Float64, r), hf)
    end
end

function MeshWave(u::AbstractVector{<:Real}, r::AbstractVector{<:Real})
    length(r) >= 2 ||
        throw(ArgumentError("MeshWave(u,r): need length(r) ≥ 2 to infer h"))
    return MeshWave(u, r, r[2] - r[1])
end


"""
    physically_normalized_waves(waves, h) -> Matrix{Float64}

Scale each column of `waves` to the physical radial normalization
`∫u² dr = Σuᵢ² h = 1`.

Every solve in the model returns its waves through this, so "the wave from a
solve" means the same thing in every basis. It exists because the raw
eigensolver hands back Euclidean-normalized columns (`Σuᵢ² = 1`) while the
oscillator path reconstructs already-physical ones — the two differ by `√h`,
and anything quadratic in `u` given the wrong convention is off by `h`.

Zero columns (below-threshold or empty channels) are left alone.
"""
function physically_normalized_waves(waves::AbstractMatrix{<:Real}, h::Real)
    out = Matrix{Float64}(waves)
    hf = float(h)
    hf > 0 || return out
    for col in axes(out, 2)
        nrm = sqrt(sum(abs2, view(out, :, col)) * hf)
        nrm > 0 && (out[:, col] ./= nrm)
    end
    return out
end


# =============================================================================
# The RadialWave interface, implemented for MeshWave by mesh quadrature.
#
# Consumers should go through this interface. Mesh implementation methods may
# read `.u`, `.r`, and `.h`; physics consumers must not. The migration from the
# older mesh-specific API is tracked in `docs/paper_algorithm_work_plan.md`.
# =============================================================================

"""
    wave_norm(w::RadialWave) -> Float64

The physical norm `integral u^2 dr`. Guaranteed to be 1 for any wave produced by
a solve; exposed so the invariant can be asserted rather than assumed.
"""
wave_norm(w::MeshWave) = sum(abs2, w.u) * w.h

function fix_outer_phase(w::MeshWave)
    peak = maximum(abs, w.u)
    index = findlast(x -> abs(x) > 0.2peak, w.u)
    (isnothing(index) || w.u[index] >= 0) && return w
    return MeshWave(-w.u, w.r, w.h)
end

"""
    radial_expect(w::RadialWave, f) -> Float64

`integral u^2 f(r) dr` with `u` normalized. `f` is called as `f(r)`.
"""
function radial_expect(w::MeshWave, f)
    nrm = wave_norm(w)
    nrm > 0 || throw(ArgumentError("radial_expect: zero-norm wave"))
    s = 0.0
    @inbounds for i in eachindex(w.r)
        s += w.u[i]^2 * f(w.r[i])
    end
    return s * w.h / nrm
end

"""
    radial_overlap(wx::RadialWave, wy::RadialWave, f) -> Float64

`integral u_x u_y f(r) dr`, each wave normalized. Both must share a mesh.
"""
function radial_overlap(wx::MeshWave, wy::MeshWave, f)
    length(wx.r) == length(wy.r) ||
        throw(ArgumentError("radial_overlap: waves live on different meshes"))
    isapprox(wx.h, wy.h; rtol = 1e-10) ||
        throw(ArgumentError("radial_overlap: mesh spacings differ"))
    all(isapprox.(wx.r, wy.r; rtol = 1e-10, atol = 1e-12)) ||
        throw(ArgumentError("radial_overlap: mesh points differ"))
    nx, ny = wave_norm(wx), wave_norm(wy)
    (nx > 0 && ny > 0) || throw(ArgumentError("radial_overlap: zero-norm wave"))
    s = 0.0
    @inbounds for i in eachindex(wx.r)
        s += wx.u[i] * wy.u[i] * f(wx.r[i])
    end
    return s * wx.h / sqrt(nx * ny)
end
