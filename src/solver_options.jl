# The two kinds of knob a spectrum calculation has, separated because the
# distinction is physical:
#
#   RadialSolver — how the radial problem is solved: which of the two methods,
#                  and how finely. Changing it must NOT move a mass. When it
#                  does, the setting was under-resolved, which is what
#                  `_warn_if_underresolved` and `_warn_if_beta_railed` detect.
#   SpinTerms    — what is in the Hamiltonian. Each switch is a paper equation,
#                  so changing it MUST move a mass.
#
# Keeping both as structs means the defaults live in exactly one place instead
# of being restated at every call site.
#
# Public API (exported from GIModel.jl):
#   RadialSolver, FiniteDifferenceSolver, OscillatorSolver, SpinTerms,
#   numerics_provenance

"""
    RadialSolver

Which method solves the radial Schrödinger equation, and how finely — the
discretization, never the physics content. Two implementations, and **the
choice of implementation lives here, not in [`GIParameters`](@ref)**: the
parameters are the model, the solver is how you solve it.

  - [`FiniteDifferenceSolver`](@ref) — the uniform-mesh solve (the default).
  - [`OscillatorSolver`](@ref) — the paper's own oscillator expansion, Eq. (A17).

Accepted by [`channel_solution`](@ref), [`central_spectrum`](@ref) and
[`compute_spectrum`](@ref), so one object describes a whole calculation. The
complete fixed-sector production path requires the relativistic FD kinetic
term; `kinetic=:nonrelativistic` is intentionally limited to central diagnostics:

    compute_spectrum(params, meson; solver = FiniteDifferenceSolver(ngrid = 900))
    compute_spectrum(params, meson; solver = OscillatorSolver(nbasis = 32))

**Changing a solver setting should never change a physical answer**, and the two
implementations should agree on every quantity. Where they do not, one of them is
under-resolved — that disagreement is a measurement, not a convention.

`RadialSolver(; kwargs...)` builds the default implementation, so
`RadialSolver(ngrid = 900)` still means what it always did.
"""
abstract type RadialSolver end

"""
    FiniteDifferenceSolver(; ngrid=450, rmax=24.0, kinetic=:relativistic,
                             eigensolver=:full, nlevels_per_channel=6)

Solve on a uniform radial mesh. Here the mesh **is** the method: every operator
is a matrix on it, so `ngrid` and `rmax` set the accuracy of the answer.

Fields:

  - `ngrid`, `rmax` — the uniform mesh `r ∈ (0, rmax]` with `ngrid` points.
    Heavy quarkonium is compact and needs points, not reach; light mesons need
    reach. Too coarse and the spin-dependent shifts degrade long before the
    eigenvalues visibly do — see `MIN_POINTS_ACROSS_STATE`.
  - `kinetic` — `:relativistic` (the model's `√(p²+m²)` kinetic term) or
    `:nonrelativistic` (`p²/2μ`, a comparator).
  - `eigensolver` — `:full` (dense `eigen`) or `:krylov`.
  - `nlevels_per_channel` — maximum radial levels kept per requested channel
    (central `L` diagnostic or physical fixed `(L,S,J)` sector).

`kinetic` and `eigensolver` are deliberately absent from
[`OscillatorSolver`](@ref), which has neither choice to make.
"""
struct FiniteDifferenceSolver{Kinetic,Eigensolver} <: RadialSolver
    ngrid::Int
    rmax::Float64
    kinetic::Symbol
    eigensolver::Symbol
    nlevels_per_channel::Int
    function FiniteDifferenceSolver(;
        ngrid::Integer = 450,
        rmax::Real = 24.0,
        kinetic::Symbol = :relativistic,
        eigensolver::Symbol = :full,
        nlevels_per_channel::Integer = 6,
    )
        ngrid >= 2 ||
            throw(ArgumentError("FiniteDifferenceSolver: ngrid must be ≥ 2, got $ngrid"))
        rmax > 0 ||
            throw(ArgumentError("FiniteDifferenceSolver: rmax must be positive, got $rmax"))
        kinetic in (:relativistic, :nonrelativistic) || throw(ArgumentError(
            "FiniteDifferenceSolver: kinetic must be :relativistic or :nonrelativistic, got `$kinetic`",
        ))
        eigensolver in (:full, :krylov) || throw(ArgumentError(
            "FiniteDifferenceSolver: eigensolver must be :full or :krylov, got `$eigensolver`",
        ))
        nlevels_per_channel >= 1 || throw(ArgumentError(
            "FiniteDifferenceSolver: nlevels_per_channel must be ≥ 1, got $nlevels_per_channel",
        ))
        return new{kinetic,eigensolver}(
            Int(ngrid), Float64(rmax), kinetic, eigensolver, Int(nlevels_per_channel),
        )
    end
end

# Defaults for the oscillator expansion. `nbasis = 24` converges every sector the
# paper uses; the β grid deliberately is not tuned per sector, spanning diffuse
# light states and compact bottomonia in one sweep. Both are now solver fields, so
# a user who outgrows them passes a bigger `OscillatorSolver` instead of editing
# this file — which is what `_warn_if_beta_railed` used to have to tell them.
const HO_DEFAULT_NBASIS = 24
const HO_BETA_GRID = collect(0.25:0.10:2.35)

"""
    OscillatorSolver(; nbasis=24, beta_grid=0.25:0.10:2.35, nlevels_per_channel=6)

Solve by expansion in harmonic-oscillator radial functions — Godfrey & Isgur's
own method, Eq. (A17). The Hamiltonian is a finite `nbasis × nbasis` matrix, and
the oscillator scale `β` is a variational parameter scanned over `beta_grid`.

**There is no mesh on this path.** `p²` has closed-form oscillator
matrix elements ([`ho_p2_matrix`](@ref)) and the smeared potential is integrated
by Gauss–Laguerre quadrature ([`ho_operator_matrix`](@ref)), so accuracy is set
by `nbasis` and `beta_grid` alone. Plotting or export code may explicitly sample
the returned [`OscillatorWave`](@ref); sampling settings are not solver state.

Fields:

  - `nbasis` — oscillator functions kept. Raising it can only lower an eigenvalue
    (the calculation is variational), so a mass that keeps falling means the
    basis was too small.
  - `beta_grid` — the `β` candidates, in GeV. One `β` is chosen per sector, the
    one minimizing the highest requested level, following the paper's convention.
    If the optimum lands on an endpoint the basis cannot represent the state and
    `_warn_if_beta_railed` says so.
  - `nlevels_per_channel` — radial levels kept per orbital channel.
"""
struct OscillatorSolver <: RadialSolver
    nbasis::Int
    beta_grid::Vector{Float64}
    nlevels_per_channel::Int
    function OscillatorSolver(;
        nbasis::Integer = HO_DEFAULT_NBASIS,
        beta_grid::AbstractVector{<:Real} = HO_BETA_GRID,
        nlevels_per_channel::Integer = 6,
    )
        nbasis >= 1 ||
            throw(ArgumentError("OscillatorSolver: nbasis must be ≥ 1, got $nbasis"))
        isempty(beta_grid) &&
            throw(ArgumentError("OscillatorSolver: beta_grid must not be empty"))
        all(>(0), beta_grid) || throw(ArgumentError(
            "OscillatorSolver: every beta must be positive, got $(collect(beta_grid))",
        ))
        issorted(beta_grid) || throw(ArgumentError(
            "OscillatorSolver: beta_grid must be sorted; `_warn_if_beta_railed` reads its endpoints",
        ))
        nlevels_per_channel >= 1 || throw(ArgumentError(
            "OscillatorSolver: nlevels_per_channel must be ≥ 1, got $nlevels_per_channel",
        ))
        return new(Int(nbasis), collect(Float64, beta_grid), Int(nlevels_per_channel))
    end
end

# `RadialSolver(...)` builds the default implementation, so every call site that
# predates the split keeps working and keeps meaning finite differences.
RadialSolver(; kwargs...) = FiniteDifferenceSolver(; kwargs...)

"""
    FiniteDifferenceSolver(base; ngrid=..., rmax=...)
    OscillatorSolver(base; nbasis=..., beta_grid=...)

Copy with fields overridden — `FiniteDifferenceSolver(solver; ngrid = 900)` for a
convergence study that changes one knob and keeps the rest.
"""
FiniteDifferenceSolver(
    base::FiniteDifferenceSolver;
    ngrid::Integer = base.ngrid,
    rmax::Real = base.rmax,
    kinetic::Symbol = base.kinetic,
    eigensolver::Symbol = base.eigensolver,
    nlevels_per_channel::Integer = base.nlevels_per_channel,
) = FiniteDifferenceSolver(;
    ngrid = ngrid,
    rmax = rmax,
    kinetic = kinetic,
    eigensolver = eigensolver,
    nlevels_per_channel = nlevels_per_channel,
)

OscillatorSolver(
    base::OscillatorSolver;
    nbasis::Integer = base.nbasis,
    beta_grid::AbstractVector{<:Real} = base.beta_grid,
    nlevels_per_channel::Integer = base.nlevels_per_channel,
) = OscillatorSolver(;
    nbasis = nbasis,
    beta_grid = beta_grid,
    nlevels_per_channel = nlevels_per_channel,
)

"""
    with_mesh(solver, ngrid, rmax) -> FiniteDifferenceSolver

Copy a finite-difference solver onto the mesh already owned by an FD operator.
There is deliberately no oscillator method: an HO solve never consumes a mesh.
"""
with_mesh(s::FiniteDifferenceSolver, ngrid::Integer, rmax::Real) =
    FiniteDifferenceSolver(s; ngrid = ngrid, rmax = rmax)

function Base.show(io::IO, ::MIME"text/plain", s::FiniteDifferenceSolver)
    print(
        io, "FiniteDifferenceSolver: ngrid = ", s.ngrid, ", rmax = ", s.rmax,
        " GeV^-1 (h = ", round(s.rmax / s.ngrid, digits = 5), "), ", s.kinetic,
        ", ", s.eigensolver, ", ", s.nlevels_per_channel, " levels/channel",
    )
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", s::OscillatorSolver)
    print(
        io, "OscillatorSolver: nbasis = ", s.nbasis, ", beta in [",
        first(s.beta_grid), ", ", last(s.beta_grid), "] GeV (",
        length(s.beta_grid), " candidates), ", s.nlevels_per_channel,
        " levels/channel",
    )
    return nothing
end

"""
    SpinTerms(; contact_hyperfine=true, fine_structure=true,
                same_j_spin_orbit=true, tensor=true)

Which spin-dependent terms enter the mass, one switch per paper equation. Unlike
[`RadialSolver`](@ref), **turning one off is meant to change the answer** — that
is how you read off what a term contributes:

    compute_spectrum(params, meson; terms = SpinTerms(tensor = false))

  - `contact_hyperfine` — the smeared spin-spin contact term (stage 2). This is
    the whole `1^1S_0` / `1^3S_1` splitting.
  - `fine_structure` — spin-orbit (vector + Thomas) and tensor shifts (stage 2).
    Gated additionally by `params.fine_structure.enabled`, so a parameter file
    that disables fine structure wins over a `true` here.
  - `same_j_spin_orbit` — the same-`J` antisymmetric spin-orbit mixing (stage 3);
    vanishes identically for equal constituent masses.
  - `tensor` — same-`J` tensor mixing (stage 3).
"""
struct SpinTerms
    contact_hyperfine::Bool
    fine_structure::Bool
    same_j_spin_orbit::Bool
    tensor::Bool
    SpinTerms(;
        contact_hyperfine::Bool = true,
        fine_structure::Bool = true,
        same_j_spin_orbit::Bool = true,
        tensor::Bool = true,
    ) = new(contact_hyperfine, fine_structure, same_j_spin_orbit, tensor)
end

SpinTerms(
    base::SpinTerms;
    contact_hyperfine::Bool = base.contact_hyperfine,
    fine_structure::Bool = base.fine_structure,
    same_j_spin_orbit::Bool = base.same_j_spin_orbit,
    tensor::Bool = base.tensor,
) = SpinTerms(;
    contact_hyperfine = contact_hyperfine,
    fine_structure = fine_structure,
    same_j_spin_orbit = same_j_spin_orbit,
    tensor = tensor,
)

"""
    numerics_provenance(solvers...) -> String

One Markdown line naming the radial solver(s) a report's numbers came from, for
the report header.

A residual report is a measurement, and a measurement without its instrument
settings cannot be checked or reproduced. `mean_abs = 6.0 MeV` says nothing about
whether that was a 450- or 900-point mesh, or the oscillator expansion at
`nbasis = 24` — and those are exactly the knobs a reader would want to rule out
before believing a residual is physics.

Pass every solver that produced numbers in the report; comparison reports pass
both.
"""
function numerics_provenance(solvers::RadialSolver...)
    isempty(solvers) && throw(ArgumentError(
        "numerics_provenance: name the solver(s) that produced the report",
    ))
    return "Numerics: " * join((_solver_phrase(s) for s in solvers), " / ") * "."
end

function _solver_phrase(solver::RadialSolver)
    text = sprint(show, MIME"text/plain"(), solver)
    name, rest = split(text, ": "; limit = 2)
    return "`" * name * "` — " * rest
end

function Base.show(io::IO, ::MIME"text/plain", t::SpinTerms)
    on = [
        name for (name, flag) in (
            ("contact_hyperfine", t.contact_hyperfine),
            ("fine_structure", t.fine_structure),
            ("same_j_spin_orbit", t.same_j_spin_orbit),
            ("tensor", t.tensor),
        ) if flag
    ]
    print(io, "SpinTerms: ", isempty(on) ? "none (central eigenvalues only)" : join(on, ", "))
    return nothing
end
