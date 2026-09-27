# The two kinds of knob a spectrum calculation has, separated because the
# distinction is physical:
#
#   RadialSolver — how the radial problem is solved: which of the two methods,
#                  and how finely. Changing it must NOT move a mass. When it
#                  does, the setting was under-resolved, which is what
#                  `_warn_if_underresolved` and adaptive HO convergence detect.
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
The defaults are fast and accurate to a fraction of an MeV for most states.
`(ngrid, rmax) = (2400, 32.0)` is a precision setting for cross-checks over
the q/s/c/b sectors. Finite differences are an independent comparator, not
part of the paper's oscillator algorithm.

Fields:

  - `ngrid`, `rmax` — the `ngrid` interior points
    `rᵢ = i rmax/(ngrid+1)`, with Dirichlet boundaries at `0` and `rmax`.
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

# Defaults for the oscillator expansion. `nbasis` is the first variational space
# considered; production solves enlarge it automatically until every requested
# eigenvalue is stable. The β grid is a hard search bracket spanning diffuse
# light states and compact bottomonia; the minimum is refined continuously inside
# it and hitting an endpoint is an error rather than a usable solution.
const HO_DEFAULT_NBASIS = 24
const HO_DEFAULT_MAX_NBASIS = 80
const HO_DEFAULT_BASIS_STEP = 8
const HO_DEFAULT_ENERGY_TOLERANCE_GEV = 1.0e-4 # 0.1 MeV
const HO_DEFAULT_BETA_TOLERANCE_GEV = 2.0e-3
const HO_BETA_GRID = collect(0.25:0.10:2.35)

"""
    OscillatorSolver(; nbasis=24, max_nbasis=80, basis_step=8,
        energy_tolerance_GeV=1e-4, beta_grid=0.25:0.10:2.35,
        beta_tolerance_GeV=2e-3, converge=true, nlevels_per_channel=6)

Solve by expansion in harmonic-oscillator radial functions — Godfrey & Isgur's
own method, Eq. (A17). The Hamiltonian is a finite `nbasis × nbasis` matrix, and
the oscillator scale `β` is a variational parameter scanned over `beta_grid`.

**There is no mesh on this path.** `p²` has closed-form oscillator
matrix elements ([`ho_p2_matrix`](@ref)) and the smeared potential is integrated
by Gauss–Laguerre quadrature ([`ho_operator_matrix`](@ref)), so accuracy is set
by the basis/beta refinement controls below. Plotting or export code may
explicitly sample the returned [`OscillatorWave`](@ref); sampling settings are
not solver state.

Fields:

  - `nbasis`, `max_nbasis`, `basis_step` — initial, maximum, and increment of
    the oscillator expansion. With `converge=true`, every requested eigenvalue
    must change by less than `energy_tolerance_GeV` on two successive
    refinements. `converge=false` performs one explicitly unchecked fixed-size
    solve for convergence studies. `nbasis` must be at least the number of
    eigenlevels requested by a solve.
  - `beta_grid` — the `β` candidates, in GeV. One `β` is chosen per sector, the
    one minimizing the highest requested level, following the paper's convention.
    The best grid cell is refined to `beta_tolerance_GeV`; an endpoint optimum
    fails because the declared bracket does not contain the variational minimum.
    A one-element grid explicitly fixes beta and therefore has no bracket to rail.
  - `energy_tolerance_GeV` — maximum basis-refinement change over all requested
    eigenvalues. The default is 0.1 MeV.
  - `nlevels_per_channel` — radial levels kept per orbital channel.
"""
struct OscillatorSolver <: RadialSolver
    nbasis::Int
    max_nbasis::Int
    basis_step::Int
    energy_tolerance_GeV::Float64
    beta_grid::Vector{Float64}
    beta_tolerance_GeV::Float64
    converge::Bool
    nlevels_per_channel::Int
    function OscillatorSolver(;
        nbasis::Integer = HO_DEFAULT_NBASIS,
        basis_step::Integer = HO_DEFAULT_BASIS_STEP,
        max_nbasis::Integer = max(
            HO_DEFAULT_MAX_NBASIS, nbasis + 2 * basis_step,
        ),
        energy_tolerance_GeV::Real = HO_DEFAULT_ENERGY_TOLERANCE_GEV,
        beta_grid::AbstractVector{<:Real} = HO_BETA_GRID,
        beta_tolerance_GeV::Real = HO_DEFAULT_BETA_TOLERANCE_GEV,
        converge::Bool = true,
        nlevels_per_channel::Integer = 6,
    )
        nbasis >= 1 ||
            throw(ArgumentError("OscillatorSolver: nbasis must be ≥ 1, got $nbasis"))
        basis_step >= 1 || throw(ArgumentError(
            "OscillatorSolver: basis_step must be ≥ 1, got $basis_step",
        ))
        max_nbasis >= nbasis || throw(ArgumentError(
            "OscillatorSolver: max_nbasis=$max_nbasis is below nbasis=$nbasis",
        ))
        energy_tolerance_GeV > 0 || throw(ArgumentError(
            "OscillatorSolver: energy_tolerance_GeV must be positive",
        ))
        isempty(beta_grid) &&
            throw(ArgumentError("OscillatorSolver: beta_grid must not be empty"))
        length(beta_grid) == 2 && throw(ArgumentError(
            "OscillatorSolver: beta_grid needs one fixed value or at least three bracket points",
        ))
        all(>(0), beta_grid) || throw(ArgumentError(
            "OscillatorSolver: every beta must be positive, got $(collect(beta_grid))",
        ))
        issorted(beta_grid) || throw(ArgumentError(
            "OscillatorSolver: beta_grid must be sorted",
        ))
        all(>(0), diff(collect(beta_grid))) || throw(ArgumentError(
            "OscillatorSolver: beta_grid must be strictly increasing",
        ))
        beta_tolerance_GeV > 0 || throw(ArgumentError(
            "OscillatorSolver: beta_tolerance_GeV must be positive",
        ))
        nlevels_per_channel >= 1 || throw(ArgumentError(
            "OscillatorSolver: nlevels_per_channel must be ≥ 1, got $nlevels_per_channel",
        ))
        return new(
            Int(nbasis),
            Int(max_nbasis),
            Int(basis_step),
            float(energy_tolerance_GeV),
            collect(Float64, beta_grid),
            float(beta_tolerance_GeV),
            converge,
            Int(nlevels_per_channel),
        )
    end
end

# Settings are values: two solvers with equal fields are the same method, even
# though `beta_grid` is a Vector (whose default struct `==` would be identity).
Base.:(==)(a::OscillatorSolver, b::OscillatorSolver) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(OscillatorSolver))
Base.hash(s::OscillatorSolver, h::UInt) =
    foldr(hash, (getfield(s, f) for f in fieldnames(OscillatorSolver)); init = hash(:OscillatorSolver, h))

# `RadialSolver(...)` builds the default implementation, so every call site that
# predates the split keeps working and keeps meaning finite differences.
RadialSolver(; kwargs...) = FiniteDifferenceSolver(; kwargs...)

"""
    FiniteDifferenceSolver(base; ngrid=..., rmax=...)
    OscillatorSolver(base; nbasis=..., max_nbasis=..., beta_grid=...)

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
    max_nbasis::Integer = base.max_nbasis,
    basis_step::Integer = base.basis_step,
    energy_tolerance_GeV::Real = base.energy_tolerance_GeV,
    beta_grid::AbstractVector{<:Real} = base.beta_grid,
    beta_tolerance_GeV::Real = base.beta_tolerance_GeV,
    converge::Bool = base.converge,
    nlevels_per_channel::Integer = base.nlevels_per_channel,
) = OscillatorSolver(;
    nbasis = nbasis,
    max_nbasis = max_nbasis,
    basis_step = basis_step,
    energy_tolerance_GeV = energy_tolerance_GeV,
    beta_grid = beta_grid,
    beta_tolerance_GeV = beta_tolerance_GeV,
    converge = converge,
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
        " GeV^-1 (h = ", round(s.rmax / (s.ngrid + 1), digits = 5), "), ", s.kinetic,
        ", ", s.eigensolver, ", ", s.nlevels_per_channel, " levels/channel",
    )
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", s::OscillatorSolver)
    basis_text = s.converge ?
        "nbasis = $(s.nbasis):$(s.basis_step):$(s.max_nbasis) adaptive, ΔE ≤ $(1000 * s.energy_tolerance_GeV) MeV" :
        "nbasis = $(s.nbasis) unchecked"
    print(
        io, "OscillatorSolver: ", basis_text, ", beta in [",
        first(s.beta_grid), ", ", last(s.beta_grid), "] GeV (",
        length(s.beta_grid), " bracket points, refined to ", s.beta_tolerance_GeV,
        " GeV), ", s.nlevels_per_channel,
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
