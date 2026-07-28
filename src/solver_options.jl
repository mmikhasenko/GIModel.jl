# The two kinds of knob a spectrum calculation has, separated because the
# distinction is physical:
#
#   RadialSolver — how well the radial problem is solved. Changing it must NOT
#                  move a mass. When it does, the grid is under-resolved, which
#                  is what `_warn_if_underresolved` (channel_solver.jl) detects.
#   SpinTerms    — what is in the Hamiltonian. Each switch is a paper equation,
#                  so changing it MUST move a mass.
#
# Keeping both as structs means the defaults live in exactly one place instead
# of being restated at every call site.
#
# Public API (exported from GIModel.jl): RadialSolver, SpinTerms

"""
    RadialSolver(; ngrid=450, rmax=24.0, kinetic=:relativistic, eigensolver=:full,
                   nlevels_per_channel=6)

Numerical settings for the radial solve — the discretization and the method,
never the physics content. Accepted by [`channel_solution`](@ref),
[`central_spectrum`](@ref) and [`compute_spectrum`](@ref), so one object
describes a whole calculation:

    solver = RadialSolver(ngrid = 900, rmax = 32.0)     # a convergence check
    compute_spectrum(params, meson; solver = solver)

Fields:

  - `ngrid`, `rmax` — the uniform radial mesh `r ∈ (0, rmax]` with `ngrid`
    points. Heavy quarkonium is compact and needs points, not reach; light
    mesons need reach. Too coarse and the spin-dependent shifts degrade long
    before the eigenvalues visibly do — see `MIN_POINTS_ACROSS_STATE`.
  - `kinetic` — `:relativistic` (the model's `√(p²+m²)` kinetic term) or
    `:nonrelativistic` (`p²/2μ`, a comparator).
  - `eigensolver` — `:full` (dense `eigen`) or `:krylov`.
  - `nlevels_per_channel` — radial levels kept per orbital channel.

**Changing a `RadialSolver` field should never change a physical answer.** If it
does, the previous setting was under-resolved.
"""
struct RadialSolver
    ngrid::Int
    rmax::Float64
    kinetic::Symbol
    eigensolver::Symbol
    nlevels_per_channel::Int
    function RadialSolver(;
        ngrid::Integer = 450,
        rmax::Real = 24.0,
        kinetic::Symbol = :relativistic,
        eigensolver::Symbol = :full,
        nlevels_per_channel::Integer = 6,
    )
        ngrid >= 2 || throw(ArgumentError("RadialSolver: ngrid must be ≥ 2, got $ngrid"))
        rmax > 0 || throw(ArgumentError("RadialSolver: rmax must be positive, got $rmax"))
        kinetic in (:relativistic, :nonrelativistic) || throw(ArgumentError(
            "RadialSolver: kinetic must be :relativistic or :nonrelativistic, got `$kinetic`",
        ))
        eigensolver in (:full, :krylov) || throw(ArgumentError(
            "RadialSolver: eigensolver must be :full or :krylov, got `$eigensolver`",
        ))
        nlevels_per_channel >= 1 || throw(ArgumentError(
            "RadialSolver: nlevels_per_channel must be ≥ 1, got $nlevels_per_channel",
        ))
        return new(Int(ngrid), Float64(rmax), kinetic, eigensolver, Int(nlevels_per_channel))
    end
end

"""
    RadialSolver(base::RadialSolver; ngrid=..., rmax=...)

Copy with fields overridden — `RadialSolver(solver; ngrid = 900)` for a
convergence study that changes one knob and keeps the rest.
"""
RadialSolver(
    base::RadialSolver;
    ngrid::Integer = base.ngrid,
    rmax::Real = base.rmax,
    kinetic::Symbol = base.kinetic,
    eigensolver::Symbol = base.eigensolver,
    nlevels_per_channel::Integer = base.nlevels_per_channel,
) = RadialSolver(;
    ngrid = ngrid,
    rmax = rmax,
    kinetic = kinetic,
    eigensolver = eigensolver,
    nlevels_per_channel = nlevels_per_channel,
)

function Base.show(io::IO, ::MIME"text/plain", s::RadialSolver)
    print(
        io, "RadialSolver: ngrid = ", s.ngrid, ", rmax = ", s.rmax,
        " GeV^-1 (h = ", round(s.rmax / s.ngrid, digits = 5), "), ", s.kinetic,
        ", ", s.eigensolver, ", ", s.nlevels_per_channel, " levels/channel",
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
