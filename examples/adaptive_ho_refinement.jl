### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook
# outside Pluto, this mock version gives every bound variable its widget default.
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ f15629ac-f5c1-4ca4-94dd-1f30a55c0fd9
begin
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()

    using GIModel
    using CairoMakie
    using PlutoUI
    using Printf
end

# ╔═╡ b5c9aa8f-17c3-48ba-b6f6-fd7bc06836d3
md"""
# How adaptive harmonic-oscillator refinement works

The Godfrey–Isgur calculation expands a radial state in harmonic-oscillator (HO)
functions. A finite expansion is an approximation, so a production answer needs
two numerical choices:

1. the oscillator scale ``\beta``, which sets the size of the basis functions;
2. the number ``N`` of functions retained.

This notebook follows one charmonium ``{}^1S_0`` channel through both choices.
It uses only GIModel's public API. There is **no spatial mesh** in any solve below.
"""

# ╔═╡ 5848082e-d88b-41ae-b263-d3d9feec2e36
md"""
## The whole controller in one picture

```text
choose initial N and a beta bracket
                |
                v
 minimize the highest requested level over beta
                |
                v
 diagonalize the complete fixed-(L,S,J) H_N(beta)
                |
                v
 compare every requested energy with the preceding N
                |
          +-----+------+
          |            |
  tolerance met?      no ----> increase N ----+
          |                                    |
         yes                                   |
          |                                    |
  twice successively? no ----------------------+
          |
         yes
          v
 return ChannelRadialSolution + OscillatorConvergence certificate
```

The common ``\beta`` minimizes the **highest requested eigenvalue**, rather than
optimizing each radial state separately. That preserves a single orthonormal
basis for the channel. Two consecutive acceptable basis enlargements protect
against a single accidental small change.
"""

# ╔═╡ 35c39bd1-1b99-4ed3-a536-477508a7b650
md"""
!!! note "What is held fixed"
    The quark masses, Hamiltonian terms, quantum numbers and requested levels do
    not change during refinement. Only ``\beta`` and ``N`` change. Thus this is
    an accuracy controller around one physics calculation, not a fit and not a
    comparison between different models.
"""

# ╔═╡ a61b7131-b1ae-4f25-a122-3a36328fe3df
md"""## Set up one complete fixed sector"""

# ╔═╡ d5c606a5-1b67-485f-8284-424933469013
params, quark_masses = load_parameters_and_quark_masses(
    joinpath(pkgdir(GIModel), "data", "parameters.provisional.toml"),
);

# ╔═╡ 40199e2b-0b48-442a-82dc-609c81aaff6f
begin
    cc_masses = ConstituentMasses(quark_masses["c"], quark_masses["c"])
    singlet_S = FineStructureMultiplet("S", 1, 0)
    nlevels = 2

    settings = OscillatorSolver(
        nbasis = 24,
        max_nbasis = 80,
        basis_step = 8,
        energy_tolerance_GeV = 1.0e-4, # 0.1 MeV
        beta_grid = collect(0.90:0.05:1.30),
        beta_tolerance_GeV = 2.0e-3,
        nlevels_per_channel = nlevels,
    )
end

# ╔═╡ d1fbb3fb-4df0-4b39-a1d9-dd0d5913ea3a
settings

# ╔═╡ 25950b14-7c04-46dd-b7b1-9f79a9b80d51
md"""
The initial basis has 24 functions. If necessary, the controller tries 32, 40,
and so on through 80. The energy tolerance is 0.1 MeV. The ``\beta`` interval is
a **bracket**, not a sampling grid: after locating the best cell, the solver
continuously refines the minimum to 0.002 GeV.
"""

# ╔═╡ b146650d-aa7a-4b3d-bf70-c45c752ac73d
certified = fixed_channel_solution(
    params,
    cc_masses,
    singlet_S;
    solver = settings,
    nlevels = nlevels,
);

# ╔═╡ 64efc734-ae5e-4834-a199-7d3d79747b7a
certificate = certified.convergence

# ╔═╡ 6a6551c1-c2aa-4b3a-a0bf-7c2616e7fbc0
md"""
## Read the result before looking inside

The production call has returned one ordinary `ChannelRadialSolution`. Its
energies and signed native `OscillatorWave`s are consumed by later spectrum and
observable stages exactly as any other radial solution is. The extra object is a
small certificate:

| certificate field | value | meaning |
|---|---:|---|
| status | `$(certificate.status)` | adaptive acceptance completed |
| final ``N`` | $(certificate.nbasis) | basis size actually used |
| final ``\beta`` | $(@sprintf "%.4f GeV" certificate.beta_GeV) | common optimized scale |
| final energy change | $(@sprintf "%.4g MeV" (1e3 * certificate.energy_delta_GeV)) | worst of the requested levels |
| final overlap defect | $(@sprintf "%.4g" certificate.max_wave_overlap_defect) | worst ``1-|\langle u_N|u_{N-\Delta N}\rangle|`` |
| refinements | $(certificate.refinements) | comparisons made after the initial solve |

The energies decide convergence. Wave overlap is recorded because later
observables can be more sensitive to wave shape, but it is not silently
substituted for the declared energy criterion.
"""

# ╔═╡ a27a9acd-2273-490d-9959-4901f2809886
md"""
!!! warning "What this certificate does—and does not—certify"
    `OscillatorConvergence` certifies the ``\beta`` search and stability under
    enlarging the HO basis. Matrix elements inside each Hamiltonian are evaluated
    by a separate adaptive Gauss–Laguerre quadrature. The current implementation
    can warn that its very strict internal `rtol=1e-10` was not reached before
    its quadrature cap; it then returns the largest-rule value. This notebook
    deliberately leaves that warning visible. It is neither a spatial mesh nor
    an FD fallback, but it is a separate numerical convergence signal and must
    not be claimed as covered by the basis certificate.
"""

# ╔═╡ ae6a80bc-dff5-4b29-9b5c-455509a6dd88
md"""
## Replay the accepted steps

The production controller intentionally returns a compact certificate rather
than retaining a second history object. For teaching, the function below
repeats each visited basis size using `converge=false`. Such a result is marked
`:unchecked`; here it is used only to expose the sequence and independently
recompute the changes between adjacent bases.
"""

# ╔═╡ f67029ec-b6a2-4f0b-b2c0-306d288144e3
function refinement_replay(params, masses, multiplet, settings, final_nbasis; nlevels)
    rows = NamedTuple[]
    previous_energies = nothing
    previous_waves = nothing
    consecutive_passes = 0

    for N in settings.nbasis:settings.basis_step:final_nbasis
        diagnostic_solver = OscillatorSolver(
            settings;
            nbasis = N,
            max_nbasis = N,
            converge = false,
        )
        solution = fixed_channel_solution(
            params,
            masses,
            multiplet;
            solver = diagnostic_solver,
            nlevels = nlevels,
        )
        waves = [radial_wave(solution, n) for n in 1:nlevels]

        if isnothing(previous_energies)
            delta_MeV = nothing
            overlap_defect = nothing
        else
            delta_MeV = 1e3 * maximum(abs.(
                solution.eigenvalues_GeV .- previous_energies,
            ))
            overlap_defect = maximum(
                max(0.0, 1.0 - min(1.0, abs(radial_overlap(
                    previous_waves[n], waves[n], _ -> 1.0,
                )))) for n in 1:nlevels
            )
            passed = delta_MeV <= 1e3 * settings.energy_tolerance_GeV
            consecutive_passes = passed ? consecutive_passes + 1 : 0
        end

        push!(rows, (
            N = N,
            beta_GeV = solution.convergence.beta_GeV,
            energies_GeV = copy(solution.eigenvalues_GeV),
            delta_MeV = delta_MeV,
            overlap_defect = overlap_defect,
            consecutive_passes = consecutive_passes,
        ))
        previous_energies = copy(solution.eigenvalues_GeV)
        previous_waves = waves
    end
    return rows
end

# ╔═╡ 8b8429b5-afaf-45a2-a55f-ad0f32ecef9a
trace = refinement_replay(
    params,
    cc_masses,
    singlet_S,
    settings,
    certificate.nbasis;
    nlevels = nlevels,
);

# ╔═╡ f47a167c-727b-4f69-8de1-b33625bb5a31
begin
    cell(x::Nothing; digits = 4) = "—"
    cell(x::Real; digits = 4) = @sprintf("%.*g", digits, x)

    lines = [
        "| ``N`` | optimized ``\\beta`` (GeV) | ``E_1`` (GeV) | ``E_2`` (GeV) | max ``\\Delta E`` (MeV) | overlap defect | consecutive passes |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in trace
        push!(lines,
            "| $(row.N) | $(cell(row.beta_GeV; digits=6)) | " *
            "$(cell(row.energies_GeV[1]; digits=9)) | " *
            "$(cell(row.energies_GeV[2]; digits=9)) | " *
            "$(cell(row.delta_MeV; digits=4)) | " *
            "$(cell(row.overlap_defect; digits=3)) | " *
            "$(row.consecutive_passes) |",
        )
    end
    Markdown.parse(join(lines, "\n"))
end

# ╔═╡ 83056bfb-e755-4251-9789-ff304c888e4c
begin
    basis_sizes = [row.N for row in trace]
    energy_matrix = reduce(hcat, [row.energies_GeV for row in trace])
    betas = [row.beta_GeV for row in trace]

    fig = Figure(size = (880, 420))
    ax_energy = Axis(
        fig[1, 1];
        xlabel = "HO basis size N",
        ylabel = "fixed-sector energy (GeV)",
        title = "The requested levels settle together",
    )
    for level in 1:nlevels
        lines!(ax_energy, basis_sizes, energy_matrix[level, :]; linewidth = 2)
        scatter!(ax_energy, basis_sizes, energy_matrix[level, :]; markersize = 9,
            label = "level $level")
    end
    axislegend(ax_energy; position = :rt)

    ax_beta = Axis(
        fig[1, 2];
        xlabel = "HO basis size N",
        ylabel = "optimized β (GeV)",
        title = "The scale is re-optimized at every N",
    )
    lines!(ax_beta, basis_sizes, betas; color = :darkorange, linewidth = 2)
    scatter!(ax_beta, basis_sizes, betas; color = :darkorange, markersize = 9)
    fig
end

# ╔═╡ 7e5189e2-43b4-4d35-bbf5-85cb1cb87cb9
md"""
The variational energies approach from above as the basis grows. ``\beta`` may
move at the same time: holding it fixed would confuse basis truncation with a
poor scale choice. The educational replay agrees with the production result to
$(round(1e6 * maximum(abs.(last(trace).energies_GeV .- certified.eigenvalues_GeV)); digits=3)) keV.
"""

# ╔═╡ 1ed95ae8-f980-4ece-beed-288f5b166408
md"""
## Guard rails: failure is part of the algorithm

The controller does not return the last available approximation and call it
converged.

* If the optimum is at either end of `beta_grid`, the bracket has not enclosed
  the minimum and the solve errors with an instruction to widen it.
* If two consecutive energy checks have not passed by `max_nbasis`, the solve
  errors and reports the last change.
* If more eigenlevels are requested than the initial basis can contain, the
  input is rejected.
* `converge=false` is permitted for studies like the replay above, but its
  certificate says `:unchecked`; it is never an adaptive production success.

These are deliberately loud outcomes. A spectrum cannot silently contain an
under-resolved HO state.
"""

# ╔═╡ 290ae600-1434-4f62-82ba-10df9929c55f
md"""
## Check your understanding

**1. What is changed during adaptive refinement?**

$(@bind q1 PlutoUI.Select([
    :unset => "Choose an answer…",
    :a => "The quark masses and spin terms",
    :b => "The HO scale β and basis size N",
    :c => "A spatial mesh spacing and cutoff",
]))

**2. Which energy is minimized when choosing one common ``\beta``?**

$(@bind q2 PlutoUI.Select([
    :unset => "Choose an answer…",
    :a => "Only the ground-state energy",
    :b => "The average of all computed energies",
    :c => "The highest requested eigenvalue",
]))

**3. When is basis convergence accepted?**

$(@bind q3 PlutoUI.Select([
    :unset => "Choose an answer…",
    :a => "After one energy change below tolerance",
    :b => "After two successive refinements below tolerance for every requested level",
    :c => "Whenever the wave overlap defect is small",
]))

**4. What role does the overlap defect play?**

$(@bind q4 PlutoUI.Select([
    :unset => "Choose an answer…",
    :a => "It replaces the energy convergence test",
    :b => "It is a recorded diagnostic for wave-sensitive follow-up",
    :c => "It selects the spatial mesh",
]))

**5. What should happen if the optimal ``\beta`` lies at the bracket edge?**

$(@bind q5 PlutoUI.Select([
    :unset => "Choose an answer…",
    :a => "Use the endpoint as the converged result",
    :b => "Fall back to finite differences",
    :c => "Fail and widen the declared bracket",
]))
"""

# ╔═╡ 598645de-fe79-4221-962d-f53cb1e33874
begin
    selected = (q1, q2, q3, q4, q5)
    correct = (:b, :c, :b, :b, :c)
    answered = count(!=(:unset), selected)
    score = count(identity, map(==, selected, correct))

    answered < length(correct) ?
    md"_Answer all five questions to see your result._" :
    score == length(correct) ?
    md"""!!! correct "5 / 5"
        You can distinguish the variational ``\beta`` search, basis convergence,
        the wave diagnostic and the two loud failure conditions.
    """ :
    md"""!!! warning "$score / $(length(correct))"
        Revisit the controller diagram and the guard-rail section. The key ideas
        are: common ``\beta`` from the highest requested level; two consecutive
        energy passes; overlap recorded but not used as a substitute; no fallback.
    """
end

# ╔═╡ 5c0cbbe2-f395-4d46-a87a-66fd3c22a84a
TableOfContents(; title = "Contents")

# ╔═╡ Cell order:
# ╟─b5c9aa8f-17c3-48ba-b6f6-fd7bc06836d3
# ╠═f15629ac-f5c1-4ca4-94dd-1f30a55c0fd9
# ╟─5848082e-d88b-41ae-b263-d3d9feec2e36
# ╟─35c39bd1-1b99-4ed3-a536-477508a7b650
# ╟─a61b7131-b1ae-4f25-a122-3a36328fe3df
# ╠═d5c606a5-1b67-485f-8284-424933469013
# ╠═40199e2b-0b48-442a-82dc-609c81aaff6f
# ╠═d1fbb3fb-4df0-4b39-a1d9-dd0d5913ea3a
# ╟─25950b14-7c04-46dd-b7b1-9f79a9b80d51
# ╠═b146650d-aa7a-4b3d-bf70-c45c752ac73d
# ╠═64efc734-ae5e-4834-a199-7d3d79747b7a
# ╟─6a6551c1-c2aa-4b3a-a0bf-7c2616e7fbc0
# ╟─a27a9acd-2273-490d-9959-4901f2809886
# ╟─ae6a80bc-dff5-4b29-9b5c-455509a6dd88
# ╠═f67029ec-b6a2-4f0b-b2c0-306d288144e3
# ╠═8b8429b5-afaf-45a2-a55f-ad0f32ecef9a
# ╠═f47a167c-727b-4f69-8de1-b33625bb5a31
# ╠═83056bfb-e755-4251-9789-ff304c888e4c
# ╟─7e5189e2-43b4-4d35-bbf5-85cb1cb87cb9
# ╟─1ed95ae8-f980-4ece-beed-288f5b166408
# ╟─290ae600-1434-4f62-82ba-10df9929c55f
# ╠═598645de-fe79-4221-962d-f53cb1e33874
# ╟─5c0cbbe2-f395-4d46-a87a-66fd3c22a84a
