#!/usr/bin/env julia
# Independent convergence certificate for the modern finite-difference
# comparator. This is deliberately not part of the 1985 HO acceptance path.

using Pkg
Pkg.activate(@__DIR__; io = devnull)

using GIModel
using LinearAlgebra
using Printf

const ROOT = dirname(@__DIR__)
const PARAMS_PATH = default_parameters_path()
const REPORT = joinpath(ROOT, "docs", "residual_reports", "fd_comparator_convergence.md")

const NLEVELS = 6
const ORIGIN_LEVELS = (1, 3, 6)
const RESOLUTION_SWEEP = (
    (600, 32.0),
    (1200, 32.0),
    (1800, 32.0),
    (2400, 32.0),
)
# h = 1/60 GeV^-1 throughout: changing the box does not change resolution.
const EXTENT_SWEEP = (
    (719, 12.0),
    (1199, 20.0),
    (1679, 28.0),
    (1919, 32.0),
    (2159, 36.0),
)
const DEFAULT_GRID = (450, 24.0)
const CERTIFIED_GRID = last(RESOLUTION_SWEEP)
const PENULTIMATE_GRID = RESOLUTION_SWEEP[end - 1]

# Internal FD gates. HO differences are reported but never enter these gates.
const MASS_TOL_MEV = 0.20
const RMS_TOL_PERCENT = 0.10
const ORIGIN_TOL_PERCENT = 0.50
const MIXED_MASS_TOL_MEV = 0.20
const COMPOSITION_DEFECT_TOL = 2.0e-5
const CUTOFF_TOL_PERCENT = 0.01
const TRANSITION_TOL_ABS = 1.0e-4

const CASES = (
    ("q", "1S0", FineStructureMultiplet("S", 1, 0)),
    ("q", "3P2", FineStructureMultiplet("P", 3, 2)),
    ("s", "1S0", FineStructureMultiplet("S", 1, 0)),
    ("s", "3P2", FineStructureMultiplet("P", 3, 2)),
    ("c", "1S0", FineStructureMultiplet("S", 1, 0)),
    ("c", "3P2", FineStructureMultiplet("P", 3, 2)),
    ("b", "1S0", FineStructureMultiplet("S", 1, 0)),
    ("b", "3P2", FineStructureMultiplet("P", 3, 2)),
)

params, mq = load_parameters_and_quark_masses(PARAMS_PATH)

fd_solver(ngrid, rmax) = FiniteDifferenceSolver(
    ngrid = ngrid,
    rmax = rmax,
    kinetic = :relativistic,
    eigensolver = :full,
    nlevels_per_channel = NLEVELS,
)

const SNAPSHOTS = Dict{Tuple{String,String,Int,Float64},NamedTuple}()

function fixed_snapshot(flavor, sector, multiplet, ngrid, rmax)
    key = (flavor, sector, ngrid, rmax)
    return get!(SNAPSHOTS, key) do
        m = mq[flavor]
        solution = fixed_channel_solution(
            params,
            ConstituentMasses(m, m),
            multiplet;
            solver = fd_solver(ngrid, rmax),
            nlevels = NLEVELS,
        )
        rms = [sqrt(radial_expect(radial_wave(solution, n), r -> r^2)) for n = 1:NLEVELS]
        origin = Dict(
            n => wavefunction_origin_smearing(
                radial_wave(solution, n), m;
                L = orbital_angular_momentum(multiplet.L_label),
            ) for n in ORIGIN_LEVELS
        )
        return (
            solution = solution,
            masses = solution.eigenvalues_GeV,
            rms = rms,
            origin = origin,
            max_norm_defect = maximum(
                abs(wave_norm(radial_wave(solution, n)) - 1) for n = 1:NLEVELS
            ),
        )
    end
end

percent_delta(value, reference) = 100 * abs(value - reference) / max(abs(reference), 1.0e-14)

function snapshot_delta(value, reference)
    return (
        mass_MeV = maximum(abs.(value.masses .- reference.masses)) * 1000,
        rms_percent = maximum(percent_delta.(value.rms, reference.rms)),
        origin_percent = maximum(
            percent_delta(value.origin[n], reference.origin[n]) for n in ORIGIN_LEVELS
        ),
    )
end

function calibration_rows()
    rows = NamedTuple[]
    sweep_rows = NamedTuple[]
    ho_rows = NamedTuple[]
    ho_solver = OscillatorSolver(nlevels_per_channel = NLEVELS)
    for (flavor, sector, multiplet) in CASES
        resolution = [
            fixed_snapshot(flavor, sector, multiplet, grid...) for grid in RESOLUTION_SWEEP
        ]
        extent = [
            fixed_snapshot(flavor, sector, multiplet, grid...) for grid in EXTENT_SWEEP
        ]
        default = fixed_snapshot(flavor, sector, multiplet, DEFAULT_GRID...)
        certified = last(resolution)
        spacing_delta = snapshot_delta(resolution[end - 1], resolution[end])
        extent_delta = snapshot_delta(extent[end - 1], extent[end])
        default_delta = snapshot_delta(default, certified)
        passes =
            spacing_delta.mass_MeV <= MASS_TOL_MEV &&
            extent_delta.mass_MeV <= MASS_TOL_MEV &&
            spacing_delta.rms_percent <= RMS_TOL_PERCENT &&
            extent_delta.rms_percent <= RMS_TOL_PERCENT &&
            spacing_delta.origin_percent <= ORIGIN_TOL_PERCENT &&
            extent_delta.origin_percent <= ORIGIN_TOL_PERCENT
        push!(rows, (;
            flavor,
            sector,
            spacing_delta,
            extent_delta,
            default_delta,
            norm_defect = certified.max_norm_defect,
            passes,
        ))

        for (kind, grids, snapshots, endpoint) in (
            ("resolution", RESOLUTION_SWEEP, resolution, last(resolution)),
            ("extent", EXTENT_SWEEP, extent, last(extent)),
        )
            for (grid, snapshot) in zip(grids, snapshots)
                push!(sweep_rows, (;
                    flavor,
                    sector,
                    kind,
                    ngrid = grid[1],
                    rmax = grid[2],
                    h = grid[2] / (grid[1] + 1),
                    delta = snapshot_delta(snapshot, endpoint),
                ))
            end
        end

        ho = fixed_channel_solution(
            params,
            ConstituentMasses(mq[flavor], mq[flavor]),
            multiplet;
            solver = ho_solver,
            nlevels = NLEVELS,
        )
        ho_rms = [sqrt(radial_expect(radial_wave(ho, n), r -> r^2)) for n = 1:NLEVELS]
        L = orbital_angular_momentum(multiplet.L_label)
        ho_origin = Dict(
            n => wavefunction_origin_smearing(radial_wave(ho, n), mq[flavor]; L = L) for
            n in ORIGIN_LEVELS
        )
        push!(ho_rows, (;
            flavor,
            sector,
            mass_MeV = maximum(abs.(certified.masses .- ho.eigenvalues_GeV)) * 1000,
            rms_percent = maximum(percent_delta.(certified.rms, ho_rms)),
            origin_percent = maximum(
                percent_delta(certified.origin[n], ho_origin[n]) for n in ORIGIN_LEVELS
            ),
        ))
    end
    return rows, sweep_rows, ho_rows
end

function origin_at_pmax(wave::MeshWave, mass, L, pmax)
    # Keep the production transform's Δp while varying only its upper limit.
    production_wave = observable_momentum_wave(wave, L; npoints = 900)
    dp = last(production_wave.p) / (length(production_wave.p) - 1)
    npoints = ceil(Int, pmax / dp) + 1
    momentum = momentum_wave(wave, L; pmax = pmax, npoints = npoints)
    kernel = p -> begin
        energy = sqrt(mass^2 + p^2)
        (p / energy)^L * mass / energy
    end
    return sqrt(2 / π) / sqrt(4π) * momentum_functional(momentum, kernel)
end

function cutoff_rows()
    rows = NamedTuple[]
    for (flavor, sector, multiplet) in CASES
        flavor in ("q", "b") || continue
        snapshot = fixed_snapshot(flavor, sector, multiplet, CERTIFIED_GRID...)
        mass = mq[flavor]
        L = orbital_angular_momentum(multiplet.L_label)
        change_45_60 = 0.0
        change_60_80 = 0.0
        for n in ORIGIN_LEVELS
            wave = radial_wave(snapshot.solution, n)
            at45 = origin_at_pmax(wave, mass, L, 45.0)
            at60 = snapshot.origin[n]
            at80 = origin_at_pmax(wave, mass, L, 80.0)
            change_45_60 = max(change_45_60, percent_delta(at45, at60))
            change_60_80 = max(change_60_80, percent_delta(at60, at80))
        end
        push!(rows, (;
            flavor,
            sector,
            change_45_60,
            change_60_80,
            passes = change_60_80 <= CUTOFF_TOL_PERCENT,
        ))
    end
    return rows
end

function coefficient_vector(spec, state, identities)
    values = Dict(
        (c.basis.n, c.basis.L_label, c.basis.multiplicity, c.basis.J, c.basis.flavors) =>
            c.coefficient for c in physical_components(spec, state)
    )
    vector = [get(values, identity, 0.0) for identity in identities]
    return vector / norm(vector)
end

function composition_defect(left_spec, left_state, right_spec, right_state)
    identities = sort!(collect(union(
        Set((c.basis.n, c.basis.L_label, c.basis.multiplicity, c.basis.J, c.basis.flavors)
            for c in physical_components(left_spec, left_state)),
        Set((c.basis.n, c.basis.L_label, c.basis.multiplicity, c.basis.J, c.basis.flavors)
            for c in physical_components(right_spec, right_state)),
    )); by = string)
    left = coefficient_vector(left_spec, left_state, identities)
    right = coefficient_vector(right_spec, right_state, identities)
    return max(0.0, 1 - abs(dot(left, right)))
end

function mixing_cases()
    return (
        (
            label = "q-s same-J P",
            mechanism = "antisymmetric_spin_orbit",
            meson = Meson(mq, :q, :s),
            levels = [
                BasisState(n, "P", mult, 1) for mult in (1, 3) for n in 1:2
            ],
        ),
        (
            label = "c-c tensor P-F",
            mechanism = "tensor_mixing",
            meson = Meson(mq, :c, :c),
            levels = [
                BasisState(n, L, 3, 2) for L in ("P", "F") for n in 1:2
            ],
        ),
    )
end

function mixing_rows()
    rows = NamedTuple[]
    for case in mixing_cases()
        penultimate = compute_spectrum(
            params,
            case.meson;
            levels = case.levels,
            solver = fd_solver(PENULTIMATE_GRID...),
        )
        certified = compute_spectrum(
            params,
            case.meson;
            levels = case.levels,
            solver = fd_solver(CERTIFIED_GRID...),
        )
        ho = compute_spectrum(
            params,
            case.meson;
            levels = case.levels,
            solver = OscillatorSolver(nlevels_per_channel = NLEVELS),
        )
        for state in certified.states
            prior = spectrum_state(penultimate, state.label)
            hostate = spectrum_state(ho, state.label)
            fd_defect = composition_defect(certified, state, penultimate, prior)
            ho_defect = composition_defect(certified, state, ho, hostate)
            fd_mass = 1000 * abs(state.mass_GeV - prior.mass_GeV)
            components = physical_components(certified, state)
            mixing_weight = 1 - maximum(abs2(c.coefficient) for c in components)
            active = any(m -> m.mechanism == case.mechanism, state.mixings)
            push!(rows, (;
                case = case.label,
                state = state.label,
                mixing_weight,
                fd_mass_MeV = fd_mass,
                fd_composition_defect = fd_defect,
                ho_mass_MeV = 1000 * (state.mass_GeV - hostate.mass_GeV),
                ho_composition_defect = ho_defect,
                passes = active && fd_mass <= MIXED_MASS_TOL_MEV &&
                         fd_defect <= COMPOSITION_DEFECT_TOL,
            ))
        end
    end
    return rows
end

function transition_snapshot(solver)
    m = mq["c"]
    masses = ConstituentMasses(m, m)
    singlet = fixed_channel_solution(
        params, masses, FineStructureMultiplet("S", 1, 0);
        solver = solver, nlevels = 2,
    )
    triplet = fixed_channel_solution(
        params, masses, FineStructureMultiplet("S", 3, 1);
        solver = solver, nlevels = 2,
    )
    p2 = fixed_channel_solution(
        params, masses, FineStructureMultiplet("P", 3, 2);
        solver = solver, nlevels = 2,
    )
    function mw(wave)
        wave isa MeshWave ? momentum_wave(wave, 0; pmax = 60.0, npoints = 900) :
                            momentum_wave(wave, 0)
    end
    m1_allowed = m1_transition_moment(
        mw(radial_wave(singlet, 1)),
        mw(radial_wave(triplet, 1)),
        m,
        m,
        ((4 / 3, m),),
    )
    m1_hindered = m1_transition_moment(
        mw(radial_wave(singlet, 1)),
        mw(radial_wave(triplet, 2)),
        m,
        m,
        ((4 / 3, m),),
    )
    e1_allowed = radial_overlap(
        radial_wave(triplet, 1), radial_wave(p2, 1), r -> r,
    )
    return (; m1_allowed, m1_hindered, e1_allowed)
end

function transition_rows()
    penultimate = transition_snapshot(fd_solver(PENULTIMATE_GRID...))
    certified = transition_snapshot(fd_solver(CERTIFIED_GRID...))
    ho = transition_snapshot(OscillatorSolver(nlevels_per_channel = NLEVELS))
    return [
        (
            observable = "allowed M1 moment",
            fd = certified.m1_allowed,
            fd_change = certified.m1_allowed - penultimate.m1_allowed,
            ho = ho.m1_allowed,
            fd_ho = certified.m1_allowed - ho.m1_allowed,
            passes = abs(certified.m1_allowed - penultimate.m1_allowed) <=
                     TRANSITION_TOL_ABS,
        ),
        (
            observable = "hindered 2S-1S M1 moment",
            fd = certified.m1_hindered,
            fd_change = certified.m1_hindered - penultimate.m1_hindered,
            ho = ho.m1_hindered,
            fd_ho = certified.m1_hindered - ho.m1_hindered,
            passes = abs(certified.m1_hindered - penultimate.m1_hindered) <=
                     TRANSITION_TOL_ABS,
        ),
        (
            observable = "allowed 1S-1P E1 radial overlap",
            fd = certified.e1_allowed,
            fd_change = certified.e1_allowed - penultimate.e1_allowed,
            ho = ho.e1_allowed,
            fd_ho = certified.e1_allowed - ho.e1_allowed,
            passes = abs(certified.e1_allowed - penultimate.e1_allowed) <=
                     TRANSITION_TOL_ABS,
        ),
    ]
end

function write_report(calibration, sweeps, cutoff, ho_rows, mixed, transitions)
    mkpath(dirname(REPORT))
    open(REPORT, "w") do io
        println(io, "# Independent FD Comparator Convergence Audit")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/audit_fd_convergence.jl`.")
        println(io, numerics_provenance(
            fd_solver(CERTIFIED_GRID...),
            OscillatorSolver(nlevels_per_channel = NLEVELS),
        ))
        println(io)
        println(io, "This certifies the modern finite-difference comparator on its own terms.")
        println(io, "It is diagnostic evidence for the package and is **not** a step or")
        println(io, "acceptance condition of the original 1985 finite-HO algorithm.")
        println(io)
        println(io, "## Protocol")
        println(io)
        println(io, "The two FD error sources are varied independently:")
        println(io)
        println(io, "1. **Resolution:** `ngrid` increases at fixed `rmax = 32 GeV^-1`.")
        println(io, "2. **Extent:** `rmax` increases at fixed `h = 1/60 GeV^-1`.")
        println(io)
        println(io, "The calibration covers six radial levels in equal-flavor q, s, c, and b")
        println(io, "`1S0` and `3P2` sectors. Thus the endpoints include diffuse light states,")
        println(io, "compact bottomonium, the resummed contact kernel, and the spin-orbit/tensor")
        println(io, "kernels. Wave tests use RMS radii for all six levels and the signed Eq. (17)")
        println(io, "smeared-origin functional for n = 1, 3, and 6.")
        println(io)
        @printf(io, "Internal gates: mass %.2f MeV; RMS radius %.2f%%; smeared origin %.2f%%.\n",
            MASS_TOL_MEV, RMS_TOL_PERCENT, ORIGIN_TOL_PERCENT)
        println(io, "Each is the largest change over the named levels between the final two")
        println(io, "points of a sweep. HO differences are never included in these gates.")
        println(io)
        println(io, "## Internal FD certificate")
        println(io)
        println(io, "| flavor | sector | spacing: mass MeV | spacing: RMS % | spacing: origin % | extent: mass MeV | extent: RMS % | extent: origin % | norm defect | result |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|:---:|")
        for row in calibration
            @printf(io, "| `%s` | `%s` | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.2e | %s |\n",
                row.flavor, row.sector,
                row.spacing_delta.mass_MeV, row.spacing_delta.rms_percent,
                row.spacing_delta.origin_percent, row.extent_delta.mass_MeV,
                row.extent_delta.rms_percent, row.extent_delta.origin_percent,
                row.norm_defect, row.passes ? "pass" : "**FAIL**")
        end
        println(io)
        println(io, "## Observable momentum cutoff")
        println(io)
        println(io, "The production FD grid historically had a Nyquist momentum just below")
        println(io, "60 GeV. FD-COMP caps finer-grid observable transforms at 60 GeV so that")
        println(io, "coordinate refinement does not enlarge and coarsen a fixed-size momentum")
        println(io, "quadrature. This independent cutoff sweep keeps Δp fixed and checks the")
        println(io, "q/b endpoint flavors at n = 1, 3, and 6. The 60→80 GeV change must remain")
        @printf(io, "below %.2f%%.\n", CUTOFF_TOL_PERCENT)
        println(io)
        println(io, "| flavor | sector | max 45→60 GeV change % | max 60→80 GeV change % | result |")
        println(io, "|---|---|---:|---:|:---:|")
        for row in cutoff
            @printf(io, "| `%s` | `%s` | %.6f | %.6f | %s |\n",
                row.flavor, row.sector, row.change_45_60, row.change_60_80,
                row.passes ? "pass" : "**FAIL**")
        end
        println(io)
        println(io, "## What the historical default costs")
        println(io)
        println(io, "The following is informational: `(450, 24)` is compared with the certified")
        println(io, "`(2400, 32)` endpoint. It explains what improves when FD is used as a")
        println(io, "precision comparator; it does not change the native-HO paper result.")
        println(io)
        println(io, "| flavor | sector | max mass difference MeV | max RMS difference % | max origin difference % |")
        println(io, "|---|---|---:|---:|---:|")
        for row in calibration
            @printf(io, "| `%s` | `%s` | %.3f | %.3f | %.3f |\n",
                row.flavor, row.sector, row.default_delta.mass_MeV,
                row.default_delta.rms_percent, row.default_delta.origin_percent)
        end
        println(io)
        println(io, "## Complete sweeps")
        println(io)
        println(io, "Every row is measured against the endpoint of its own sweep.")
        println(io)
        println(io, "| kind | flavor | sector | ngrid | rmax | h | max mass MeV | max RMS % | max origin % |")
        println(io, "|---|---|---|---:|---:|---:|---:|---:|---:|")
        for row in sweeps
            @printf(io, "| %s | `%s` | `%s` | %d | %.0f | %.6f | %.4f | %.4f | %.4f |\n",
                row.kind, row.flavor, row.sector, row.ngrid, row.rmax, row.h,
                row.delta.mass_MeV, row.delta.rms_percent, row.delta.origin_percent)
        end
        println(io)
        println(io, "## Converged FD versus converged native HO")
        println(io)
        println(io, "These are diagnostic method differences, not pass/fail thresholds.")
        println(io)
        println(io, "| flavor | sector | max mass difference MeV | max RMS difference % | max origin difference % |")
        println(io, "|---|---|---:|---:|---:|")
        for row in ho_rows
            @printf(io, "| `%s` | `%s` | %.3f | %.3f | %.3f |\n",
                row.flavor, row.sector, row.mass_MeV, row.rms_percent, row.origin_percent)
        end
        println(io)
        println(io, "## Full mixing pipeline")
        println(io)
        println(io, "The q-s block exercises unequal-mass antisymmetric spin-orbit mixing; the")
        println(io, "c-c block exercises the physically important P/F cross-orbital tensor")
        println(io, "mixing directly. Composition defects are")
        println(io, "`1 - |<c_a|c_b>|`, so an irrelevant whole-state sign cannot fake a failure.")
        println(io, "The non-leading weight `1 - max |c_i|^2` demonstrates that the named")
        println(io, "mechanism is active rather than passing through an unmixed state.")
        println(io)
        println(io, "| block | assigned state | non-leading weight | final FD step MeV | FD composition defect | FD-HO MeV | FD-HO composition defect | result |")
        println(io, "|---|---|---:|---:|---:|---:|---:|:---:|")
        for row in mixed
            @printf(io, "| %s | `%s` | %.3e | %.4f | %.2e | %+.3f | %.2e | %s |\n",
                row.case, row.state, row.mixing_weight, row.fd_mass_MeV,
                row.fd_composition_defect,
                row.ho_mass_MeV, row.ho_composition_defect,
                row.passes ? "pass" : "**FAIL**")
        end
        println(io)
        println(io, "## Cross-state wave functionals")
        println(io)
        println(io, "The allowed and cancellation-sensitive M1 moments use the two separately")
        println(io, "solved spin-distorted S waves; the E1 entry is their cross-L radial integral.")
        println(io, "Absolute changes are shown because a relative error is misleading for the")
        @printf(io, "near-zero hindered amplitude. Every final FD step must be below %.1e.\n",
            TRANSITION_TOL_ABS)
        println(io)
        println(io, "| observable | certified FD | final FD step | native HO | FD-HO | result |")
        println(io, "|---|---:|---:|---:|---:|:---:|")
        for row in transitions
            @printf(io, "| %s | %+.8f | %+.2e | %+.8f | %+.2e | %s |\n",
                row.observable, row.fd, row.fd_change, row.ho, row.fd_ho,
                row.passes ? "pass" : "**FAIL**")
        end
        println(io)
        println(io, "## Conclusion")
        println(io)
        println(io, "FD is now a separately certified modern comparator over eigenvalues,")
        println(io, "short- and long-distance wave functionals, cross-state overlaps, and the")
        println(io, "post-sector mixing pipeline. The production default remains a fast report")
        println(io, "setting; precision cross-checks should use the certified grid above.")
        println(io, "No public mesh access, fallback, or additional wave/spectrum holder was")
        println(io, "introduced. The comparison calls the same representation-independent")
        println(io, "physics consumers as the native-HO path.")
    end
end

function main()
    calibration, sweeps, ho_rows = calibration_rows()
    cutoff = cutoff_rows()
    mixed = mixing_rows()
    transitions = transition_rows()
    write_report(calibration, sweeps, cutoff, ho_rows, mixed, transitions)
    failures = count(row -> !row.passes, calibration) +
               count(row -> !row.passes, cutoff) +
               count(row -> !row.passes, mixed) +
               count(row -> !row.passes, transitions)
    println("wrote ", REPORT)
    println("FD calibration sectors: ", length(calibration), "; mixed states: ", length(mixed))
    failures == 0 || error("FD-COMP failed $failures internal convergence gates; inspect $REPORT")
end

main()
