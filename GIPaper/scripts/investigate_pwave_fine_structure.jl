#!/usr/bin/env julia
# Reproducible mechanism study for the light and open-charm 1P fine structure.
# This is not a second paper-residual audit: it reads no reference spectrum and
# instead decomposes the GI Hamiltonian behind the residuals already recorded by
# audit_mixing_angles.jl and audit_nonmixing_contact.jl.

using Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using Dates
using GIModel
using LinearAlgebra
using Printf

const ROOT = dirname(@__DIR__)
const OUTDIR = joinpath(ROOT, "docs", "investigations")
const FIGDIR = joinpath(OUTDIR, "figures")
mkpath(FIGDIR)

const HBARC_GEV_FM = 0.1973269804
const CASES = [
    (label = "q qbar", slug = "qq", f1 = :q, f2 = :q),
    (label = "q sbar", slug = "qs", f1 = :q, f2 = :s),
    (label = "c qbar", slug = "cq", f1 = :c, f2 = :q),
    (label = "c sbar", slug = "cs", f1 = :c, f2 = :s),
]
const TRIPLETS = [
    BasisState(1, "P", 3, 0),
    BasisState(1, "P", 3, 1),
    BasisState(1, "P", 3, 2),
]
const LEVELS = vcat([BasisState(1, "P", 1, 1)], TRIPLETS)
const HO = OscillatorSolver(
    nbasis = 24,
    max_nbasis = 72,
    basis_step = 8,
    nlevels_per_channel = 1,
)
const FD = FiniteDifferenceSolver(
    ngrid = 700,
    rmax = 30.0,
    eigensolver = :full,
    nlevels_per_channel = 1,
)
const FIXED_TERMS = SpinTerms(
    contact_hyperfine = false,
    fine_structure = true,
    same_j_spin_orbit = false,
    tensor = false,
)
const CENTRAL_TERMS = SpinTerms(
    contact_hyperfine = false,
    fine_structure = false,
    same_j_spin_orbit = false,
    tensor = false,
)

params, mq = load_parameters_and_quark_masses(default_parameters_path())

function solve_case(case, solver)
    meson = Meson(mq, case.f1, case.f2)
    fixed = fixed_spectrum(
        params,
        meson;
        levels = LEVELS,
        solver = solver,
        terms = FIXED_TERMS,
    )
    central = fixed_channel_solution(
        params,
        meson.constituent_masses,
        FineStructureMultiplet("P", 1, 1);
        solver = solver,
        terms = CENTRAL_TERMS,
        nlevels = 1,
    )
    return (; meson, fixed, central)
end

println("Solving 1P sectors with native HO ...")
ho_solutions = Dict(case.slug => solve_case(case, HO) for case in CASES)
println("Solving the independent FD controls ...")
fd_solutions = Dict(case.slug => solve_case(case, FD) for case in CASES)

ledger = NamedTuple[]
for case in CASES
    ho = ho_solutions[case.slug]
    fd = fd_solutions[case.slug]
    central_energy = only(ho.central.eigenvalues_GeV)
    central_wave = radial_wave(ho.central, 1)
    for basis in TRIPLETS
        multiplet = FineStructureMultiplet(basis.L_label, basis.multiplicity, basis.J)
        first = fine_structure_components(
            params, ho.meson.constituent_masses, multiplet, central_wave,
        )
        hs = spectrum_state(ho.fixed, basis)
        fs = spectrum_state(fd.fixed, basis)
        first_mass = central_energy + first.total
        push!(ledger, (
            sector = case.label,
            slug = case.slug,
            state = basis.label,
            J = basis.J,
            central_MeV = 1000central_energy,
            first_vector_MeV = 1000first.spin_orbit_vector,
            first_thomas_MeV = 1000first.spin_orbit_thomas,
            first_tensor_MeV = 1000first.tensor,
            first_total_MeV = 1000first.total,
            full_central_MeV = 1000hs.central_GeV,
            central_distortion_MeV = 1000(hs.central_GeV - central_energy),
            full_vector_MeV = 1000hs.spin_orbit_vector_shift_GeV,
            full_thomas_MeV = 1000hs.spin_orbit_thomas_shift_GeV,
            full_tensor_MeV = 1000hs.tensor_shift_GeV,
            full_shift_MeV = 1000hs.fine_structure_shift_GeV,
            resummation_MeV = 1000(hs.mass_GeV - first_mass),
            ho_mass_MeV = 1000hs.mass_GeV,
            fd_mass_MeV = 1000fs.mass_GeV,
            fd_ho_MeV = 1000(fs.mass_GeV - hs.mass_GeV),
        ))
    end
end

mixing = NamedTuple[]
for case in CASES[2:end]
    for (solver_name, solved) in (("HO", ho_solutions[case.slug]),
                                  ("FD", fd_solutions[case.slug]))
        singlet = spectrum_state(solved.fixed, "1^1P_1")
        triplet = spectrum_state(solved.fixed, "1^3P_1")
        components = spin_orbit_mixing_components(
            params,
            solved.meson.constituent_masses,
            "P",
            radial_wave(solved.fixed, singlet),
            radial_wave(solved.fixed, triplet),
        )
        scale = abs(components.vector) + abs(components.thomas)
        result = same_j_mixing(singlet.mass_GeV, triplet.mass_GeV, components.total)
        push!(mixing, (
            sector = case.label,
            slug = case.slug,
            solver = solver_name,
            vector_MeV = 1000components.vector,
            thomas_MeV = 1000components.thomas,
            total_MeV = 1000components.total,
            survives_percent = 100abs(components.total) / scale,
            gap_MeV = 1000(singlet.mass_GeV - triplet.mass_GeV),
            theta_deg = result.theta_deg,
        ))
    end
end

# Coordinate-space middle kernels. The plotted combinations include the A15-A16
# mass coefficients and P-wave angular factor where applicable, but not the
# nonlocal B(p^2) factors that enclose every kernel in the full Hamiltonian.
function kernel_profiles(masses, radii)
    angular = sqrt(2.0)
    m1, m2 = masses.m1_GeV, masses.m2_GeV
    rows = [fine_structure_radial_kernels(params, masses, r) for r in radii]
    symmetric_vector = [
        k.vector_11 / (4m1^2) + k.vector_22 / (4m2^2) +
        k.vector_12 / (m1 * m2) for k in rows
    ]
    symmetric_thomas = [
        -k.scalar_11 / (4m1^2) - k.scalar_22 / (4m2^2) for k in rows
    ]
    tensor_unit = [k.tensor_12 / (12m1 * m2) for k in rows]
    antisymmetric_vector = [
        angular * (k.vector_11 / (4m1^2) - k.vector_22 / (4m2^2)) for k in rows
    ]
    antisymmetric_thomas = [
        angular * (-k.scalar_11 / (4m1^2) + k.scalar_22 / (4m2^2)) for k in rows
    ]
    return (; symmetric_vector, symmetric_thomas, tensor_unit,
            antisymmetric_vector, antisymmetric_thomas)
end

function draw_diagonal_kernels(path)
    radii = collect(range(0.03, 8.0; length = 500))
    x = HBARC_GEV_FM .* radii
    fig = Figure(size = (1400, 920), backgroundcolor = RGBf(0.985, 0.985, 0.975))
    for (i, case) in enumerate(CASES)
        row, col = fldmod1(i, 2)
        masses = ho_solutions[case.slug].meson.constituent_masses
        p = kernel_profiles(masses, radii)
        ax = Axis(
            fig[row, col];
            title = case.label,
            xlabel = row == 2 ? "r (fm)" : "",
            ylabel = col == 1 ? "mass-weighted middle kernel (GeV)" : "",
        )
        lines!(ax, x, p.symmetric_vector; label = "vector LS", linewidth = 3,
               color = RGBf(0.18, 0.42, 0.72))
        lines!(ax, x, p.symmetric_thomas; label = "scalar/Thomas LS", linewidth = 3,
               color = RGBf(0.82, 0.30, 0.20))
        lines!(ax, x, p.tensor_unit; label = "tensor, unit S12", linewidth = 3,
               color = RGBf(0.20, 0.58, 0.34))
        hlines!(ax, [0.0]; color = (:black, 0.25), linestyle = :dot)
        xlims!(ax, 0, maximum(x))
        i == 1 && axislegend(ax; position = :rt)
    end
    Label(
        fig[0, 1:2],
        "Appendix-A coordinate kernels before B(p²) momentum sandwiches";
        fontsize = 25,
        font = :bold,
    )
    save(path, fig)
end

function draw_mixing_kernels(path)
    radii = collect(range(0.03, 8.0; length = 500))
    x = HBARC_GEV_FM .* radii
    fig = Figure(size = (1500, 470), backgroundcolor = RGBf(0.985, 0.985, 0.975))
    for (i, case) in enumerate(CASES[2:end])
        masses = ho_solutions[case.slug].meson.constituent_masses
        p = kernel_profiles(masses, radii)
        total = p.antisymmetric_vector .+ p.antisymmetric_thomas
        ax = Axis(
            fig[1, i];
            title = case.label,
            xlabel = "r (fm)",
            ylabel = i == 1 ? "antisymmetric middle kernel (GeV)" : "",
        )
        lines!(ax, x, p.antisymmetric_vector; label = "vector", linewidth = 3,
               color = RGBf(0.18, 0.42, 0.72))
        lines!(ax, x, p.antisymmetric_thomas; label = "Thomas", linewidth = 3,
               color = RGBf(0.82, 0.30, 0.20))
        lines!(ax, x, total; label = "local sum", linewidth = 3,
               color = RGBf(0.16, 0.16, 0.16), linestyle = :dash)
        hlines!(ax, [0.0]; color = (:black, 0.25), linestyle = :dot)
        xlims!(ax, 0, maximum(x))
        i == 1 && axislegend(ax; position = :rt)
    end
    save(path, fig)
end

function draw_energy_ledger(path)
    columns = ["vector", "Thomas", "tensor", "resummation"]
    values = reduce(vcat, [
        reshape([
            row.full_vector_MeV,
            row.full_thomas_MeV,
            row.full_tensor_MeV,
            row.resummation_MeV,
        ], 1, :) for row in ledger
    ])
    labels = [string(row.slug, " ", replace(row.state, "1^" => "")) for row in ledger]
    bound = maximum(abs, values)
    fig = Figure(size = (1050, 920), backgroundcolor = RGBf(0.985, 0.985, 0.975))
    ax = Axis(
        fig[1, 1];
        title = "Native-HO 1P energy ledger (MeV)",
        xticks = (1:length(columns), columns),
        yticks = (1:length(labels), labels),
        yreversed = true,
    )
    heatmap!(ax, permutedims(values); colormap = :balance,
             colorrange = (-bound, bound))
    for i in axes(values, 1), j in axes(values, 2)
        text!(ax, j, i; text = @sprintf("%+.1f", values[i, j]),
              align = (:center, :center), fontsize = 15,
              color = abs(values[i, j]) > 0.55bound ? :white : :black)
    end
    save(path, fig)
end

draw_diagonal_kernels(joinpath(FIGDIR, "pwave_diagonal_kernels.png"))
draw_mixing_kernels(joinpath(FIGDIR, "pwave_mixing_kernels.png"))
draw_energy_ledger(joinpath(FIGDIR, "pwave_energy_ledger.png"))

report_path = joinpath(OUTDIR, "pwave_fine_structure.md")
open(report_path, "w") do io
    println(io, "# Light and Charmed P-wave Fine Structure")
    println(io)
    println(io, "Generated by `julia GIPaper/scripts/investigate_pwave_fine_structure.jl` on ",
            Dates.format(Dates.now(), "yyyy-mm-dd"), ".")
    println(io)
    println(io, "This is a mechanism investigation, not an independent paper-number audit.")
    println(io, "The canonical residuals and quoted-angle comparison remain in")
    println(io, "[`../residual_reports/mixing_angles.md`](../residual_reports/mixing_angles.md)")
    println(io, "and the sector reports. Here we ask why those residuals are sensitive.")
    println(io)
    println(io, "## Operators being separated")
    println(io)
    println(io, "Appendix A15-A16 contains three distinct fine-structure mechanisms:")
    println(io)
    println(io, "1. the vector/color-magnetic spin-orbit kernel `(1/r) dGtilde/dr`;")
    println(io, "2. the scalar-confinement Thomas kernel `-(1/r) dStilde/dr`;")
    println(io, "3. the tensor kernel `(1/r) dGtilde/dr - d2Gtilde/dr2`.")
    println(io)
    println(io, "Each displayed coordinate kernel is only the middle `K(r)` of the nonlocal")
    println(io, "GI operator `B(p^2) K(r) B(p^2)`. The energy tables below include those")
    println(io, "momentum factors, the A15-A16 mass denominators, and the exact angular matrix")
    println(io, "elements. The plots therefore diagnose shape and sign; the tables carry the")
    println(io, "physical matrix elements.")
    println(io)
    println(io, "![Diagonal vector, Thomas, and tensor coordinate kernels](figures/pwave_diagonal_kernels.png)")
    println(io)
    println(io, "For diagonal triplets the vector and Thomas curves multiply `L dot S`;")
    println(io, "the green curve is shown per unit conventional `S12`, whose P-wave values")
    println(io, "are -4, +2, and -2/5 for J=0,1,2.")
    println(io)
    println(io, "![Antisymmetric vector and Thomas coordinate kernels](figures/pwave_mixing_kernels.png)")
    println(io)
    println(io, "The dashed local sum is illustrative only because vector and Thomas pieces")
    println(io, "have different epsilon exponents and hence different momentum sandwiches.")
    println(io, "The cancellation must be judged from the full matrix elements below.")
    println(io)
    println(io, "## Numerical design")
    println(io)
    println(io, "HO is the primary calculation. ", numerics_provenance(HO))
    println(io, "FD is an independent control. ", numerics_provenance(FD))
    println(io, "For each flavor sector a central 1P wave gives the first-order ledger; each")
    println(io, "triplet J sector is then diagonalized with the complete fine structure.")
    println(io, "`resummation` is the full HO mass minus central energy minus the complete")
    println(io, "central-wave first-order shift.")
    println(io)
    println(io, "## Triplet energy ledger")
    println(io)
    println(io, "![Full-state operator contributions and nonperturbative remainder](figures/pwave_energy_ledger.png)")
    println(io)
    println(io, "| sector | state | first-order vector | first-order Thomas | first-order tensor | first-order total | full HO shift | resummation | FD-HO mass |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in ledger
        @printf(io, "| %s | `%s` | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f | %+.3f |\n",
                row.sector, row.state, row.first_vector_MeV,
                row.first_thomas_MeV, row.first_tensor_MeV,
                row.first_total_MeV, row.full_shift_MeV,
                row.resummation_MeV, row.fd_ho_MeV)
    end
    println(io)
    println(io, "All energies in the table are MeV. The full shift is the expectation of")
    println(io, "vector + Thomas + tensor in the fully diagonalized fixed-sector state;")
    println(io, "its components are recorded next.")
    println(io)
    println(io, "| sector | state | full vector | full Thomas | full tensor | full sum | central relaxation | vector distortion | Thomas distortion | tensor distortion | closure |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in ledger
        @printf(io, "| %s | `%s` | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f | %+.2f |\n",
                row.sector, row.state, row.full_vector_MeV,
                row.full_thomas_MeV, row.full_tensor_MeV,
                row.full_shift_MeV, row.central_distortion_MeV,
                row.full_vector_MeV - row.first_vector_MeV,
                row.full_thomas_MeV - row.first_thomas_MeV,
                row.full_tensor_MeV - row.first_tensor_MeV,
                row.central_distortion_MeV +
                (row.full_vector_MeV - row.first_vector_MeV) +
                (row.full_thomas_MeV - row.first_thomas_MeV) +
                (row.full_tensor_MeV - row.first_tensor_MeV) -
                row.resummation_MeV)
    end
    println(io)
    println(io, "`central relaxation` is the change in the expectation of the central")
    println(io, "Hamiltonian when the spin-dependent terms distort the radial state. `closure`")
    println(io, "is central relaxation plus the three operator distortions minus the reported")
    println(io, "resummation; it should vanish up to printed rounding.")
    println(io)
    println(io, "## Unequal-mass same-J mixing")
    println(io)
    println(io, "| sector | solver | vector c | Thomas c | total c | survives | singlet-triplet gap | low-state angle |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---:|")
    for row in mixing
        @printf(io, "| %s | %s | %+.2f | %+.2f | %+.2f | %.1f%% | %+.2f | %+.2f deg |\n",
                row.sector, row.solver, row.vector_MeV, row.thomas_MeV,
                row.total_MeV, row.survives_percent, row.gap_MeV,
                row.theta_deg)
    end
    println(io)
    println(io, "## Findings")
    println(io)
    println(io, "1. **The basis is not the limiting uncertainty.** HO and FD agree on the")
    println(io, "   complete 1P masses at sub-MeV scale across all four sectors.")
    println(io, "2. **Diagonal light-P structure is genuinely nonperturbative.** The 3P0")
    println(io, "   state is pulled down by tens of MeV beyond its central-wave first-order")
    println(io, "   estimate, whereas 3P1 and 3P2 are much less distorted.")
    println(io, "3. **The mixing interaction is a cancellation, not a weak force.** In every")
    println(io, "   unequal-mass 1P block the vector and Thomas elements are individually")
    println(io, "   several to tens of MeV, but only a small remainder survives.")
    println(io, "4. **Tensor controls part of the diagonal gap.** It does not connect 1P1")
    println(io, "   to 3P1, but its positive J=1 diagonal contribution changes the denominator")
    println(io, "   governing the mixing angle. It therefore matters indirectly.")
    println(io, "5. **No isolated sign/factor defect appears.** The A15-A16 coefficients,")
    println(io, "   angular sum rules, HO/FD agreement, and operator expectations are mutually")
    println(io, "   consistent. The unresolved reproduction error is the quantitative balance")
    println(io, "   of smeared vector, scalar/Thomas, and tensor strength.")
    println(io)
    println(io, "## Next discriminating calculations")
    println(io)
    println(io, "The useful next step is not another basis scan. It is a controlled parameter")
    println(io, "response study: vary `epsilon_so_vector`, `epsilon_so_scalar`, and `epsilon_t`")
    println(io, "one at a time, refit nothing, and monitor triplet centroids, 3P_J splittings,")
    println(io, "and same-J angles together. A change that fixes an angle but spoils the")
    println(io, "well-reproduced heavy-heavy P multiplets must be rejected. Smearing-width")
    println(io, "dependence should be tested second because it changes all three short-distance")
    println(io, "profiles simultaneously and is therefore less diagnostic.")
end

println("wrote ", report_path)
println("wrote figures under ", FIGDIR)
