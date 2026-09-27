include(joinpath(@__DIR__, "common.jl"))

if !isdefined(@__MODULE__, :SPECTRUM_COMPUTED)
    const SPECTRUM_COMPUTED = Ref{Union{Nothing,Dict{String,String}}}(nothing)
    const L_ORDER = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)
    const SPECTRUM_FILES = [
        "reference_spectrum_isovector.csv", "reference_spectrum_strange.csv",
        "reference_spectrum_isoscalar.csv", "reference_spectrum_charmonium.csv",
        "reference_spectrum_charmed.csv", "reference_spectrum_bottomonium.csv",
        "reference_spectrum_b_flavored.csv",
    ]
    const FIGURE_PANELS = Dict(
        "fig_iii" => Set(["isovector"]),
        "fig_iv" => Set(["strange"]),
        "fig_v" => Set(["isoscalar"]),
        "fig_vi" => Set(["charmonium"]),
        "fig_vii" => Set(["charmed", "charmed_strange"]),
        "fig_viii" => Set(["bottomonium"]),
        "fig_ix" => Set(["bottom_light", "bottom_strange", "bottom_charm"]),
    )

    row_key(row) = (String(row.sector), String(row.state), row.n,
        row.multiplicity, String(row.L), row.J)

    function jpc(row)
        l = get(L_ORDER, String(row.L), 0)
        spin = row.multiplicity == 3 ? 1 : 0
        parity = iseven(l + 1) ? '+' : '-'
        if isapprox(row.m1_GeV, row.m2_GeV; rtol=0, atol=1e-12)
            charge = iseven(l + spin) ? '+' : '-'
            return string(row.J, parity, charge)
        end
        return string(row.J, parity)
    end

    function apply_available_mixing(rows, params)
        bykey = Dict(row_key(row) => i for (i, row) in enumerate(rows))
        predicted = [row.predicted_GeV for row in rows]
        touched = falses(length(rows))
        solve_cache = Dict{Tuple{ConstituentMasses,String},Any}()
        groups = Dict{Tuple{String,Int,String,Int,Float64,Float64},Vector{Any}}()
        for row in rows
            l = get(L_ORDER, String(row.L), -1)
            l > 0 && row.J == l || continue
            key = (String(row.sector), row.n, String(row.L), row.J,
                row.m1_GeV, row.m2_GeV)
            push!(get!(groups, key, Any[]), row)
        end
        for group in values(groups)
            singlets = [row for row in group if row.multiplicity == 1]
            triplets = [row for row in group if row.multiplicity == 3]
            for singlet in singlets, triplet in triplets
                masses = ConstituentMasses(singlet.m1_GeV, singlet.m2_GeV)
                solution = get!(solve_cache, (masses, String(singlet.L))) do
                    channel_solution(params, masses, L_ORDER[String(singlet.L)];
                        nlevels=6,
                        solver=FiniteDifferenceSolver(kinetic=:relativistic))
                end
                singlet.n <= length(solution.waves) || continue
                wave = radial_wave(solution, singlet.n)
                offdiag = spin_orbit_mixing_components(
                    params, masses, String(singlet.L), wave).total
                isapprox(offdiag, 0.0; atol=1e-12, rtol=0) && continue
                mixed = same_j_mixing(
                    singlet.predicted_GeV, triplet.predicted_GeV, offdiag)
                si, ti = bykey[row_key(singlet)], bykey[row_key(triplet)]
                if singlet.reference_GeV <= triplet.reference_GeV
                    predicted[si], predicted[ti] = mixed.masses
                else
                    predicted[si], predicted[ti] = reverse(mixed.masses)
                end
                touched[si] = touched[ti] = true
            end
        end
        return [merge(row, (predicted_GeV=predicted[i], mixed=touched[i]))
                for (i, row) in enumerate(rows)]
    end

    function compute_spectra()
        !isnothing(SPECTRUM_COMPUTED[]) && return SPECTRUM_COMPUTED[]
        params, masses = load_parameters_and_quark_masses(default_parameters_path())
        rows = NamedTuple[]
        for filename in SPECTRUM_FILES
            reference = GIPaper.load_reference_spectrum(joinpath(GIPAPER_DIR, "data", filename))
            compared = compare_reference(params, masses, reference;
                contact_hyperfine=true,
                use_fine_structure=params.fine_structure.enabled,
                kinetic=:relativistic)
            append!(rows, apply_available_mixing(compared, params))
        end
        outputs = Dict{String,String}()
        for (stem, panels) in FIGURE_PANELS
            lines = ["sector\tstate\tn\tmultiplicity\tL\tJ\tJPC\tcomputed_GeV\tgi_GeV\tmixed"]
            for row in rows
                String(row.sector) in panels || continue
                push!(lines, join((row.sector, row.state, row.n, row.multiplicity,
                    row.L, row.J, jpc(row), row.predicted_GeV,
                    row.reference_GeV, row.mixed), '\t'))
            end
            outputs[stem] = write_paper_table("$stem.tsv", join(lines, '\n'))
        end
        SPECTRUM_COMPUTED[] = outputs
        return outputs
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && compute_spectra()
