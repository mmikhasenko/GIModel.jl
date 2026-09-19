#!/usr/bin/env julia
# Ten-panel GI meson spectrum.  Colored horizontal bars are newly calculated
# model eigenvalues; open circles and error bars come from fetch_pdg_mesons.py.

using Pkg
Pkg.activate(@__DIR__)

using CairoMakie
using CSV
using DataFrames
using GIModel
using LaTeXStrings

const SIMPLIFIED = "--simplified" in ARGS
const PAPER_ROOT = dirname(@__DIR__)
const PDG_TABLE = joinpath(PAPER_ROOT, "data", "pdg_mesons_2026.csv")
const OUTDIR = joinpath(@__DIR__, "spectrum_plots")
const OUTPUT_STEM = SIMPLIFIED ? "ten_meson_spectra_simplified_pdg2026" :
    "ten_meson_spectra_pdg2026"
const FULL_CALCULATED_TABLE = joinpath(OUTDIR, "ten_meson_calculated_spectrum.csv")
const CALCULATED_TABLE = joinpath(OUTDIR,
    SIMPLIFIED ? "ten_meson_calculated_spectrum_simplified.csv" :
        "ten_meson_calculated_spectrum.csv")

const CLASS_ORDER = [
    "light isoscalar", "light isovector", "charmonium", "bottomonium", "Bc",
    "kaons", "D", "B", "Ds", "Bs",
]

const CLASS_TITLES = Dict(
    "light isoscalar" => "light isoscalar  (n n̄, s s̄)",
    "light isovector" => "light isovector  (n n̄)",
    "kaons" => "kaons  (n s̄)",
    "D" => "D  (c n̄)",
    "Ds" => "Dₛ  (c s̄)",
    "charmonium" => "charmonium  (c c̄)",
    "bottomonium" => "bottomonium  (b b̄)",
    "Bc" => "B꜀  (b c̄)",
    "B" => "B  (b n̄)",
    "Bs" => "Bₛ  (b s̄)",
)

const FLAVORS = Dict(
    "light isovector" => [(:q, :q, "n")],
    "light isoscalar" => [(:q, :q, "n"), (:s, :s, "s")],
    "kaons" => [(:q, :s, "")],
    "D" => [(:c, :q, "")],
    "Ds" => [(:c, :s, "")],
    "B" => [(:b, :q, "")],
    "Bs" => [(:b, :s, "")],
    "Bc" => [(:b, :c, "")],
    "charmonium" => [(:c, :c, "")],
    "bottomonium" => [(:b, :b, "")],
)

const SELF_CONJUGATE = Set([
    "light isoscalar", "light isovector", "charmonium", "bottomonium",
])

# A spectrum envelope like the reference figure: three S radials, two P/D
# radials, and the leading F/G multiplets.
function plotted_levels()
    levels = BasisState[]
    envelope = SIMPLIFIED ?
        (("S", 3), ("P", 2), ("D", 1), ("F", 1)) :
        (("S", 3), ("P", 2), ("D", 2), ("F", 1), ("G", 1))
    for (L, nmax) in envelope
        append!(levels, spectrum_levels(nmax; L_labels = (L,)))
    end
    SIMPLIFIED && filter!(state -> state.J <= 4, levels)
    return levels
end

const FAMILY_COLORS = Dict(
    (1, "S") => "#111111",
    (2, "S") => "#D55E00",
    (3, "S") => "#E69F00",
    (1, "P") => "#0072B2",
    (2, "P") => "#56B4E9",
    (1, "D") => "#009E73",
    (2, "D") => "#44AA99",
    (1, "F") => "#CC79A7",
    (1, "G") => "#7B2CBF",
)

parity(L::AbstractString) = iseven(Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)[String(L)] + 1) ? "+" : "-"
charge_conjugation(L::AbstractString, multiplicity::Integer) =
    iseven(Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)[String(L)] + (multiplicity == 3 ? 1 : 0)) ? "+" : "-"

function quantum_key(J, P, C, class)
    j = Int(J)
    p = SIMPLIFIED && j >= 3 ? "*" : String(P)
    self_conjugate = class in SELF_CONJUGATE
    collapse_c = SIMPLIFIED && self_conjugate && (
        j >= 3 || (j == 1 && String(P) == "+") || (j == 2 && String(P) == "-")
    )
    c = self_conjugate ? (collapse_c ? "*" : String(C)) : ""
    return (j, p, c)
end

function key_label(key)
    J, P, C = key
    P == "*" && C == "*" && return latexstring(J)
    p = P == "*" ? "\\pm" : P
    c = C == "*" ? "\\pm" : C
    return isempty(C) ? latexstring(J, "^{", p, "}") :
        latexstring(J, "^{", p, c, "}")
end

sign_order(sign) = sign == "-" ? 0 : sign == "+" ? 1 : 2
key_sort(key) = (key[1], sign_order(key[2]), sign_order(key[3]))

function model_label(row)
    base = latexstring(row.n, "^{", row.multiplicity, "}", row.L, "_{", row.J, "}")
    isempty(row.flavor_tag) && return base
    return latexstring(row.n, "^{", row.multiplicity, "}", row.L, "_{", row.J,
        ",\\,", row.flavor_tag, "\\bar{", row.flavor_tag, "}}")
end

function calculate_spectra()
    params, mq = load_parameters_and_quark_masses(default_parameters_path())
    levels = plotted_levels()
    solver = FiniteDifferenceSolver(ngrid = 450, rmax = 24.0, kinetic = :relativistic)
    cache = Dict{Tuple{Symbol,Symbol},MixedSpectrum}()
    rows = NamedTuple[]

    for class in CLASS_ORDER
        for (f1, f2, flavor_tag) in FLAVORS[class]
            spectrum = get!(cache, (f1, f2)) do
                println("computing ", f1, " ", f2, "̄ spectrum …")
                compute_spectrum(params, Meson(mq, f1, f2); levels = levels, solver = solver)
            end
            self_conjugate = class in SELF_CONJUGATE
            for state in spectrum.states
                P = parity(state.L)
                C = self_conjugate ? charge_conjugation(state.L, state.multiplicity) : ""
                push!(rows, (
                    class = class,
                    flavor1 = String(f1),
                    flavor2 = String(f2),
                    flavor_tag = flavor_tag,
                    n = state.n,
                    multiplicity = state.multiplicity,
                    L = String(state.L),
                    J = state.J,
                    P = P,
                    C = C,
                    mass_MeV = 1000state.mass_GeV,
                    label = String(state.label),
                ))
            end
        end
    end
    return rows
end

function model_row_in_scope(row)
    SIMPLIFIED || return true
    radial_limit = Dict("S" => 3, "P" => 2, "D" => 1, "F" => 1)
    L = String(row.L)
    return row.J <= 4 && haskey(radial_limit, L) && row.n <= radial_limit[L]
end

function pdg_rows(table::DataFrame, class::String)
    self_conjugate = class in SELF_CONJUGATE
    rows = NamedTuple[]
    for row in eachrow(table)
        String(row.class) == class || continue
        ismissing(row.mass_MeV) && continue
        ismissing(row.J) && continue
        Jtext = String(row.J)
        occursin(r"^\d+$", Jtext) || continue
        Jvalue = parse(Int, Jtext)
        SIMPLIFIED && Jvalue > 4 && continue
        ismissing(row.P) && continue
        P = String(row.P)
        P in ("+", "-") || continue
        C = ismissing(row.C) ? "" : String(row.C)
        self_conjugate && !(C in ("+", "-")) && continue
        ep = ismissing(row.error_plus_MeV) ? 0.0 : Float64(row.error_plus_MeV)
        em = ismissing(row.error_minus_MeV) ? 0.0 : Float64(row.error_minus_MeV)
        push!(rows, (
            key = quantum_key(Jvalue, P, C, class),
            mass_MeV = Float64(row.mass_MeV),
            error_plus_MeV = ep,
            error_minus_MeV = em,
            particle = String(row.particle),
            charge_e = Int(round(Float64(row.charge_e))),
            isospin = ismissing(row.I) ? "" : String(row.I),
            observed = row.observed === true || lowercase(string(row.observed)) == "true",
            multiplet = String(row.multiplet_pdgids),
        ))
    end
    return rows
end

function isospin_offset(row)
    if row.isospin == "1"
        return 0.15 * row.charge_e
    elseif row.isospin == "1/2"
        return row.charge_e == 0 ? -0.11 : 0.11
    end
    return 0.0
end

"""Convert PDG ASCII particle names to compact mathematical labels."""
function pdg_latex(raw::AbstractString)
    s = replace(String(raw), "()" => "")
    charge = ""
    m = match(r"([+\-0])$", s)
    if !isnothing(m)
        charge = m.captures[1]
        s = s[1:(end - 1)]
    end

    for symbol in ("D", "K", "B")
        prefix = string(symbol, "bar")
        if startswith(s, prefix)
            s = string("\\bar{", symbol, "}", s[(length(prefix) + 1):end])
            break
        end
    end
    greek = (
        "Upsilon" => "\\Upsilon", "omega" => "\\omega", "rho" => "\\rho",
        "eta" => "\\eta", "psi" => "\\psi", "chi" => "\\chi",
        "phi" => "\\phi", "pi" => "\\pi",
    )
    if startswith(s, "J/psi")
        s = string("J/\\psi", s[6:end])
    else
        for (ascii, tex) in greek
            if startswith(s, ascii)
                s = string(tex, s[(length(ascii) + 1):end])
                break
            end
        end
    end
    s = replace(s, "^'" => "^{\\prime}", "^*" => "^{*}")
    s = replace(s, r"_[A-Za-z0-9]+" => token -> string("_{", token[2:end], "}"))
    if !isempty(charge)
        modifier = ""
        if occursin("^{*}", s)
            s = replace(s, "^{*}" => "")
            modifier = "*"
        elseif occursin("^{\\prime}", s)
            s = replace(s, "^{\\prime}" => "")
            modifier = "\\prime"
        end
        s = string(s, "^{", modifier, charge, "}")
    end
    return latexstring(s)
end

function draw_panel!(fig, slot, class, model, pdg)
    model_keys = [quantum_key(r.J, r.P, r.C, class) for r in model]
    if SIMPLIFIED
        model_ceiling = Dict{Tuple{Int,String,String},Float64}()
        for (row, key) in zip(model, model_keys)
            model_ceiling[key] = max(get(model_ceiling, key, -Inf), row.mass_MeV)
        end
        pdg = [row for row in pdg if
            haskey(model_ceiling, row.key) && row.mass_MeV <= model_ceiling[row.key]]
    end
    keys = sort!(unique(vcat(model_keys, [r.key for r in pdg])); by = key_sort)
    xof = Dict(key => i for (i, key) in enumerate(keys))

    all_masses = vcat([r.mass_MeV for r in model], [r.mass_MeV for r in pdg])
    ylo, yhi = extrema(all_masses)
    pad = max(65.0, 0.055(yhi - ylo))

    ax = Axis(
        fig[slot...];
        title = CLASS_TITLES[class],
        titlealign = :left,
        titlesize = SIMPLIFIED ? 19 : 16,
        xticks = (1:length(keys), key_label.(keys)),
        xticklabelrotation = 0.0,
        xticklabelsize = SIMPLIFIED ? 15 : 11,
        yticklabelsize = SIMPLIFIED ? 13 : 11,
        yminorticksvisible = true,
        xgridvisible = false,
        ygridvisible = false,
        xticksmirrored = true,
        yticksmirrored = true,
        xtickalign = 1,
        ytickalign = 1,
        topspinevisible = true,
        rightspinevisible = true,
    )
    xlims!(ax, 0.45, length(keys) + 0.55)
    ylims!(ax, ylo - pad, yhi + pad)

    for x in 1:length(keys)
        vlines!(ax, x; color = (:gray70, 0.18), linewidth = 0.6)
    end

    for row in model
        key = quantum_key(row.J, row.P, row.C, class)
        x = xof[key]
        linestyle = class == "light isoscalar" && row.flavor_tag == "s" ? :dash : :solid
        lines!(ax, [x - 0.32, x + 0.32], fill(row.mass_MeV, 2);
            color = FAMILY_COLORS[(row.n, row.L)],
            linewidth = SIMPLIFIED ? 2.7 : 2.0, linestyle)
        text!(ax, x - 0.31, row.mass_MeV + 0.010(yhi - ylo);
            text = model_label(row), align = (:left, :bottom),
            fontsize = SIMPLIFIED ? 9.2 : 7.0,
            color = (:black, 0.82))
    end

    if !isempty(pdg)
        xs = Float64[xof[row.key] + isospin_offset(row) for row in pdg]
        ys = Float64[row.mass_MeV for row in pdg]
        errorbars!(ax, xs, ys,
            Float64[row.error_minus_MeV for row in pdg],
            Float64[row.error_plus_MeV for row in pdg];
            color = (:gray35, 0.34), whiskerwidth = 2.5, linewidth = 0.45)
        observed = findall(row -> row.observed, pdg)
        unobserved = findall(row -> !row.observed, pdg)
        isempty(observed) || scatter!(ax, xs[observed], ys[observed];
            marker = :circle, markersize = SIMPLIFIED ? 6.2 : 5.2, color = :black,
            strokecolor = :black, strokewidth = 0.8)
        isempty(unobserved) || scatter!(ax, xs[unobserved], ys[unobserved];
            marker = :circle, markersize = SIMPLIFIED ? 6.2 : 5.2, color = :white,
            strokecolor = :black, strokewidth = 0.9)

        # Alternating offsets keep dense same-J columns readable while retaining
        # labels for every plotted PDG family, as in the reference figure.
        ranks = Dict{Tuple{Int,String,String},Int}()
        for row in sort(pdg; by = r -> (r.key, r.mass_MeV))
            rank = get(ranks, row.key, 0)
            ranks[row.key] = rank + 1
            dx = iseven(rank) ? 0.09 : -0.09
            align = iseven(rank) ? (:left, :center) : (:right, :center)
            dy = ((rank % 3) - 1) * 0.008(yhi - ylo)
            text!(ax, xof[row.key] + isospin_offset(row) + dx, row.mass_MeV + dy;
                text = pdg_latex(row.particle), align,
                fontsize = SIMPLIFIED ? 6.8 : 5.4,
                color = (:gray45, 0.62))
        end
    end
    return ax
end

isfile(PDG_TABLE) || error(
    "missing $PDG_TABLE; run `python3 GIPaper/scripts/fetch_pdg_mesons.py` first",
)
mkpath(OUTDIR)

set_theme!(merge(theme_latexfonts(), Theme(fontsize = 13, linewidth = 1.0)))
pdg = CSV.read(PDG_TABLE, DataFrame; stringtype = String)
reuse = lowercase(get(ENV, "GI_REUSE_SPECTRUM", "false")) == "true"
reuse_table = if SIMPLIFIED && isfile(FULL_CALCULATED_TABLE)
    FULL_CALCULATED_TABLE
elseif reuse && isfile(CALCULATED_TABLE)
    CALCULATED_TABLE
else
    nothing
end
if !isnothing(reuse_table)
    println("reusing ", reuse_table)
    model = [merge(NamedTuple(row), (
        flavor_tag = coalesce(row.flavor_tag, ""),
        C = coalesce(row.C, ""),
    )) for row in eachrow(CSV.read(reuse_table, DataFrame; stringtype = String))]
else
    model = calculate_spectra()
    CSV.write(CALCULATED_TABLE, DataFrame(model))
end
filter!(model_row_in_scope, model)

fig = Figure(size = (1900, 900), backgroundcolor = :white)
for (index, class) in enumerate(CLASS_ORDER)
    slot = ((index - 1) ÷ 5 + 1, (index - 1) % 5 + 1)
    model_class = [row for row in model if row.class == class]
    draw_panel!(fig, slot, class, model_class, pdg_rows(pdg, class))
end

Label(fig[:, 0], L"\mathrm{Mass\ (MeV)}"; rotation = pi / 2,
    fontsize = SIMPLIFIED ? 27 : 24)
Label(fig[0, :],
    SIMPLIFIED ?
        "Godfrey–Isgur low-lying meson spectrum — calculated levels, PDG 2026 and recent results" :
        "Godfrey–Isgur meson spectrum — calculated levels, PDG 2026 and recent results";
    fontsize = SIMPLIFIED ? 25 : 23, font = :bold)
Label(fig[3, :],
    SIMPLIFIED ?
        "colored bars: calculated n²ˢ⁺¹Lⱼ     allowed C partners combined     J=3,4: P,C fully merged     ● observed     ○ inferred partner     whiskers: PDG uncertainty" :
        "colored bars: calculated  n²ˢ⁺¹Lⱼ     ●: charge state observed     ○: isospin partner not charge-confirmed     whiskers: PDG uncertainty";
    fontsize = SIMPLIFIED ? 14 : 13)
rowgap!(fig.layout, 8)
colgap!(fig.layout, 6)

png_path = joinpath(OUTDIR, string(OUTPUT_STEM, ".png"))
pdf_path = joinpath(OUTDIR, string(OUTPUT_STEM, ".pdf"))
save(png_path, fig; px_per_unit = 2)
save(pdf_path, fig)
isfile(CALCULATED_TABLE) && println("spectrum table ", CALCULATED_TABLE)
println("wrote ", png_path)
println("wrote ", pdf_path)
