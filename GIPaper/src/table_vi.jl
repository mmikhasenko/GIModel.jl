"""
    load_table_vi([path])

Read canonical Table VI references. Identity is the explicit
(multipole, parent, daughter) tuple, independent of CSV order. Preserve the
printed prediction alongside its numeric value and approximation qualifiers.
Formula strings are transcriptions, not executable expressions.
"""
function load_table_vi(path::AbstractString = joinpath(
    paper_data_dir(), "raw", "digitized_tables", "table_vi_photon_decays.csv"))
    rows = map(CSV.File(path; types = String)) do row
        printed = String(row.predicted)
        parenthetical = startswith(printed, "(")
        approximate = occursin("~", printed)
        number = parse(Float64, replace(printed, "(" => "", ")" => "", "~" => ""))
        (
            id = (Symbol(row.multipole), String(row.parent), String(row.daughter)),
            decay = String(row.decay),
            predicted = number,
            predicted_text = printed,
            parenthetical = parenthetical,
            approximate = approximate,
            formula = String(row.formula),
            page = parse(Int, row.page),
            footnotes = ismissing(row.footnotes) ? "" : String(row.footnotes),
        )
    end
    length(unique(r.id for r in rows)) == length(rows) ||
        throw(ArgumentError("Duplicate Table VI transition identity"))
    length(unique(r.decay for r in rows)) == length(rows) ||
        throw(ArgumentError("Duplicate Table VI decay label"))
    return rows
end

"""
    load_table_vi_states()

Canonical spectroscopic identities and explicit photon-kinematics inputs for
Table VI. Kinematics use the pinned PDG registry. Missing masses mean an
experimental assignment is unavailable; historical inputs are archival only.
The `mixed` flag selects a composed isoscalar state; other entries use their
native fixed-channel state. Paper state names are labels, not inferred quark
charges (in particular the historical B*- label is retained).
"""
function load_table_vi_states(path::AbstractString = joinpath(
    paper_data_dir(), "table_vi_states.csv"))
    rows = map(CSV.File(path; types = String)) do row
        (
            id = String(row.id),
            flavors = (Symbol(row.flavor1), Symbol(row.flavor2)),
            basis = BasisState(parse(Int, row.n), row.L,
                parse(Int, row.multiplicity), parse(Int, row.J);
                flavors = (Symbol(row.flavor1), Symbol(row.flavor2))),
            mass_GeV = experimental_mass("VI", String(row.id)),
            historical_mass_GeV = historical_mass("VI", String(row.id)),
            mass_source = mass_input("VI", String(row.id)).status,
            historical_mass_source = ismissing(row.mass_source) ? "" : String(row.mass_source),
            mixed = parse(Bool, row.mixed),
        )
    end
    length(unique(r.id for r in rows)) == length(rows) ||
        throw(ArgumentError("duplicate Table VI state"))
    return Dict(r.id => r for r in rows)
end
