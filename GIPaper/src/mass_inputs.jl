for file in ("pdg-2026.csv", "assignments.csv", "averages.csv")
    Base.include_dependency(joinpath(paper_data_dir(),"mass_inputs",file))
end

"""Pinned experimental kinematic inputs, separate from GI model eigenvalues.

Unknown aliases are errors. A registered but unassigned state returns `nothing`;
there is deliberately no model-mass or historical-mass fallback.
"""
function load_mass_inputs()
    dir = joinpath(paper_data_dir(), "mass_inputs")
    masses = Dict{String,NamedTuple}()
    for r in CSV.File(joinpath(dir, "pdg-2026.csv"); types=String)
        masses[r.key] = (key=String(r.key), particle=String(r.particle),
            mass_GeV=parse(Float64,r.mass_GeV), pdg_id=String(r.pdg_id),
            pdg_display=String(r.pdg_display), edition=String(r.edition),
            status=String(r.status), note=ismissing(r.note) ? "" : String(r.note),
            source_url=String(r.source_url))
    end
    averages = Dict(String(r.key)=>split(String(r.components),";")
        for r in CSV.File(joinpath(dir,"averages.csv");types=String))
    for (key, parts) in averages
        masses[key] = (key=key, particle="(" * join(parts," + ") * ") / $(length(parts))",
            mass_GeV=sum(masses[p].mass_GeV for p in parts)/length(parts),
            pdg_id=join([masses[p].pdg_id for p in parts],";"), pdg_display="",
            edition="2026",status="derived experimental input",
            note="Arithmetic isospin average for charge-unspecified channels; components: " * join(parts,", "),
            source_url="https://pdg.lbl.gov/2026/api/index.html")
    end
    assignments = Dict{Tuple{String,String},NamedTuple}()
    for r in CSV.File(joinpath(dir,"assignments.csv");types=String)
        id=(String(r.context),String(r.label))
        haskey(assignments,id) && error("Duplicate mass assignment $id")
        key=ismissing(r.key) ? "" : String(r.key)
        note=ismissing(r.note) ? "" : String(r.note)
        entry=isempty(key) ? (key="",particle="unassigned",mass_GeV=nothing,
            pdg_id="",pdg_display="",edition="2026",status="unavailable",note="",source_url="") : masses[key]
        assignments[id]=merge(entry,(context=id[1],label=id[2],assignment_note=note,))
    end
    return assignments
end

const _MASS_INPUTS = load_mass_inputs()
mass_input(context::AbstractString,label::AbstractString) = _MASS_INPUTS[(String(context),String(label))]
experimental_mass(context::AbstractString,label::AbstractString) = mass_input(context,label).mass_GeV

"""Archived audit mass for comparison only; never an experimental fallback."""
function historical_mass(context::AbstractString, label::AbstractString)
    matches = [r for r in CSV.File(joinpath(paper_data_dir(),"mass_inputs","historical.csv"))
        if r.context == context && r.label == label]
    isempty(matches) && return nothing
    length(matches)==1 || error("Ambiguous historical mass for $context/$label")
    return Float64(only(matches).mass_GeV)
end
