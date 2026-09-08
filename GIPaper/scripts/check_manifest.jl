#!/usr/bin/env julia
# =============================================================================
# Paper-manifest checker — the anti-drift gate
# =============================================================================
# Globs docs/paper_manifest/*.toml, validates that every cited code symbol /
# test / report / page-image actually exists, rolls status up per section, and
# writes docs/residual_reports/manifest_check.md. Exits nonzero on any broken
# link so it can gate CI (scripts/verify_project.sh).
#
# Uses only Julia stdlib (TOML, Printf) — no project activation needed.

using TOML
using Printf

const ROOT = dirname(dirname(dirname(abspath(@__FILE__))))   # scripts -> GIPaper -> repo
const MANIFEST_DIR = joinpath(ROOT, "docs", "paper_manifest")
const REPORT = joinpath(ROOT, "GIPaper", "docs", "residual_reports", "manifest_check.md")
const TEST_DIRS = [
    joinpath(ROOT, "test"),
    joinpath(ROOT, "GIPaper", "test"),
]

const STATUSES = Set(["reproduced", "partial", "implemented", "folded", "context", "missing", "todo"])
const KINDS = Set(["equation", "table", "figure", "section", "input", "fit"])
const EFFORTS = Set(["small", "medium", "large", "research"])
# statuses that count as open work (surfaced in the derived "Remaining work" list)
const OPEN_STATUS = ("missing", "partial", "todo")
const EFFORT_RANK = Dict("large" => 1, "research" => 2, "medium" => 3, "small" => 4)

# repo-relative path -> absolute; check existence
exists_rel(p) = isfile(joinpath(ROOT, p))

# "file:symbol" or "file" -> (file, symbol|nothing)
function split_code_ref(ref)
    parts = split(ref, ':')
    length(parts) == 1 ? (parts[1], nothing) : (join(parts[1:end-1], ':'), parts[end])
end

function load_units()
    units = Dict{String,Any}[]
    for f in sort(readdir(MANIFEST_DIR))
        endswith(f, ".toml") || continue
        data = TOML.parsefile(joinpath(MANIFEST_DIR, f))
        for u in get(data, "unit", Any[])
            u["_file"] = f
            push!(units, u)
        end
    end
    return units
end

function main()
    units = load_units()
    ids = Set(u["id"] for u in units)
    test_sources = join(
        (
            read(joinpath(dir, file), String)
            for dir in TEST_DIRS if isdir(dir)
            for file in sort(readdir(dir)) if endswith(file, ".jl")
        ),
        "\n",
    )
    # cache file contents for symbol checks
    filecache = Dict{String,String}()
    readcached(p) = get!(filecache, p) do
        exists_rel(p) ? read(joinpath(ROOT, p), String) : ""
    end

    problems = Tuple{String,String}[]   # (unit id, message)
    seen = Set{String}()
    for u in units
        id = u["id"]
        id in seen && push!(problems, (id, "duplicate id (also in another fragment)"))
        push!(seen, id)
        get(u, "status", "") in STATUSES || push!(problems, (id, "unknown status `$(get(u,"status",""))`"))
        get(u, "kind", "") in KINDS || push!(problems, (id, "unknown kind `$(get(u,"kind",""))`"))
        eff = get(u, "effort", "")
        isempty(eff) || eff in EFFORTS || push!(problems, (id, "unknown effort `$eff` (small|medium|large|research)"))

        # code[]: file must exist; symbol (if given) must occur in it
        for ref in get(u, "code", String[])
            file, sym = split_code_ref(ref)
            if !exists_rel(String(file))
                push!(problems, (id, "code file missing: `$file`"))
            elseif sym !== nothing && !occursin(String(sym), readcached(String(file)))
                push!(problems, (id, "symbol `$sym` not found in `$file`"))
            end
        end
        # tests[]: name must occur in runtests.jl
        for t in get(u, "tests", String[])
            occursin(String(t), test_sources) || push!(
                problems,
                (id, "test `$t` not found under core or GIPaper test/"),
            )
        end
        # report / anchor: file must exist (if given)
        for key in ("report", "anchor")
            p = get(u, key, "")
            isempty(p) || exists_rel(String(p)) || push!(problems, (id, "$key missing: `$p`"))
        end
        # pages: multipage span "A-B" -> every page image in the range must exist
        pg = get(u, "pages", "")
        if !isempty(pg)
            m = match(r"^(\d+)-(\d+)$", String(pg))
            if m === nothing
                push!(problems, (id, "malformed pages `$pg` (expected \"A-B\")"))
            else
                a, b = parse(Int, m.captures[1]), parse(Int, m.captures[2])
                a <= b || push!(problems, (id, "pages `$pg` is not ascending"))
                for n in a:b
                    img = "paper/vision_ocr/page_images/" * @sprintf("page-%03d.png", n)
                    exists_rel(img) || push!(problems, (id, "pages span image missing: `$img`"))
                end
            end
        end
        # parent must be "" or an existing id
        par = get(u, "parent", "")
        (isempty(par) || par in ids) || push!(problems, (id, "parent `$par` is not a known unit id"))
    end

    # --- rollup by section ---------------------------------------------------
    sections = [u for u in units if get(u, "kind", "") == "section"]
    childstatus(secid) = [get(u, "status", "todo") for u in units if get(u, "parent", "") == secid]
    # Derived section status = the weakest SUBSTANTIVE child. `folded`/`context`
    # are "done by other means / not-applicable" and don't drag a section down;
    # a section with only those rolls up to folded (or context if purely
    # discussion), never to a green "reproduced".
    function derived(statuses)
        isempty(statuses) && return "todo"
        core = filter(s -> s in ("reproduced", "implemented", "partial", "missing", "todo"), statuses)
        isempty(core) && return all(==("context"), statuses) ? "context" : "folded"
        all(s -> s in ("missing", "todo"), core) && return "todo"
        any(s -> s in ("missing", "todo", "partial"), core) && return "partial"
        any(==("implemented"), core) && return "implemented"
        return "reproduced"
    end

    open(REPORT, "w") do io
        println(io, "# Paper Manifest Check")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/check_manifest.jl`. Validates that every")
        println(io, "cited code symbol / test / report / page image in `docs/paper_manifest/*.toml`")
        println(io, "resolves, and rolls status up per section.")
        println(io)
        println(io, @sprintf("**%d units across %d fragments; %d broken links.**",
            length(units), length(unique(u["_file"] for u in units)), length(problems)))
        println(io)
        println(io, "## Broken links")
        println(io)
        if isempty(problems)
            println(io, "None — every reference resolves.")
        else
            for (id, msg) in problems
                println(io, "- `", id, "`: ", msg)
            end
        end
        println(io)
        println(io, "## Section rollup")
        println(io)
        println(io, "| section | units | reproduced | partial | implemented | folded | context | todo/missing | derived |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|:-:|")
        cnt(ss, s) = count(==(s), ss)
        for sec in sort(sections, by = s -> s["id"])
            ss = childstatus(sec["id"])
            println(io, @sprintf("| %s | %d | %d | %d | %d | %d | %d | %d | %s |",
                sec["label"], length(ss), cnt(ss, "reproduced"), cnt(ss, "partial"),
                cnt(ss, "implemented"), cnt(ss, "folded"), cnt(ss, "context"),
                cnt(ss, "todo") + cnt(ss, "missing"), derived(ss)))
        end

        # --- derived "Remaining work" list (projection of open statuses) -----
        seclabel = Dict(s["id"] => get(s, "label", s["id"]) for s in sections)
        openu = [u for u in units if get(u, "kind", "") != "section" &&
                                     get(u, "status", "") in OPEN_STATUS]
        sort!(openu, by = u -> (get(EFFORT_RANK, get(u, "effort", ""), 9),
                                get(u, "status", ""), u["id"]))
        println(io)
        println(io, "## Remaining work")
        println(io)
        println(io, "Projected from every `missing`/`partial`/`todo` unit — not a separate list, ",
                    "so it cannot drift from the statuses above.")
        println(io)
        if isempty(openu)
            println(io, "None — every unit is reproduced/implemented/folded/context.")
        else
            println(io, "| unit | section | status | effort | next step |")
            println(io, "|---|---|:-:|:-:|---|")
            for u in openu
                println(io, @sprintf("| %s | %s | %s | %s | %s |",
                    get(u, "label", u["id"]), get(seclabel, get(u, "parent", ""), ""),
                    get(u, "status", ""), get(u, "effort", "—"),
                    get(u, "next", get(u, "statement", ""))))
            end
        end
    end

    @printf("manifest: %d units, %d fragments, %d broken links -> %s\n",
        length(units), length(unique(u["_file"] for u in units)), length(problems), REPORT)
    for (id, msg) in problems
        println("  BROKEN  ", id, ": ", msg)
    end
    exit(isempty(problems) ? 0 : 1)
end

main()
