#!/usr/bin/env julia
# Execute the existing table pipelines and preserve their numerical inputs.
# This is an audit adapter: it does not implement a second rate formula.
using Pkg
Pkg.activate(@__DIR__; io=devnull)
using TOML, SHA, CSV, GIModel, QuarkModelTransitions, GIPaper
using GIModel.QuarkModelTransitions: ALPHA_EM, ELECTROMAGNETIC_DEFAULTS, G_FERMI_GEV, HBARC_FM2, LeadingS0, NUCLEON_MASS_GEV, STRONG_DECAY_DEFAULTS

const TRACE_ROOT = dirname(dirname(@__DIR__))
const TRACE_DIR = joinpath(dirname(@__DIR__), "docs", "input_traces")

function fingerprint_sources()
    files = String[]
    for dir in ("src", "QuarkModelTransitions/src", "GIPaper/src", "GIPaper/data/mass_inputs")
        for (base, _, names) in walkdir(joinpath(TRACE_ROOT,dir)), name in names
            push!(files, relpath(joinpath(base,name),TRACE_ROOT))
        end
    end
    append!(files, ["GIPaper/data/table_policy.toml", "data/parameters.provisional.toml", "GIPaper/data/table_vi_states.csv",
        "GIPaper/scripts/trace_rate_inputs.jl", "GIPaper/scripts/reproduce_table_v.jl",
        "GIPaper/scripts/audit_table_vi_photon_decays.jl", "GIPaper/scripts/audit_table_vii.jl",
        "GIPaper/data/raw/digitized_tables/table_v_strong_decays.csv",
        "GIPaper/data/raw/digitized_tables/table_vi_photon_decays.csv",
        "GIPaper/data/raw/digitized_tables/table_vii_annihilation_em.csv"])
    return Dict(f=>bytes2hex(sha256(read(joinpath(TRACE_ROOT,f)))) for f in sort(files))
end
const INITIAL_SOURCE_HASHES = fingerprint_sources()

trace_value(x::Union{Nothing,Missing}) = "unavailable"
trace_value(x::Symbol) = String(x)
trace_value(x::AbstractString) = String(x)
trace_value(x::Union{Real,Bool}) = x
trace_value(x::Tuple) = [trace_value(v) for v in x]
trace_value(x::AbstractVector) = [trace_value(v) for v in x]
trace_value(x::AbstractDict) = Dict(string(k)=>trace_value(v) for (k,v) in x)
trace_value(x::NamedTuple) = Dict(string(k)=>trace_value(v) for (k,v) in pairs(x))
trace_value(x) = Dict(string(k)=>trace_value(getfield(x,k)) for k in fieldnames(typeof(x)))

function write_trace(name, data)
    mkpath(TRACE_DIR)
    path = joinpath(TRACE_DIR, name * ".toml")
    open(path, "w") do io
        TOML.print(io, trace_value(data); sorted=true)
    end
    println("wrote input trace ", path)
end

mass_trace(context, label) = trace_value(GIPaper.mass_input(context, label))

module TableV
include("reproduce_table_v.jl")
end
module TableVII
include("audit_table_vii.jl")
end
# This script has an intentional top-level pipeline and writes its normal report.
module TableVI
include("audit_table_vi_photon_decays.jl")
end

function table_v_trace(run)
    rows = map(enumerate(CSV.File(TableV.TABLE))) do (index, row)
        section, parent = String(row.section), String(row.parent)
        context = parent == "delta" && section == "1^3F_4" ? "V:1^3F_4" : "V"
        channels = GIPaper.load_table_v((row,);
            heavy_fraction_for = r -> TableV.heavy_fraction_of(String(r.parent)))
        mix = get(TableV.MIX, section, nothing)
        Dict(
            "row"=>index, "decay"=>String(row.decay), "section"=>section,
            "parent"=>mass_trace(context, parent),
            "daughter1"=>mass_trace("V", String(row.daughter1)),
            "daughter2"=>mass_trace("V", String(row.daughter2)),
            "operator"=>isempty(channels) ? "not implemented for this row" : trace_value(only(channels)),
            "external_mixing"=>isnothing(mix) ? "none" : trace_value(mix),
            "comparison_only"=>Dict("predicted"=>string(row.amp_MeV),
                "experiment"=>string(row.experiment), "realistic_factor"=>string(row.realistic_factor)),
            "formula_reference"=>string(row.formula), "footnotes"=>string(row.footnotes),
        )
    end
    Dict("common_inputs"=>"common.toml", "convention"=>"LeadingS0",
        "spectrum_wavefunctions_used"=>false,
        "calibration"=>Dict("paper_policy"=>TableV.TABLE_POLICY, "operator_defaults"=>STRONG_DECAY_DEFAULTS,
            "rho_q_GeV"=>run.q_rho, "B_q_GeV"=>run.q_B,
            "mass_inputs"=>[mass_trace("V", s) for s in ("rho","pi","B","omega")]),
        "derived_model"=>run.model, "canonical_rows"=>rows, "results"=>run.rows,
        "mixing_rule"=>"Q1=cos(theta)*singlet+sin(theta)*triplet; Q2=-sin(theta)*singlet+cos(theta)*triplet; missing partner=0; evaluate both at the physical parent's momentum")
end

function component_trace(component)
    w = component.wave
    Dict("basis"=>trace_value(component.basis), "coefficient"=>component.coefficient,
        "wave"=>w isa OscillatorWave ?
            Dict("representation"=>"HO", "L"=>w.L, "beta_GeV"=>w.beta,
                "coefficients"=>w.coefficients) :
            Dict("representation"=>"FD", "r_GeV_inv"=>w.r, "u"=>w.u))
end

function table_vi_trace()
    rows = map(TableVI.results) do result
        ref = result.reference
        multipole, parent, daughter = ref.id
        ps, p = TableVI.resolved(parent); ds, d = TableVI.resolved(daughter)
        pc = physical_components(ps,p); dc = physical_components(ds,d)
        Dict("id"=>ref.id, "decay"=>ref.decay, "computed"=>result.computed,
            "q_GeV"=>result.q_GeV, "parent_mass"=>mass_trace("VI", parent),
            "daughter_mass"=>mass_trace("VI", daughter),
            "parent_components"=>component_trace.(pc), "daughter_components"=>component_trace.(dc),
            "magnetic_exponent"=>multipole == :M1 ? ELECTROMAGNETIC_DEFAULTS.magnetic_exponent : "unused",
            "electric_exponent"=>multipole != :M1 || "c" in split(ref.footnotes,",") ? ELECTROMAGNETIC_DEFAULTS.electric_exponent : "unused",
            "recoil_c"=>"c" in split(ref.footnotes,","),
            "recoil_g_beta_GeV"=>"g" in split(ref.footnotes,",") ? ELECTROMAGNETIC_DEFAULTS.recoil_beta_GeV : "unused",
            "supplementary_mu_N"=>result.supplementary,
            "isovector_parent"=>parent in TableVI.ISOVECTORS,
            "isovector_daughter"=>daughter in TableVI.ISOVECTORS,
            "formula_reference"=>ref.formula, "footnotes"=>ref.footnotes,
            "comparison_only"=>Dict("paper"=>ref.predicted,"historical_q_GeV"=>result.q_historical_GeV))
    end
    Dict("common_inputs"=>"common.toml", "solver"=>TableVI.solver,
        "isoscalar_policy"=>Dict("pseudoscalar"=>"PaperP1Annihilation",
            "pseudoscalar_basis"=>TableVI.p1_basis,
            "smearing"=>PaperP1Annihilation().smearing,
            "vector_amplitude"=>TableVI.params.annihilation.s1_A,
            "tensor_amplitude"=>TableVI.params.annihilation.a_3p2,
            "annihilation_radial_levels"=>[1]), "rows"=>rows)
end

function table_vii_trace(run)
    function record(config, result, profile)
        label = config.decay
        Dict("row_configuration"=>trace_value(config), "result"=>trace_value(result),
            "mass"=>profile == "charge_radius" ? "not used" : mass_trace("VII",label),
            "wave_profile"=>profile)
    end
    glu = [record(c,r,"fixed sector, at least four radial levels") for (c,r) in zip(TableVII.GLUONIC_ROWS,run.glu)]
    lep = [record(c,r,c.kind in (:Vp_V,:Pp_A1) ? "central, four radial levels" : "fixed S sector, four radial levels") for (c,r) in zip(TableVII.LEPTONIC_ROWS,run.lep)]
    gg = [record(c,r,c.kind == :P ? "fixed S singlet, four radial levels" : "central P, four radial levels; ideal flavor") for (c,r) in zip(TableVII.TWO_PHOTON_ROWS,run.gg)]
    cr = [record(c,r,"charge_radius") for (c,r) in zip(TableVII.CHARGE_RADIUS_ROWS,run.cr)]
    mixed = [Dict("label"=>r.label,"mass"=>mass_trace("VII",r.label),
        "result"=>trace_value(r), "basis"=>["1qq","1ss","2qq","2ss"],
        "q_eff_qq"=>TableVII.QNS,"q_eff_ss"=>1/9,"annihilation"=>"PaperP1Annihilation") for r in run.ggm]
    # Keep every canonical row visible, including unsupported top states.
    implemented = Set(r.label for results in (run.glu, run.lep, run.gg, run.ggm, run.cr) for r in results)
    canonical = [Dict("subtable"=>String(r.subtable), "decay"=>String(r.decay),
        "implemented"=>String(r.decay) in implemented,
        "comparison_only"=>Dict("predicted"=>string(r.predicted)),
        "status"=>String(r.decay) in implemented ? "see result; may lack experimental mass" : "no constituent top mass")
        for r in CSV.File(TableVII.TABLE)]
    for (configs, results) in ((TableVII.GLUONIC_ROWS,run.glu), (TableVII.LEPTONIC_ROWS,run.lep),
                              (TableVII.TWO_PHOTON_ROWS,run.gg), (TableVII.GG_MIXED,run.ggm),
                              (TableVII.CHARGE_RADIUS_ROWS,run.cr))
        @assert length(configs) == length(results)
    end
    Dict("common_inputs"=>"common.toml", "canonical_rows"=>canonical, "solver"=>run.solver_ho,
        "observable_npoints"=>TableVII.NPTS,
        "gluonic"=>glu,"leptonic"=>lep,"two_photon"=>gg,"mixed_two_photon"=>mixed,"charge_radius"=>cr,
        "derived_dilepton_widths"=>run.wid,
        "experimental_validation_only"=>Dict("anchors"=>TableVII.VALIDATION_INPUTS,
            "dilepton_anchors"=>TableVII.DILEPTON_VALIDATION_ROWS,
            "muon"=>mass_trace("D7","mu"), "wpi_from_measured_f_GeV"=>run.wpi_exp,
            "wpi_from_lifetime_GeV"=>run.wpi_pdg,"fpsi_from_measured_width"=>run.fpsi_exp),
        "not_computed"=>"hypothetical top rows: no constituent top mass; registered missing experimental masses remain unavailable")
end

function main()
    v = TableV.main()
    vii = TableVII.main()
    fingerprint_sources() == INITIAL_SOURCE_HASHES || error("Inputs changed during trace generation; rerun after concurrent edits settle")
    common = Dict(
        "schema_version"=>1, "julia_version"=>string(VERSION),
        "paper_table_policy"=>GIPaper.load_table_policy(),
        "spectrum_configuration"=>TOML.parsefile(default_parameters_path()),
        "running_coupling"=>Dict("weights"=>GIModel.ALPHA_COEFFS,"gammas_GeV"=>GIModel.ALPHA_GAMMAS,
            "historical_Lambda_MeV_not_used"=>200),
        "electromagnetic_defaults"=>ELECTROMAGNETIC_DEFAULTS,
        "constants"=>Dict("alpha_em"=>ALPHA_EM,"fermi_GeV_minus2"=>G_FERMI_GEV,
            "hbarc_squared_GeV2_fm2"=>HBARC_FM2,"nucleon_mass_GeV"=>NUCLEON_MASS_GEV,
            "up_charge"=>2/3,"down_charge"=>-1/3),
        "numerics"=>Dict("HO_operator_rtol"=>1e-10,"HO_narrow_spin_rtol"=>1e-8,
            "HO_operator_initial_order"=>64,"HO_operator_max_order"=>8192,
            "HO_momentum_quadgk_rtol"=>1e-9,"FD_default_momentum_pmax_GeV"=>30.0,
            "FD_default_momentum_npoints"=>1501,"FD_observable_pmax_cap_GeV"=>60.0),
        "source_sha256"=>INITIAL_SOURCE_HASHES)
    write_trace("table_v", table_v_trace(v))
    write_trace("table_vi", table_vi_trace())
    write_trace("table_vii", table_vii_trace(vii))
    common["trace_sha256"] = Dict(name=>bytes2hex(sha256(read(joinpath(TRACE_DIR,name))))
        for name in ("table_v.toml", "table_vi.toml", "table_vii.toml"))
    write_trace("common", common)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
