"""
    load_table_policy([path])

Read GIPaper's table-specific prescriptions and comparison-only validation anchors.
Units, source and role travel with the values. These are separate from generic
QuarkModelTransitions operator defaults and from the experimental mass registry.
"""
function load_table_policy(path::AbstractString = joinpath(paper_data_dir(), "table_policy.toml"))
    policy = TOML.parsefile(path)
    for table in ("table_v", "table_vi", "table_vii")
        haskey(policy, table) || throw(ArgumentError("missing policy section $table"))
    end
    seen = Set{String}()
    for row in policy["table_v"]["mixing"]
        row["section"] in seen && throw(ArgumentError("duplicate mixing section $(row["section"])"))
        push!(seen, row["section"])
        isfinite(row["theta_deg"]) || throw(ArgumentError("nonfinite mixing angle"))
        row["singlet"] != row["triplet"] || throw(ArgumentError("mixing partners must differ"))
        isempty(row["source"]) && throw(ArgumentError("missing mixing provenance"))
    end
    for key in ("rho_amplitude", "B_amplitude")
        isfinite(policy["table_v"][key]) || throw(ArgumentError("nonfinite calibration amplitude"))
    end
    policy["table_v"]["convention"] == "LeadingS0" || throw(ArgumentError("unsupported Table V convention"))
    isfinite(policy["table_vi"]["supplementary_mu_N"]) || throw(ArgumentError("nonfinite supplementary moment"))
    seen = Set{String}()
    for row in policy["table_vii"]["dilepton_validation"]
        row["label"] in seen && throw(ArgumentError("duplicate dilepton validation label"))
        push!(seen, row["label"])
        isfinite(row["width_GeV"]) && row["width_GeV"] > 0 || throw(ArgumentError("invalid validation width"))
    end
    return policy
end
