# Adapter from the paper's canonical Table V schema to GIModel decay channels.
# The CSV schema and its provenance belong to GIPaper; the decay algebra and
# channel type belong to GIModel.

const TABLE_V_AMPLITUDE_CLASSES = Dict(
    "A" => :A,
    "Aprime" => :Aprime,
    "Adoubleprime" => :Adoubleprime,
    "A0" => :A0,
    "S" => :S,
    "D" => :D,
    "P" => :P,
    "Ac" => :A_c,
    "Sc" => :S_c,
)

"""
    load_table_v(rows; heavy_fraction_for = _ -> nothing) -> Vector{DecayChannel}

Convert rows in GIPaper's canonical `table_v_strong_decays.csv` schema to
GIModel [`DecayChannel`](@ref)s. Rows marked `mixing_only` or `unlisted` are
skipped.

`heavy_fraction_for(row)` supplies `r = m_Q/(m_Q + m_q)`. The paper table does
not encode quark content, so callers loading unequal-mass `Ac`/`Sc` rows must
provide this resolver; GIModel will reject those channels when it is absent.
"""
function load_table_v(rows; heavy_fraction_for = _ -> nothing)
    channels = DecayChannel[]
    for row in rows
        class_name = String(row.amp_class)
        haskey(TABLE_V_AMPLITUDE_CLASSES, class_name) || continue
        push!(channels, DecayChannel(
            String(row.parent),
            String(row.daughter1),
            String(row.daughter2),
            Float64(row.coefficient),
            TABLE_V_AMPLITUDE_CLASSES[class_name],
            Int(row.qbar_power);
            label = String(row.decay),
            section = String(row.section),
            heavy_fraction = heavy_fraction_for(row),
        ))
    end
    return channels
end
