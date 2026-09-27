include(joinpath(@__DIR__, "common.jl"))

# Table IV prints reduced partial-wave amplitude forms. Their calibration
# amplitudes A and S0 are fitted by the Table V reproduction; regenerate that
# report first (scripts/reproduce_table_v.jl, or trace_rate_inputs.jl).
function report_number(text, pattern)
    match_result = match(pattern, text)
    isnothing(match_result) && error("calibration value not found in Table V report")
    return parse(Float64, match_result.captures[1])
end

function compute_table_iv()
    report = read(joinpath(RESIDUAL_REPORTS_DIR, "table_v_reproduction.md"), String)
    A = report_number(report, r"`A = ([0-9.]+)`")
    S0 = report_number(report, r"`S0 = ([0-9.]+)`")
    qbar, heavy_fraction = 1.0, 0.881
    rows = [
        ("A, A', A'', A0", "A", A),
        ("S", "S0 - (1/2) A qbar^2", S0 - 0.5A*qbar^2),
        ("D", "S0 - (3/10) A qbar^2", S0 - 0.3A*qbar^2),
        ("P", "S0 - (3/4) A qbar^2", S0 - 0.75A*qbar^2),
        ("A_c", "A", A),
        ("S_c", "S0 - r A qbar^2", S0 - heavy_fraction*A*qbar^2),
    ]
    lines = ["class\tformula\tcomputed_at_qbar_1\tA\tS0\tr"]
    append!(lines, [join((row..., A, S0, heavy_fraction), '\t') for row in rows])
    return write_paper_table("table_iv.tsv", join(lines, '\n'))
end

abspath(PROGRAM_FILE) == (@__FILE__) && compute_table_iv()
