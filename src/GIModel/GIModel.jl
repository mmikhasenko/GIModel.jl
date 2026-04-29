module GIModel

using LinearAlgebra
using Printf
using TOML

export GIParameters,
    ReferenceState,
    load_parameters,
    load_reference_spectrum,
    solve_sector,
    compare_sector,
    write_residual_report,
    parse_quark_masses,
    reduced_mass,
    channel_solution,
    fine_structure_split,
    LdotS,
    tensor_triplet_LJ,
    spin_dot,
    CentralPotentialPath,
    central_potential_path

const ALPHA_COEFFS = (0.25, 0.15, 0.20)
const ALPHA_GAMMAS = (0.5, sqrt(10.0) / 2, sqrt(1000.0) / 2)
const L_SYMBOLS = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)

struct GIParameters
    masses::Dict{String, Float64}
    b::Float64
    c::Float64
    sigma0::Float64
    smearing_s::Float64
    appendix_a_smearing::Bool
    epsilon_c::Float64
    epsilon_t::Float64
    epsilon_so_vector::Float64
    epsilon_so_scalar::Float64
    fine_structure::Bool
    k_spin_orbit::Float64
    k_tensor::Float64
    coulomb_1d_smear::Bool
end

struct ReferenceState
    sector::String
    quark_content::String
    composition::String
    n::Int
    multiplicity::Int
    L::String
    J::Int
    mass_GeV::Float64
    confidence::String
end

function load_parameters(path::AbstractString)
    raw = TOML.parsefile(path)
    masses = Dict(
        "u" => raw["masses"]["m_ud_avg_MeV"] / 1000,
        "d" => raw["masses"]["m_ud_avg_MeV"] / 1000,
        "q" => raw["masses"]["m_ud_avg_MeV"] / 1000,
        "s" => raw["masses"]["m_s_MeV"] / 1000,
        "c" => raw["masses"]["m_c_MeV"] / 1000,
        "b" => raw["masses"]["m_b_MeV"] / 1000,
    )
    rf = get(raw, "relativistic_factors", nothing)
    eps_c = isnothing(rf) ? 0.0 : get(rf, "epsilon_c", 0.0)
    eps_t = isnothing(rf) ? 0.0 : get(rf, "epsilon_t", 0.0)
    eps_v = isnothing(rf) ? 0.0 : get(rf, "epsilon_so_vector", 0.0)
    eps_s = isnothing(rf) ? 0.0 : get(rf, "epsilon_so_scalar", 0.0)
    fs = get(raw, "fine_structure", nothing)
    fine_on = isnothing(fs) ? true : get(fs, "enabled", true)
    k_so = isnothing(fs) ? 0.5 : get(fs, "k_spin_orbit", 0.5)
    k_tn = isnothing(fs) ? 0.4 : get(fs, "k_tensor", 0.4)
    GIParameters(
        masses,
        raw["potential"]["b_GeV2"],
        raw["potential"]["c_MeV"] / 1000,
        raw["relativistic_smearing"]["sigma0_GeV"],
        raw["relativistic_smearing"]["s"],
        get(raw["potential"], "appendix_a_smearing", false),
        float(eps_c),
        float(eps_t),
        float(eps_v),
        float(eps_s),
        fine_on,
        float(k_so),
        float(k_tn),
        get(raw["potential"], "coulomb_1d_smear", false),
    )
end

function parse_csv_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inquote = false
    i = firstindex(line)
    while i <= lastindex(line)
        ch = line[i]
        if ch == '"'
            if inquote && i < lastindex(line) && line[nextind(line, i)] == '"'
                print(buf, '"')
                i = nextind(line, i)
            else
                inquote = !inquote
            end
        elseif ch == ',' && !inquote
            push!(fields, String(take!(buf)))
        else
            print(buf, ch)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    fields
end

function load_reference_spectrum(path::AbstractString)
    lines = readlines(path)
    header = parse_csv_line(first(lines))
    index = Dict(name => i for (i, name) in enumerate(header))
    states = ReferenceState[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        row = parse_csv_line(line)
        push!(
            states,
            ReferenceState(
                row[index["sector"]],
                row[index["quark_content"]],
                row[index["composition_raw"]],
                parse(Int, row[index["n"]]),
                parse(Int, row[index["multiplicity"]]),
                row[index["L"]],
                parse(Int, row[index["J"]]),
                parse(Float64, row[index["mass_GeV"]]),
                row[index["confidence"]],
            ),
        )
    end
    states
end

function erf_approx(x::Real)
    signx = x < 0 ? -1.0 : 1.0
    z = abs(float(x))
    t = 1 / (1 + 0.3275911 * z)
    poly =
        (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t
    signx * (1 - poly * exp(-z^2))
end

function erf_approx_prime(x::Real)
    z = abs(float(x))
    p = 0.3275911
    t = 1 / (1 + p * z)
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    q = a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))
    qp = a2 + t * (2a3 + t * (3a4 + t * (4a5)))
    poly = t * q
    poly_p = q + t * qp
    dt_dz = -p * t^2
    expfac = exp(-z^2)
    expfac * (2z * poly - poly_p * dt_dz)
end

alpha_s_r(r::Real) = sum(a * erf_approx(g * r) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

function central_potential(r::Real, params::GIParameters)
    params.b * r - (4 / 3) * alpha_s_r(r) / r + params.c
end

"""Coulomb piece G(r) = -4 α_s / (3 r) from the text; independent of the linear + constant term S(r) = b r + c."""
function static_coulomb_G(r::Real, params::GIParameters)
    ri = max(float(r), 1.0e-12)
    -(4 / 3) * alpha_s_r(ri) / ri
end

function static_confinement_S(r::Real, params::GIParameters)
    params.b * float(r) + params.c
end

#
# 3D isotropic Gaussian smearing of a spherically symmetric radial function V (|r|),
# Appendix A, Eqs. (A7)–(A8), PDF p. 36. Same σ as contact_smearing_sigma (A9), Table II.
# R > 0: one-dimensional form from the angle-integrated convolution (e.g. difference of
# Gaussians). R → 0: direct radial integral with isotropic 3D Gaussian.
#
function smear_3d_radial(v::AbstractVector{<:Real}, r::AbstractVector{<:Real}, σ::Real)
    n = length(r)
    n == 0 && return eltype(r)[]
    n == 1 && return v
    h = r[2] - r[1]
    w = fill(h, n)
    w[1] = h / 2
    w[n] = h / 2
    σf = max(float(σ), 1.0e-12)
    R0 = 0.25 * h
    out = similar(r, Float64)
    for i in eachindex(r)
        R = r[i]
        s = 0.0
        if R < R0
            for j in eachindex(r)
                rp = r[j]
                # ρ(r) = σ^3 / π^(3/2) exp(-σ^2 r^2), with σ in GeV and r in GeV^-1.
                ρ = σf^3 / (π^(3 / 2)) * exp(-(σf * rp)^2)
                s += 4 * π * rp^2 * w[j] * ρ * v[j]
            end
        else
            for j in eachindex(r)
                rp = r[j]
                # Angle-integrated convolution of a 3D isotropic Gaussian:
                # f̃(R) = (σ / (√π R)) ∫ dr' r' [e^{-σ^2 (R-r')^2} - e^{-σ^2 (R+r')^2}] f(r')
                pre = σf / (sqrt(π) * R)
                s += w[j] * pre * rp * (exp(-(σf * (R - rp))^2) - exp(-(σf * (R + rp))^2)) * v[j]
            end
        end
        out[i] = s
    end
    return out
end

# Experimental. Not the GI (A12)–(A13) smeared potential used in the paper: (A7)–(A8)
# with the Table II width applied to the pointwise G and S (Eqs. (11)–(13) orient.)
# is numerically uncontrolled on a fixed radial line when $\sigma$ is O(1): the 3D
# convolution weights the large-$r$ region by volume and can remove the $1/r$ well.
# Enable only for research; production defaults keep `appendix_a_smearing = false`.
function smeared_central_values(params::GIParameters, m1::Real, m2::Real, r::AbstractVector{<:Real})
    σ = contact_smearing_sigma(params, m1, m2)
    n = length(r)
    n < 2 && return [central_potential(ri, params) for ri in r]
    h = r[2] - r[1]
    rmax0 = r[end]
    # Kernel tail: exp(-(σ Δr)^2) at Δr = 8/σ gives exp(-64), effectively zero.
    n_tail = σ > 0 ? max(0, Int(ceil(8 / (σ * h)))) : 0
    r_ext = n_tail > 0 ? vcat(r, collect(range(rmax0 + h, rmax0 + n_tail * h; step = h))) : r
    g0 = [static_coulomb_G(ri, params) for ri in r_ext]
    s0 = [static_confinement_S(ri, params) for ri in r_ext]
    vsum = smear_3d_radial(g0, r_ext, σ) .+ smear_3d_radial(s0, r_ext, σ)
    return vsum[1:n]
end

include("radial_1d_coulomb_smear.jl")

function reduced_mass(m1::Real, m2::Real)
    m1 * m2 / (m1 + m2)
end

function radial_grid(ngrid::Integer, rmax::Real)
    h = rmax / (ngrid + 1)
    collect(h:h:(ngrid * h)), h
end

function p2_operator(m::Real, L::Integer, r::AbstractVector, h::Real)
    ngrid = length(r)
    diagonal = similar(r)
    offdiag = fill(-1 / h^2, ngrid - 1)
    for i in eachindex(r)
        diagonal[i] = 2 / h^2 + L * (L + 1) / r[i]^2
    end
    SymTridiagonal(diagonal, offdiag)
end

function potential_diagonal(params::GIParameters, m1::Real, m2::Real, r::AbstractVector)
    if params.appendix_a_smearing
        return smeared_central_values(params, m1, m2, r)
    end
    if params.coulomb_1d_smear
        return coulomb_1d_smeared_central_values(params, m1, m2, r)
    end
    return [central_potential(ri, params) for ri in r]
end

function nonrelativistic_hamiltonian(params::GIParameters, m1::Real, m2::Real, L::Integer; ngrid::Integer = 900, rmax::Real = 24.0)
    mu = reduced_mass(m1, m2)
    r, h = radial_grid(ngrid, rmax)
    diagonal = similar(r)
    offdiag = fill(-1 / (2 * mu * h^2), ngrid - 1)
    vdiag = potential_diagonal(params, m1, m2, r)
    for i in eachindex(r)
        ri = r[i]
        diagonal[i] = 1 / (mu * h^2) + L * (L + 1) / (2 * mu * ri^2) + vdiag[i]
    end
    SymTridiagonal(diagonal, offdiag), r
end

function sqrt_kinetic_matrix(p2::SymTridiagonal, m::Real)
    fact = eigen(p2)
    fact.vectors * Diagonal(sqrt.(max.(fact.values, 0) .+ m^2)) * fact.vectors'
end

function relativistic_hamiltonian(params::GIParameters, m1::Real, m2::Real, L::Integer; ngrid::Integer = 450, rmax::Real = 24.0)
    r, h = radial_grid(ngrid, rmax)
    p2 = p2_operator(m1, L, r, h)
    kinetic = sqrt_kinetic_matrix(p2, m1) + sqrt_kinetic_matrix(p2, m2)
    Symmetric(kinetic + Diagonal(potential_diagonal(params, m1, m2, r))), r
end

function channel_solution(
    params::GIParameters,
    m1::Real,
    m2::Real,
    L::Integer;
    nlevels::Integer = 6,
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
)
    hamiltonian, r =
        if kinetic == :relativistic
            relativistic_hamiltonian(params, m1, m2, L; ngrid = ngrid, rmax = rmax)
        elseif kinetic == :nonrelativistic
            nonrelativistic_hamiltonian(params, m1, m2, L; ngrid = ngrid, rmax = rmax)
        else
            error("unknown kinetic mode: $kinetic")
        end
    fact = eigen(hamiltonian)
    values = fact.values
    if kinetic == :relativistic
        values[1:nlevels], fact.vectors[:, 1:nlevels], r
    else
        values[1:nlevels] .+ (m1 + m2), fact.vectors[:, 1:nlevels], r
    end
end

function solve_channel(args...; kwargs...)
    values, _vectors, _r = channel_solution(args...; kwargs...)
    values
end

function contact_smearing_sigma(params::GIParameters, m1::Real, m2::Real)
    mass_factor = 4 * m1 * m2 / (m1 + m2)^2
    reduced_twice = 2 * m1 * m2 / (m1 + m2)
    sqrt(params.sigma0^2 * (0.5 + 0.5 * mass_factor^4) + params.smearing_s^2 * reduced_twice^2)
end

function spin_dot(multiplicity::Integer)
    S = (multiplicity - 1) / 2
    0.5 * (S * (S + 1) - 1.5)
end

function contact_hyperfine_shift(params::GIParameters, m1::Real, m2::Real, L::String, multiplicity::Integer, vector::AbstractVector, r::AbstractVector)
    L == "S" || return 0.0
    multiplicity in (1, 3) || return 0.0
    length(r) >= 2 || return 0.0
    sigma = contact_smearing_sigma(params, m1, m2)
    h = r[2] - r[1]
    expectation = radial_expect_udr(
        vector,
        r,
        h,
        (ri, i) -> begin
            delta_sigma = sigma^3 / (π^(3 / 2)) * exp(-(sigma * ri)^2)
            alpha_s_r(ri) * delta_sigma
        end,
    )
    (1.0 + params.epsilon_c) * (32 * π / (9 * m1 * m2)) * expectation *
    spin_dot(multiplicity)
end

include("masses_from_content.jl")
include("spin_fine_structure.jl")
include("appendix_a_status.jl")

function solve_sector(params::GIParameters, flavor::String; maxn::Integer = 6, ngrid::Integer = 450, rmax::Real = 24.0, kinetic::Symbol = :relativistic)
    m = params.masses[flavor]
    results = Dict{Tuple{Int, String}, Float64}()
    for (symbol, L) in L_SYMBOLS
        levels = solve_channel(params, m, m, L; nlevels = maxn, ngrid = ngrid, rmax = rmax, kinetic = kinetic)
        for n in 1:length(levels)
            results[(n, symbol)] = levels[n]
        end
    end
    results
end

function compare_sector(
    params::GIParameters,
    reference::Vector{ReferenceState},
    flavor::String;
    ngrid::Integer = 450,
    rmax::Real = 24.0,
    kinetic::Symbol = :relativistic,
    contact_hyperfine::Bool = true,
    use_fine_structure::Bool = true,
)
    m_fallback = params.masses[flavor]
    function masses_for(s::ReferenceState)
        try
            return parse_quark_masses(params, String(s.sector), String(s.quark_content))
        catch
            return m_fallback, m_fallback
        end
    end
    channel_cache = Dict{Tuple{Float64, Float64, String}, Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}}}()
    for state in reference
        m1, m2 = masses_for(state)
        cache_key = (round(m1, sigdigits = 12), round(m2, sigdigits = 12), state.L)
        if !haskey(channel_cache, cache_key)
            Lval = L_SYMBOLS[state.L]
            channel_cache[cache_key] = channel_solution(
                params, m1, m2, Lval;
                nlevels = 6, ngrid = ngrid, rmax = rmax, kinetic = kinetic,
            )
        end
    end
    rows = NamedTuple[]
    for state in reference
        m1, m2 = masses_for(state)
        cache_key = (round(m1, sigdigits = 12), round(m2, sigdigits = 12), state.L)
        haskey(channel_cache, cache_key) || continue
        values, vectors, r = channel_cache[cache_key]
        state.n <= length(values) || continue
        h = r[2] - r[1]
        central = values[state.n]
        contact_shift = 0.0
        spin_orbit_shift = 0.0
        tensor_shift = 0.0
        fine_structure_shift = 0.0
        if contact_hyperfine
            contact_shift = contact_hyperfine_shift(params, m1, m2, state.L, state.multiplicity, vectors[:, state.n], r)
        end
        if use_fine_structure && params.fine_structure
            comp = fine_structure_components(
                params, m1, m2, state.L, state.multiplicity, state.J,
                collect(vectors[:, state.n]), collect(r), h;
                enabled = true, k_spin_orbit = params.k_spin_orbit, k_tensor = params.k_tensor,
            )
            spin_orbit_shift = comp.spin_orbit
            tensor_shift = comp.tensor
            fine_structure_shift = comp.total
        end
        predicted = central + contact_shift + fine_structure_shift
        push!(
            rows,
            (
                sector = state.sector,
                state = state.composition,
                L = state.L,
                n = state.n,
                J = state.J,
                multiplicity = state.multiplicity,
                reference_GeV = state.mass_GeV,
                central_GeV = central,
                contact_shift_GeV = contact_shift,
                spin_orbit_shift_GeV = spin_orbit_shift,
                tensor_shift_GeV = tensor_shift,
                fine_structure_shift_GeV = fine_structure_shift,
                predicted_GeV = predicted,
                residual_MeV = 1000 * (predicted - state.mass_GeV),
                confidence = state.confidence,
            ),
        )
    end
    rows
end

function write_residual_report(
    path::AbstractString,
    title::AbstractString,
    rows;
    kinetic::Symbol = :relativistic,
    contact_hyperfine::Bool = true,
    appendix_a_smearing::Bool = false,
    coulomb_1d_smear::Bool = false,
    appendix_a_central::Union{Nothing, Bool} = nothing,
    use_fine_structure::Bool = true,
)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io)
        if !isnothing(appendix_a_central)
            appendix_a_smearing = appendix_a_central
        end
        hyperfine_note = contact_hyperfine ? "with smeared S-wave contact hyperfine" : "without S-wave contact hyperfine"
        fs_note = use_fine_structure ? " first-order L·S (vector+Thomas) and OGE-tensor; " : " no first-order L·S/tensor; "
        central_note = if appendix_a_smearing
            "experimental (A7)–(A8)-style 3D isotropic smearing of pointwise Coulomb G and confinement S (Table II σ₀, s), "
        elseif coulomb_1d_smear
            "1D Gaussian renormalization of G(r) only (pointwise S); same σ as contact (A9); not the full (A12)–(A13) expansion, "
        else
            "pointwise Coulomb + linear + constant (no Appendix A or 1D G smear), "
        end
        println(
            io,
            "Model: finite-difference + `$kinetic` kinetic, $hyperfine_note,$fs_note",
            "GI Table II `b`, `c`, masses, `ε` factors, and Fig. 2 `α_s(r)`;",
            " ",
            central_note,
            "see `src/GIModel/`.",
        )
        println(io)
        println(io, "| state | reference GeV | baseline GeV | residual MeV | confidence |")
        println(io, "|---|---:|---:|---:|---|")
        for row in rows
            label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
            println(
                io,
                @sprintf(
                    "| `%s` | %.3f | %.3f | %+7.1f | %s |",
                    label,
                    row.reference_GeV,
                    row.predicted_GeV,
                    row.residual_MeV,
                    row.confidence,
                ),
            )
        end

        if !isempty(rows) && hasproperty(rows[1], :central_GeV)
            println(io)
            println(io, "## Contribution Breakdown (diagnostic)")
            println(io)
            println(io, "All shifts below are relative to the central FD eigenvalue (the spin-independent Hamiltonian on the current mesh).")
            println(io)
            println(io, "| state | central GeV | contact MeV | L·S MeV | tensor MeV | total shift MeV | predicted GeV |")
            println(io, "|---|---:|---:|---:|---:|---:|---:|")
            for row in rows
                label = @sprintf("%d^%d%s_%d", row.n, row.multiplicity, row.L, row.J)
                total_shift = row.contact_shift_GeV + row.fine_structure_shift_GeV
                println(
                    io,
                    @sprintf(
                        "| `%s` | %.3f | %+7.1f | %+7.1f | %+7.1f | %+7.1f | %.3f |",
                        label,
                        row.central_GeV,
                        1000 * row.contact_shift_GeV,
                        1000 * row.spin_orbit_shift_GeV,
                        1000 * row.tensor_shift_GeV,
                        1000 * total_shift,
                        row.predicted_GeV,
                    ),
                )
            end
        end

        residuals = [abs(row.residual_MeV) for row in rows]
        if !isempty(residuals)
            println(io)
            println(io, @sprintf("Mean absolute residual: %.1f MeV.", sum(residuals) / length(residuals)))
            println(io, @sprintf("Max absolute residual: %.1f MeV.", maximum(residuals)))
        end
        groups = Dict{Tuple{Int, String}, Vector{eltype(rows)}}()
        for row in rows
            push!(get!(groups, (row.n, row.L), eltype(rows)[]), row)
        end
        println(io)
        println(io, "## Spin-Averaged Diagnostics")
        println(io)
        println(io, "Weighted by `2J+1` within each available `(n, L)` group.")
        println(io)
        println(io, "| multiplet | states | reference GeV | baseline GeV | residual MeV |")
        println(io, "|---|---:|---:|---:|---:|")
        for key in sort(collect(keys(groups)); by = x -> (x[2], x[1]))
            group = groups[key]
            weights = [2 * row.J + 1 for row in group]
            weight_sum = sum(weights)
            ref = sum(w * row.reference_GeV for (w, row) in zip(weights, group)) / weight_sum
            pred = sum(w * row.predicted_GeV for (w, row) in zip(weights, group)) / weight_sum
            println(io, @sprintf("| `%d%s` | %d | %.3f | %.3f | %+7.1f |", key[1], key[2], length(group), ref, pred, 1000 * (pred - ref)))
        end
        println(io)
        println(io, "This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.")
    end
end

end
