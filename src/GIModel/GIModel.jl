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
    write_residual_report

const ALPHA_COEFFS = (0.25, 0.15, 0.20)
const ALPHA_GAMMAS = (0.5, sqrt(10.0) / 2, sqrt(1000.0) / 2)
const L_SYMBOLS = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4)

struct GIParameters
    masses::Dict{String, Float64}
    b::Float64
    c::Float64
    sigma0::Float64
    smearing_s::Float64
end

struct ReferenceState
    sector::String
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
    GIParameters(
        masses,
        raw["potential"]["b_GeV2"],
        raw["potential"]["c_MeV"] / 1000,
        raw["relativistic_smearing"]["sigma0_GeV"],
        raw["relativistic_smearing"]["s"],
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

alpha_s_r(r::Real) = sum(a * erf_approx(g * r) for (a, g) in zip(ALPHA_COEFFS, ALPHA_GAMMAS))

function central_potential(r::Real, params::GIParameters)
    params.b * r - (4 / 3) * alpha_s_r(r) / r + params.c
end

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

function potential_diagonal(params::GIParameters, r::AbstractVector)
    [central_potential(ri, params) for ri in r]
end

function nonrelativistic_hamiltonian(params::GIParameters, m1::Real, m2::Real, L::Integer; ngrid::Integer = 900, rmax::Real = 24.0)
    mu = reduced_mass(m1, m2)
    r, h = radial_grid(ngrid, rmax)
    diagonal = similar(r)
    offdiag = fill(-1 / (2 * mu * h^2), ngrid - 1)
    for i in eachindex(r)
        ri = r[i]
        diagonal[i] = 1 / (mu * h^2) + L * (L + 1) / (2 * mu * ri^2) + central_potential(ri, params)
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
    Symmetric(kinetic + Diagonal(potential_diagonal(params, r))), r
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
    sigma = contact_smearing_sigma(params, m1, m2)
    expectation = 0.0
    for i in eachindex(r)
        delta_sigma = sigma^3 / (π^(3 / 2)) * exp(-(sigma * r[i])^2)
        expectation += abs2(vector[i]) * alpha_s_r(r[i]) * delta_sigma
    end
    (32 * π / (9 * m1 * m2)) * expectation * spin_dot(multiplicity)
end

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
)
    m = params.masses[flavor]
    channel_cache = Dict{String, Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}}}()
    for symbol in keys(L_SYMBOLS)
        channel_cache[symbol] = channel_solution(params, m, m, L_SYMBOLS[symbol]; nlevels = 6, ngrid = ngrid, rmax = rmax, kinetic = kinetic)
    end
    rows = NamedTuple[]
    for state in reference
        haskey(channel_cache, state.L) || continue
        values, vectors, r = channel_cache[state.L]
        state.n <= length(values) || continue
        predicted = values[state.n]
        if contact_hyperfine
            predicted += contact_hyperfine_shift(params, m, m, state.L, state.multiplicity, vectors[:, state.n], r)
        end
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
                predicted_GeV = predicted,
                residual_MeV = 1000 * (predicted - state.mass_GeV),
                confidence = state.confidence,
            ),
        )
    end
    rows
end

function write_residual_report(path::AbstractString, title::AbstractString, rows; kinetic::Symbol = :relativistic, contact_hyperfine::Bool = true)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ", title)
        println(io)
        hyperfine_note = contact_hyperfine ? "with smeared S-wave contact hyperfine" : "without spin-dependent terms"
        println(io, "Baseline model: radial finite-difference solver with `$kinetic` kinetic energy, $hyperfine_note, GI Table II masses, `b`, `c`, and Fig. 2 running Coulomb ansatz.")
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
        residuals = [abs(row.residual_MeV) for row in rows]
        if !isempty(residuals)
            println(io)
            println(io, @sprintf("Mean absolute residual: %.1f MeV.", sum(residuals) / length(residuals)))
            println(io, @sprintf("Max absolute residual: %.1f MeV.", maximum(residuals)))
        end
        println(io)
        println(io, "This is a diagnostic baseline, not the final GI Hamiltonian. Large residuals are expected until the full smeared potential, tensor/spin-orbit terms, and mixing are added.")
    end
end

end
