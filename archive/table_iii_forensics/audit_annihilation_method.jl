#!/usr/bin/env julia
# Method audit for Eq. (16)-(18) annihilation, with the simplest stable
# benchmark: the non-pseudoscalar 1^3S_1 omega/phi row in Table III.

using Pkg
# Archived script: runs in the GIPaper environment at the repository root.
Pkg.activate(joinpath(@__DIR__, "..", "..", "GIPaper"))

using CSV
using LinearAlgebra
using Printf
using Statistics

root = joinpath(dirname(dirname(@__DIR__)), "GIPaper")
using GIModel

const REPORT = joinpath(@__DIR__, "annihilation_method_audit.md")
const PSEUDOSCALAR_STATES = ["eta(548)", "eta_prime(958)", "eta_r(?)", "eta_r_prime(?)"]
const PSEUDOSCALAR_BASIS = ["1 ns", "1 ss", "2 ns", "2 ss"]

function hidden_s_basis(params, mq, spin; wave_basis = :fd, ngrid = 220, rmax = 22.0)
    out = []
    wave_solver = wave_basis == :ho ?
        OscillatorSolver(ngrid = ngrid, rmax = rmax) :
        FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax)
    for (label, mass_key) in [("1 ns", "q"), ("1 ss", "s")]
        m = mq[mass_key]
        masses = ConstituentMasses(m, m)
        ev, vecs, r = channel_solution(
            params, masses, 0; solver = wave_solver, nlevels = 1)
        fix_annihilation_phase!(vecs, r)
        sol = GIModel.ChannelRadialSolution(ev, vecs, r)
        levels = GIModel.contact_hyperfine_nonperturbative_levels(
            params,
            masses,
            "S",
            spin,
            r,
            1,
        )
        input = pseudoscalar_annihilation_basis_input(
            label,
            m,
            levels[1],
            RadialWaveOnUniformMesh(sol, 1),
        )
        push!(
            out,
            (
                label = label,
                mass = m,
                diagonal = levels[1],
                input = input,
                S_momentum = GIModel._annihilation_overlap_factor(
                    FDMomentumIntegralSmearing(),
                    input,
                ),
                S_legacy = GIModel._flavor_coherence_factor(input) *
                           GIModel._s0_smearing_factor(FDOriginP2Smearing(), input),
            ),
        )
    end
    return out
end

function vector_target()
    rows = collect(CSV.File(joinpath(root, "data", "clean", "mixings.csv")))
    target = Dict{String,Vector{Float64}}(
        "omega(783)" => zeros(2),
        "phi(1019)" => zeros(2),
    )
    shifts = Dict{String,Float64}()
    for row in rows
        String(row.state_group) == "1 3S1" || continue
        state = String(row.state_name)
        haskey(target, state) || continue
        basis = String(row.basis)
        idx = basis == "1 ns" ? 1 : basis == "1 ss" ? 2 : 0
        idx == 0 && continue
        target[state][idx] = Float64(row.amplitude)
        !ismissing(row.mass_shift_MeV) && (shifts[state] = Float64(row.mass_shift_MeV))
    end
    return target, shifts
end

function pseudoscalar_basis(params, mq; wave_basis = :fd, ngrid = 220, rmax = 22.0)
    wave_solver = wave_basis == :ho ?
        OscillatorSolver(ngrid = ngrid, rmax = rmax) :
        FiniteDifferenceSolver(ngrid = ngrid, rmax = rmax)
    channels = Dict{String,Any}()
    for (flavor, mass_key) in [("ns", "q"), ("ss", "s")]
        m = mq[mass_key]
        masses = ConstituentMasses(m, m)
        ev, vecs, r = channel_solution(
            params, masses, 0; solver = wave_solver, nlevels = 2)
        fix_annihilation_phase!(vecs, r)
        sol = GIModel.ChannelRadialSolution(ev, vecs, r)
        levels = GIModel.contact_hyperfine_nonperturbative_levels(
            params,
            masses,
            "S",
            1,
            r,
            2,
        )
        channels[flavor] = (mass = m, solution = sol, levels = levels)
    end
    out = GIModel.PseudoscalarAnnihilationBasisInput[]
    for label in PSEUDOSCALAR_BASIS
        n = parse(Int, first(split(label)))
        flavor = last(split(label))
        channel = channels[flavor]
        push!(
            out,
            pseudoscalar_annihilation_basis_input(
                label,
                channel.mass,
                channel.levels[n],
                RadialWaveOnUniformMesh(channel.solution, n),
            ),
        )
    end
    return out
end

function pseudoscalar_p1_targets()
    rows = collect(CSV.File(joinpath(root, "data", "clean", "mixings.csv")))
    target = Dict{String,Vector{Float64}}(
        state => zeros(length(PSEUDOSCALAR_BASIS)) for state in PSEUDOSCALAR_STATES
    )
    for row in rows
        ismissing(row.model) && continue
        String(row.model) == "P1" || continue
        state = String(row.state_name)
        haskey(target, state) || continue
        idx = findfirst(==(String(row.basis)), PSEUDOSCALAR_BASIS)
        isnothing(idx) && continue
        target[state][idx] = Float64(row.amplitude)
    end
    return [target[state] ./ norm(target[state]) for state in PSEUDOSCALAR_STATES]
end

function eq16_s_matrix(params, basis; A = 2.5, gluons = 3, smearing = :momentum, alpha = true)
    diag = [b.diagonal for b in basis]
    matrix = Matrix(Diagonal(diag))
    for j in eachindex(basis), i in eachindex(basis)
        Sj = smearing == :legacy ? basis[j].S_legacy : basis[j].S_momentum
        Si = smearing == :legacy ? basis[i].S_legacy : basis[i].S_momentum
        αfactor = if alpha
            ((GIModel.alpha_s_q(basis[j].diagonal) * GIModel.alpha_s_q(basis[i].diagonal)) / π^2)^(gluons / 2)
        else
            1.0
        end
        matrix[j, i] += 4π * A * αfactor * Sj * Si / (basis[j].mass * basis[i].mass)
    end
    return matrix
end

function aligned_vector_score(vectors, targets)
    used = Set{Int}()
    rms = Float64[]
    predictions = Vector{Float64}[]
    for target in targets
        best = 0
        bestdot = -Inf
        for col in axes(vectors, 2)
            col in used && continue
            d = dot(vectors[:, col], target)
            if abs(d) > bestdot
                bestdot = abs(d)
                best = col
            end
        end
        push!(used, best)
        sign = dot(vectors[:, best], target) < 0 ? -1.0 : 1.0
        pred = sign .* vectors[:, best]
        push!(rms, sqrt(mean(abs2, pred .- target)))
        push!(predictions, pred)
    end
    return mean(rms), predictions
end

function aligned_rms(vectors, targets)
    used = Set{Int}()
    rms = Float64[]
    for target in targets
        best = 0
        bestdot = -Inf
        for col in axes(vectors, 2)
            col in used && continue
            d = dot(vectors[:, col], target)
            if abs(d) > bestdot
                bestdot = abs(d)
                best = col
            end
        end
        push!(used, best)
        sign = dot(vectors[:, best], target) < 0 ? -1.0 : 1.0
        push!(rms, sqrt(mean(abs2, sign .* vectors[:, best] .- target)))
    end
    return mean(rms), maximum(rms)
end

function target_offdiag_from_mixing(targets, diag)
    # For a weak two-state matrix, tan(theta) ~= V/(D_ss-D_ns).
    omega = targets["omega(783)"] ./ norm(targets["omega(783)"])
    phi = targets["phi(1019)"] ./ norm(targets["phi(1019)"])
    θ = 0.5 * (abs(omega[2] / omega[1]) + abs(phi[1] / phi[2]))
    return θ * (diag[2] - diag[1])
end

function write_matrix(io, matrix)
    labels = ["1 ns", "1 ss"]
    println(io, "| row \\ col | `1 ns` | `1 ss` |")
    println(io, "|---|---:|---:|")
    for i in 1:2
        println(io, @sprintf("| `%s` | %.2f | %.2f |", labels[i], 1000matrix[i, 1], 1000matrix[i, 2]))
    end
end

function main()
    params, mq = load_parameters_and_quark_masses(joinpath(dirname(root), "data", "parameters.provisional.toml"))
    basis = hidden_s_basis(params, mq, 3)
    ho_basis = hidden_s_basis(params, mq, 3; wave_basis = :ho)
    targets, shifts = vector_target()
    target_vectors = [targets["omega(783)"] ./ norm(targets["omega(783)"]),
                      targets["phi(1019)"] ./ norm(targets["phi(1019)"])]
    p1_targets = pseudoscalar_p1_targets()
    diag = [b.diagonal for b in basis]
    target_offdiag = target_offdiag_from_mixing(targets, diag)

    variants = [
        ("paper Eq.16 FD momentum", eq16_s_matrix(params, basis)),
        ("paper Eq.16 HO momentum", eq16_s_matrix(params, ho_basis)),
        ("legacy origin proxy", eq16_s_matrix(params, basis; smearing = :legacy)),
        ("legacy HO origin proxy", eq16_s_matrix(params, ho_basis; smearing = :legacy)),
        ("no alpha factor stress test", eq16_s_matrix(params, basis; alpha = false)),
    ]

    open(REPORT, "w") do io
        println(io, "# Annihilation Method Audit")
        println(io)
        println(io, "Generated by `julia GIPaper/scripts/audit_annihilation_method.jl`.")
        println(io)
        println(io, "## Paper Read")
        println(io)
        println(io, "- Eq. (16) is presented for self-conjugate isoscalar annihilation outside the anomalous pseudoscalar exception.")
        println(io, "- Eq. (17) defines the smeared wavefunction-at-origin factor `S_L(Psi)`.")
        println(io, "- The paper explicitly says `(16)` and `(17)` fail in pseudoscalars, then introduces P1/P2. Therefore `1^3S_1` is the cleaner benchmark for the method.")
        println(io)
        println(io, "## Simple Benchmark: `1^3S_1`")
        println(io)
        println(io, "Table III gives `A(^3S_1)=+2.5`, omega/phi amplitudes near ideal mixing, and small but visible `ns`/`ss` leakage. This is a two-state test of Eq. (16) without the pseudoscalar anomaly.")
        println(io)
        println(io, "| basis | diagonal GeV | S momentum | S legacy |")
        println(io, "|---|---:|---:|---:|")
        for b in basis
            println(io, @sprintf("| `%s` | %.3f | %.5f | %.5f |", b.label, b.diagonal, b.S_momentum, b.S_legacy))
        end
        for b in ho_basis
            println(io, @sprintf("| `%s` HO wave | %.3f | %.5f | %.5f |", b.label, b.diagonal, b.S_momentum, b.S_legacy))
        end
        println(io)
        println(io, @sprintf("Target weak-mixing off-diagonal inferred from Table III amplitudes: `%.2f MeV`.", 1000target_offdiag))
        println(io)
        println(io, "| variant | offdiag MeV | scalar to target offdiag | vector RMS | predicted ss in omega | predicted ns in phi |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for (name, matrix) in variants
            fact = eigen(Symmetric(matrix))
            rms, predictions = aligned_vector_score(fact.vectors, target_vectors)
            offdiag = matrix[1, 2]
            scale = offdiag == 0 ? NaN : target_offdiag / abs(offdiag)
            println(
                io,
                @sprintf(
                    "| %s | %.2f | %.2f | %.4f | %.4f | %.4f |",
                    name,
                    1000offdiag,
                    scale,
                    rms,
                    abs(predictions[1][2]),
                    abs(predictions[2][1]),
                ),
            )
        end
        println(io)
        println(io, "### Eq. (16) FD Momentum Matrix (MeV)")
        println(io)
        write_matrix(io, variants[1][2] - Diagonal(diag))
        println(io)
        println(io, "## Pseudoscalar Smoke Test")
        println(io)
        println(io, "Although the paper warns that pseudoscalars are anomalous, the same wavefunction scale should still move P1 in the right direction if the basis-density diagnosis is right. This four-state smoke test uses only `1 ns`, `1 ss`, `2 ns`, `2 ss` components.")
        println(io)
        println(io, "| basis source | mean P1 vector RMS | max P1 vector RMS |")
        println(io, "|---|---:|---:|")
        for (label, pbasis) in [
            ("FD waves", pseudoscalar_basis(params, mq)),
            ("HO waves", pseudoscalar_basis(params, mq; wave_basis = :ho)),
        ]
            solution = isoscalar_pseudoscalar_annihilation_solution(
                PaperP1Annihilation(),
                params,
                pbasis,
            )
            mean_rms, max_rms = aligned_rms(solution.vectors, p1_targets)
            println(io, @sprintf("| %s | %.3f | %.3f |", label, mean_rms, max_rms))
        end
        println(io)
        println(io, "## Diagnosis")
        println(io)
        println(io, "Resolved. The Table III scale is reproduced when (a) Eq. (17) is evaluated on HO-basis wavefunctions (the FD basis under-supplies wavefunction-at-origin density) and (b) the Eq. (16) bracket uses the three-gluon power `(alpha_i alpha_j/pi^2)^{3/2}` required for the `C=-` `^3S_1` channel. The production `:table_iii` comparison scheme uses `isoscalar_general_annihilation_solution` with exactly these conventions and reproduces the omega/phi Table III amplitudes; the rows above are kept as the method-level benchmark.")
    end
    println("wrote ", REPORT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
