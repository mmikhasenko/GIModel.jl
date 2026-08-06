# GI Hamiltonian / smearing / fine-structure parameters from the parameters TOML,
# grouped by model aspect (one substruct per TOML section).
#
# Public API (exported from GIModel.jl):
#   GIParameters, ConfinementPotential, RelativisticSmearing, RelativisticFactors,
#   FineStructure, AnnihilationAmplitudes, CentralPotentialMethod and its
#   singletons, load_parameters

# `GIParameters` used to carry a `Basis` type parameter (`FiniteDifferenceBasis` /
# `HarmonicOscillatorBasis`) that no field ever used: it existed only to dispatch
# the radial solve. That put a statement about *how you intend to discretize* in
# the struct that holds the *model*, and it meant a physics object had to be
# rebuilt to change a numerical method. The choice now lives where it belongs, on
# the solver — see [`RadialSolver`](@ref) and its two implementations.

"""
    ConfinementPotential(; b, c)

Linear-plus-constant confinement from Table II: slope `b` in GeV² and offset
`c` in GeV (the TOML stores `c_MeV`). The Coulomb part comes from the running
coupling, not from these constants.
"""
Base.@kwdef struct ConfinementPotential
    b::Float64
    c::Float64
end

"""
    RelativisticSmearing(; sigma0, s)

Appendix A (A9) universal smearing width inputs: `sigma0` in GeV and the
dimensionless `s`, combined per quark-mass pair by `contact_smearing_sigma`.
"""
Base.@kwdef struct RelativisticSmearing
    sigma0::Float64
    s::Float64
end

"""
    CentralPotentialMethod

Which construction evaluates the spin-independent central potential. One of
[`PointwiseCentral`](@ref), [`Coulomb1DSmearing`](@ref),
[`AppendixASmearing3D`](@ref), [`AppendixADerivativeG`](@ref),
[`AppendixAClosedForm`](@ref), [`AppendixAMomentumSandwich`](@ref).
Selected by the `central` key of the `[potential]` TOML section; the methods
are mutually exclusive by construction (no precedence rules).
"""
abstract type CentralPotentialMethod end

"""
    PointwiseCentral

Raw diagnostic baseline: `V = b r - 4α_s/(3r) + c` evaluated pointwise on the
FD mesh (Eqs. (11)–(13) orientation).
"""
struct PointwiseCentral <: CentralPotentialMethod end

"""
    Coulomb1DSmearing

Research comparator: 1D Gaussian renormalization of `G(r)` on the radial grid
only; `S(r) = br + c` kept pointwise. Same σ as contact (A9); not (A12)–(A13).
"""
struct Coulomb1DSmearing <: CentralPotentialMethod end

"""
    AppendixASmearing3D

Experimental (A7)–(A8)-style 3D Gaussian smearing of pointwise `G` and `S`
with the contact-hyperfine σ (A9). Can wipe the small-r Coulomb well on a
radial grid; kept as a comparator, not the active reproduction path.
"""
struct AppendixASmearing3D <: CentralPotentialMethod end

"""
    AppendixADerivativeG

Research comparator: first derivative-expansion term for the Gaussian-smeared
Coulomb `G(r)` (`G + ∇²G/(4σ²)`), with `S(r)` pointwise.
"""
struct AppendixADerivativeG <: CentralPotentialMethod end

"""
    AppendixAClosedForm

Closed-form Gaussian smearing of GI `G(r)` and `S(r)` via the checked
(A12)–(A14) τ_k and smeared-linear formulas. Diagonal only: omits the central
Coulomb momentum sandwich `G' = A(p) G̃ A(p)`.
"""
struct AppendixAClosedForm <: CentralPotentialMethod end

"""
    AppendixAMomentumSandwich

Active reproduction path: closed-form `G̃`, `S̃` plus the nonlocal Coulomb
momentum sandwich `A(p) G̃ A(p)`, `A(p) = sqrt(1 + p²/(E₁E₂))`, on the FD basis.
"""
struct AppendixAMomentumSandwich <: CentralPotentialMethod end

const CENTRAL_POTENTIAL_METHODS = Dict{String,CentralPotentialMethod}(
    "pointwise" => PointwiseCentral(),
    "coulomb_1d_smear" => Coulomb1DSmearing(),
    "appendix_a_smearing" => AppendixASmearing3D(),
    "appendix_a_derivative_g" => AppendixADerivativeG(),
    "appendix_a_closed_form" => AppendixAClosedForm(),
    "appendix_a_momentum_sandwich" => AppendixAMomentumSandwich(),
)

"""
    central_potential_method(name::AbstractString) -> CentralPotentialMethod

Resolve the `central` key of the `[potential]` TOML section to its
[`CentralPotentialMethod`](@ref) singleton; throws `ArgumentError` for unknown names.
"""
function central_potential_method(name::AbstractString)
    haskey(CENTRAL_POTENTIAL_METHODS, String(name)) || throw(ArgumentError(
        "unknown central potential method `$name`; expected one of: " *
        join(sort(collect(keys(CENTRAL_POTENTIAL_METHODS))), ", "),
    ))
    return CENTRAL_POTENTIAL_METHODS[String(name)]
end

"""
    RelativisticFactors(; epsilon_c=0, epsilon_t=0, epsilon_so_vector=0, epsilon_so_scalar=0,
                         contact_momentum_sandwich=false,
                         fine_structure_momentum_sandwich=false,
                         fine_structure_smeared_kernels=false)

Post-(A14) relativistic ε factors and the switches choosing how they are
applied: `contact_momentum_sandwich` wraps the smeared contact operator in the
GI Hermitian momentum-factor sandwich; `fine_structure_momentum_sandwich` does
the same for first-order spin-orbit/tensor expectations;
`fine_structure_smeared_kernels` uses derivatives of the closed-form smeared
`G̃`/`S̃` for those radial kernels instead of the pointwise running-coupling forms.
"""
Base.@kwdef struct RelativisticFactors
    epsilon_c::Float64 = 0.0
    epsilon_t::Float64 = 0.0
    epsilon_so_vector::Float64 = 0.0
    epsilon_so_scalar::Float64 = 0.0
    contact_momentum_sandwich::Bool = false
    fine_structure_momentum_sandwich::Bool = false
    fine_structure_smeared_kernels::Bool = false
end

"""
    FineStructure(; enabled=true, k_spin_orbit=0.5, k_tensor=0.4)

Vector + Thomas spin-orbit and OGE-tensor first-order corrections on the FD
radial `u(r)`: master switch and the global scales aligning with the paper's
HO result.
"""
Base.@kwdef struct FineStructure
    enabled::Bool = true
    k_spin_orbit::Float64 = 0.5
    k_tensor::Float64 = 0.4
end

"""
    AnnihilationAmplitudes(; p1_A_np=0.5, p1_m_eta=0.548, p2_A_np=0.55, p2_M0=1.17,
                            s1_A=2.5, a_3p2=-0.8)

Table III pseudoscalar/vector/tensor annihilation constants: Eq. (18a) `P1`
(`p1_A_np`, `p1_m_eta` in GeV), Eq. (18b) `P2` (`p2_A_np`, zero at `p2_M0` GeV),
and the Eq. (16) channel amplitudes `A(^3S_1) = s1_A`, `A(^3P_2) = a_3p2`.
"""
Base.@kwdef struct AnnihilationAmplitudes
    p1_A_np::Float64 = 0.5
    p1_m_eta::Float64 = 0.548
    p2_A_np::Float64 = 0.55
    p2_M0::Float64 = 1.17
    s1_A::Float64 = 2.5
    a_3p2::Float64 = -0.8
end

"""
    GIParameters

Model parameters grouped by aspect, one field per TOML section:

  - `potential::ConfinementPotential` — `b`, `c` constants.
  - `central::CentralPotentialMethod` — which central-potential construction runs.
  - `smearing::RelativisticSmearing` — (A9) σ₀ and s.
  - `factors::RelativisticFactors` — ε factors and sandwich/kernel switches.
  - `fine_structure::FineStructure` — spin-orbit/tensor master switch and scales.
  - `annihilation::AnnihilationAmplitudes` — Table III constants.

Construct via [`load_parameters`](@ref) (TOML) or keywords:
`GIParameters(potential = ConfinementPotential(b = 0.18, c = -0.253), ...)`.

This is the **model**, and nothing more: every field is a number the paper
quotes. How the resulting Schrödinger equation gets solved is a separate choice,
carried by the [`RadialSolver`](@ref) you pass to
[`channel_solution`](@ref) or [`compute_spectrum`](@ref).
"""
Base.@kwdef struct GIParameters
    potential::ConfinementPotential
    central::CentralPotentialMethod = PointwiseCentral()
    smearing::RelativisticSmearing
    factors::RelativisticFactors = RelativisticFactors()
    fine_structure::FineStructure = FineStructure()
    annihilation::AnnihilationAmplitudes = AnnihilationAmplitudes()
end

function gi_parameters_from_raw(raw)::GIParameters
    pot = raw["potential"]
    rf = get(raw, "relativistic_factors", Dict{String,Any}())
    fs = get(raw, "fine_structure", Dict{String,Any}())
    ann = get(raw, "annihilation", Dict{String,Any}())
    return GIParameters(
        potential = ConfinementPotential(
            b = pot["b_GeV2"],
            c = pot["c_MeV"] / 1000,
        ),
        central = central_potential_method(get(pot, "central", "pointwise")),
        smearing = RelativisticSmearing(
            sigma0 = raw["relativistic_smearing"]["sigma0_GeV"],
            s = raw["relativistic_smearing"]["s"],
        ),
        factors = RelativisticFactors(
            epsilon_c = float(get(rf, "epsilon_c", 0.0)),
            epsilon_t = float(get(rf, "epsilon_t", 0.0)),
            epsilon_so_vector = float(get(rf, "epsilon_so_vector", 0.0)),
            epsilon_so_scalar = float(get(rf, "epsilon_so_scalar", 0.0)),
            contact_momentum_sandwich = get(rf, "contact_momentum_sandwich", false),
            fine_structure_momentum_sandwich = get(
                rf,
                "fine_structure_momentum_sandwich",
                false,
            ),
            fine_structure_smeared_kernels = get(
                rf,
                "fine_structure_smeared_kernels",
                false,
            ),
        ),
        fine_structure = FineStructure(
            enabled = get(fs, "enabled", true),
            k_spin_orbit = float(get(fs, "k_spin_orbit", 0.5)),
            k_tensor = float(get(fs, "k_tensor", 0.4)),
        ),
        annihilation = AnnihilationAmplitudes(
            p1_A_np = float(get(ann, "p1_A_np", 0.5)),
            p1_m_eta = float(get(ann, "p1_m_eta_GeV", 0.548)),
            p2_A_np = float(get(ann, "p2_A_np", 0.55)),
            p2_M0 = float(get(ann, "p2_M0_GeV", 1.17)),
            s1_A = float(get(ann, "s1_A", 2.5)),
            a_3p2 = float(get(ann, "a_3p2", -0.8)),
        ),
    )
end

"""
    load_parameters(path) -> GIParameters

Read a solver-parameter TOML file (`data/parameters.provisional.toml` layout:
`[potential]` with a `central` method name, `[relativistic_smearing]`,
`[relativistic_factors]`, `[fine_structure]`, optional `[annihilation]`).
Missing switches default to `false`/paper values. Use
[`load_parameters_and_quark_masses`](@ref) to also get the `[masses]` table.
The file describes the model only; pick the radial method with a
[`RadialSolver`](@ref) at the call site.
"""
function load_parameters(path::AbstractString)
    return gi_parameters_from_raw(TOML.parsefile(path))
end
