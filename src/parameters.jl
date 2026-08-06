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
struct ConfinementPotential{T<:Real}
    b::T
    c::T
end

function ConfinementPotential(; b::Real, c::Real)
    promoted = promote(b, c)
    return ConfinementPotential(promoted...)
end

ConfinementPotential(
    base::ConfinementPotential;
    b::Real = base.b,
    c::Real = base.c,
) = ConfinementPotential(; b = b, c = c)

"""
    RelativisticSmearing(; sigma0, s)

Appendix A (A9) universal smearing width inputs: `sigma0` in GeV and the
dimensionless `s`, combined per quark-mass pair by `contact_smearing_sigma`.
"""
struct RelativisticSmearing{T<:Real}
    sigma0::T
    s::T
end

function RelativisticSmearing(; sigma0::Real, s::Real)
    promoted = promote(sigma0, s)
    return RelativisticSmearing(promoted...)
end

RelativisticSmearing(
    base::RelativisticSmearing;
    sigma0::Real = base.sigma0,
    s::Real = base.s,
) = RelativisticSmearing(; sigma0 = sigma0, s = s)

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
struct RelativisticFactors{T<:Real}
    epsilon_c::T
    epsilon_t::T
    epsilon_so_vector::T
    epsilon_so_scalar::T
    contact_momentum_sandwich::Bool
    fine_structure_momentum_sandwich::Bool
    fine_structure_smeared_kernels::Bool
end

function RelativisticFactors(;
    epsilon_c::Real = 0.0,
    epsilon_t::Real = 0.0,
    epsilon_so_vector::Real = 0.0,
    epsilon_so_scalar::Real = 0.0,
    contact_momentum_sandwich::Bool = false,
    fine_structure_momentum_sandwich::Bool = false,
    fine_structure_smeared_kernels::Bool = false,
)
    ec, et, esov, esos = promote(
        epsilon_c, epsilon_t, epsilon_so_vector, epsilon_so_scalar,
    )
    return RelativisticFactors(
        ec,
        et,
        esov,
        esos,
        contact_momentum_sandwich,
        fine_structure_momentum_sandwich,
        fine_structure_smeared_kernels,
    )
end

RelativisticFactors(
    base::RelativisticFactors;
    epsilon_c::Real = base.epsilon_c,
    epsilon_t::Real = base.epsilon_t,
    epsilon_so_vector::Real = base.epsilon_so_vector,
    epsilon_so_scalar::Real = base.epsilon_so_scalar,
    contact_momentum_sandwich::Bool = base.contact_momentum_sandwich,
    fine_structure_momentum_sandwich::Bool = base.fine_structure_momentum_sandwich,
    fine_structure_smeared_kernels::Bool = base.fine_structure_smeared_kernels,
) = RelativisticFactors(;
    epsilon_c = epsilon_c,
    epsilon_t = epsilon_t,
    epsilon_so_vector = epsilon_so_vector,
    epsilon_so_scalar = epsilon_so_scalar,
    contact_momentum_sandwich = contact_momentum_sandwich,
    fine_structure_momentum_sandwich = fine_structure_momentum_sandwich,
    fine_structure_smeared_kernels = fine_structure_smeared_kernels,
)

"""
    FineStructure(; enabled=true, k_spin_orbit=0.5, k_tensor=0.4)

Vector + Thomas spin-orbit and OGE-tensor first-order corrections on the FD
radial `u(r)`: master switch and the global scales aligning with the paper's
HO result.
"""
struct FineStructure{T<:Real}
    enabled::Bool
    k_spin_orbit::T
    k_tensor::T
end

function FineStructure(;
    enabled::Bool = true,
    k_spin_orbit::Real = 0.5,
    k_tensor::Real = 0.4,
)
    spin_orbit, tensor = promote(k_spin_orbit, k_tensor)
    return FineStructure(enabled, spin_orbit, tensor)
end

FineStructure(
    base::FineStructure;
    enabled::Bool = base.enabled,
    k_spin_orbit::Real = base.k_spin_orbit,
    k_tensor::Real = base.k_tensor,
) = FineStructure(;
    enabled = enabled,
    k_spin_orbit = k_spin_orbit,
    k_tensor = k_tensor,
)

"""
    AnnihilationAmplitudes(; p1_A_np=0.5, p1_m_eta=0.548, p2_A_np=0.55, p2_M0=1.17,
                            s1_A=2.5, a_3p2=-0.8)

Table III pseudoscalar/vector/tensor annihilation constants: Eq. (18a) `P1`
(`p1_A_np`, `p1_m_eta` in GeV), Eq. (18b) `P2` (`p2_A_np`, zero at `p2_M0` GeV),
and the Eq. (16) channel amplitudes `A(^3S_1) = s1_A`, `A(^3P_2) = a_3p2`.
"""
struct AnnihilationAmplitudes{T<:Real}
    p1_A_np::T
    p1_m_eta::T
    p2_A_np::T
    p2_M0::T
    s1_A::T
    a_3p2::T
end

function AnnihilationAmplitudes(;
    p1_A_np::Real = 0.5,
    p1_m_eta::Real = 0.548,
    p2_A_np::Real = 0.55,
    p2_M0::Real = 1.17,
    s1_A::Real = 2.5,
    a_3p2::Real = -0.8,
)
    values = promote(p1_A_np, p1_m_eta, p2_A_np, p2_M0, s1_A, a_3p2)
    return AnnihilationAmplitudes(values...)
end

AnnihilationAmplitudes(
    base::AnnihilationAmplitudes;
    p1_A_np::Real = base.p1_A_np,
    p1_m_eta::Real = base.p1_m_eta,
    p2_A_np::Real = base.p2_A_np,
    p2_M0::Real = base.p2_M0,
    s1_A::Real = base.s1_A,
    a_3p2::Real = base.a_3p2,
) = AnnihilationAmplitudes(;
    p1_A_np = p1_A_np,
    p1_m_eta = p1_m_eta,
    p2_A_np = p2_A_np,
    p2_M0 = p2_M0,
    s1_A = s1_A,
    a_3p2 = a_3p2,
)

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
struct GIParameters{
    P<:ConfinementPotential,
    C<:CentralPotentialMethod,
    S<:RelativisticSmearing,
    R<:RelativisticFactors,
    F<:FineStructure,
    A<:AnnihilationAmplitudes,
}
    potential::P
    central::C
    smearing::S
    factors::R
    fine_structure::F
    annihilation::A
end

function GIParameters(;
    potential::ConfinementPotential,
    central::CentralPotentialMethod = PointwiseCentral(),
    smearing::RelativisticSmearing,
    factors::RelativisticFactors = RelativisticFactors(),
    fine_structure::FineStructure = FineStructure(),
    annihilation::AnnihilationAmplitudes = AnnihilationAmplitudes(),
)
    return GIParameters(
        potential,
        central,
        smearing,
        factors,
        fine_structure,
        annihilation,
    )
end

GIParameters(
    base::GIParameters;
    potential::ConfinementPotential = base.potential,
    central::CentralPotentialMethod = base.central,
    smearing::RelativisticSmearing = base.smearing,
    factors::RelativisticFactors = base.factors,
    fine_structure::FineStructure = base.fine_structure,
    annihilation::AnnihilationAmplitudes = base.annihilation,
) = GIParameters(;
    potential = potential,
    central = central,
    smearing = smearing,
    factors = factors,
    fine_structure = fine_structure,
    annihilation = annihilation,
)

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
