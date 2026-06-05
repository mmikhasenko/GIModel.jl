# GI Hamiltonian / smearing / fine-structure switches from the parameters TOML.
#
# Public API (exported from GIModel.jl): GIParameters, load_parameters

abstract type GIBasis end

"""
    FiniteDifferenceBasis

Current radial-coordinate finite-difference basis. This is the default basis
for `GIParameters(...)`; specialized GI/HO-basis paths can be added by defining
methods on `GIParameters{<:GIBasis}` without adding more runtime switches.
"""
struct FiniteDifferenceBasis <: GIBasis end

"""
    HarmonicOscillatorBasis

Paper-style oscillator expansion path. The implementation uses oscillator
radial basis functions on the same diagnostic mesh, projects `p^2` and radial
operators into a finite oscillator space, scans the oscillator scale `beta`,
and reconstructs eigenvectors onto the mesh for common reporting.
"""
struct HarmonicOscillatorBasis <: GIBasis end

struct GIParameters{Basis<:GIBasis}
    b::Float64
    c::Float64
    sigma0::Float64
    smearing_s::Float64
    appendix_a_smearing::Bool
    appendix_a_derivative_g::Bool
    appendix_a_closed_form::Bool
    appendix_a_momentum_sandwich::Bool
    contact_momentum_sandwich::Bool
    epsilon_c::Float64
    epsilon_t::Float64
    epsilon_so_vector::Float64
    epsilon_so_scalar::Float64
    fine_structure_momentum_sandwich::Bool
    fine_structure_smeared_kernels::Bool
    fine_structure::Bool
    k_spin_orbit::Float64
    k_tensor::Float64
    coulomb_1d_smear::Bool
    annihilation_p1_A_np::Float64
    annihilation_p1_m_eta::Float64
    annihilation_p2_A_np::Float64
    annihilation_p2_M0::Float64
    annihilation_s1_A::Float64
end

function GIParameters(args...)
    return GIParameters{FiniteDifferenceBasis}(args...)
end

basis_type(::GIParameters{Basis}) where {Basis<:GIBasis} = Basis

function with_basis(params::GIParameters, ::Type{Basis}) where {Basis<:GIBasis}
    return GIParameters{Basis}(
        params.b,
        params.c,
        params.sigma0,
        params.smearing_s,
        params.appendix_a_smearing,
        params.appendix_a_derivative_g,
        params.appendix_a_closed_form,
        params.appendix_a_momentum_sandwich,
        params.contact_momentum_sandwich,
        params.epsilon_c,
        params.epsilon_t,
        params.epsilon_so_vector,
        params.epsilon_so_scalar,
        params.fine_structure_momentum_sandwich,
        params.fine_structure_smeared_kernels,
        params.fine_structure,
        params.k_spin_orbit,
        params.k_tensor,
        params.coulomb_1d_smear,
        params.annihilation_p1_A_np,
        params.annihilation_p1_m_eta,
        params.annihilation_p2_A_np,
        params.annihilation_p2_M0,
        params.annihilation_s1_A,
    )
end

function gi_parameters_from_raw(raw)::GIParameters
    rf = get(raw, "relativistic_factors", nothing)
    eps_c = isnothing(rf) ? 0.0 : get(rf, "epsilon_c", 0.0)
    eps_t = isnothing(rf) ? 0.0 : get(rf, "epsilon_t", 0.0)
    eps_v = isnothing(rf) ? 0.0 : get(rf, "epsilon_so_vector", 0.0)
    eps_s = isnothing(rf) ? 0.0 : get(rf, "epsilon_so_scalar", 0.0)
    fs = get(raw, "fine_structure", nothing)
    fine_on = isnothing(fs) ? true : get(fs, "enabled", true)
    k_so = isnothing(fs) ? 0.5 : get(fs, "k_spin_orbit", 0.5)
    k_tn = isnothing(fs) ? 0.4 : get(fs, "k_tensor", 0.4)
    ann = get(raw, "annihilation", nothing)
    p1_A = isnothing(ann) ? 0.5 : get(ann, "p1_A_np", 0.5)
    p1_meta = isnothing(ann) ? 0.548 : get(ann, "p1_m_eta_GeV", 0.548)
    p2_A = isnothing(ann) ? 0.55 : get(ann, "p2_A_np", 0.55)
    p2_M0 = isnothing(ann) ? 1.17 : get(ann, "p2_M0_GeV", 1.17)
    s1_A = isnothing(ann) ? 2.5 : get(ann, "s1_A", 2.5)
    return GIParameters(
        raw["potential"]["b_GeV2"],
        raw["potential"]["c_MeV"] / 1000,
        raw["relativistic_smearing"]["sigma0_GeV"],
        raw["relativistic_smearing"]["s"],
        get(raw["potential"], "appendix_a_smearing", false),
        get(raw["potential"], "appendix_a_derivative_g", false),
        get(raw["potential"], "appendix_a_closed_form", false),
        get(raw["potential"], "appendix_a_momentum_sandwich", false),
        get(raw["relativistic_factors"], "contact_momentum_sandwich", false),
        float(eps_c),
        float(eps_t),
        float(eps_v),
        float(eps_s),
        get(raw["relativistic_factors"], "fine_structure_momentum_sandwich", false),
        get(raw["relativistic_factors"], "fine_structure_smeared_kernels", false),
        fine_on,
        float(k_so),
        float(k_tn),
        get(raw["potential"], "coulomb_1d_smear", false),
        float(p1_A),
        float(p1_meta),
        float(p2_A),
        float(p2_M0),
        float(s1_A),
    )
end

function load_parameters(path::AbstractString)
    return gi_parameters_from_raw(TOML.parsefile(path))
end
