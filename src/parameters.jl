# GI Hamiltonian / smearing / fine-structure switches from the parameters TOML.
#
# Public API (exported from GIModel.jl): GIParameters, load_parameters

struct GIParameters
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
    )
end

function load_parameters(path::AbstractString)
    return gi_parameters_from_raw(TOML.parsefile(path))
end
