# Public API (exported from GIModel.jl):
#   central_potential_values

"""
    central_potential_values(params, m1, m2, r; method = params.central)
    central_potential_values(params, masses::ConstituentMasses, r; kwargs...)

Diagonal central potential `V(r)` (GeV) on the mesh `r` for the given
constituent masses, evaluated through the [`CentralPotentialMethod`](@ref)
`method`. Pass `method` explicitly to compare constructions on identical
inputs (as the Appendix-A audit scripts do). Note [`AppendixAMomentumSandwich`](@ref)
returns the same closed-form diagonal as [`AppendixAClosedForm`](@ref); the
nonlocal sandwich part lives in the Hamiltonian builder, not in these
pointwise values.
"""
function central_potential_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector;
    method::CentralPotentialMethod = params.central,
)
    return central_potential_values(method, params, m1, m2, r)
end

function central_potential_values(
    params::GIParameters,
    masses::ConstituentMasses,
    r::AbstractVector;
    kwargs...,
)
    return central_potential_values(params, masses.m1_GeV, masses.m2_GeV, r; kwargs...)
end

central_potential_values(
    ::PointwiseCentral,
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector,
) = [central_potential(ri, params) for ri in r]

central_potential_values(
    ::AppendixASmearing3D,
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector,
) = smeared_central_values(params, m1, m2, r)

central_potential_values(
    ::Coulomb1DSmearing,
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector,
) = coulomb_1d_smeared_central_values(params, m1, m2, r)

central_potential_values(
    ::AppendixADerivativeG,
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector,
) = appendix_a_derivative_central_values(params, m1, m2, r)

central_potential_values(
    ::Union{AppendixAClosedForm,AppendixAMomentumSandwich},
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector,
) = appendix_a_closed_central_values(params, m1, m2, r)

function potential_diagonal(params::GIParameters, m1::Real, m2::Real, r::AbstractVector)
    return central_potential_values(params, m1, m2, r)
end

function potential_diagonal(params::GIParameters, masses::ConstituentMasses, r::AbstractVector)
    return potential_diagonal(params, masses.m1_GeV, masses.m2_GeV, r)
end
