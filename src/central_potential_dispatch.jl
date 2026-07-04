# Public API (exported from GIModel.jl):
#   central_potential_mode, central_potential_values

"""
    central_potential_mode(params) -> Symbol

Which central-potential evaluation path the parameter switches select, in
precedence order: `:appendix_a_momentum_sandwich`, `:appendix_a_closed_form`,
`:appendix_a_derivative_g`, `:appendix_a_3d_a7a8`, `:coulomb_1d`, else
`:pointwise`. See [`CentralPotentialPath`](@ref) for the audit-facing summary.
"""
function central_potential_mode(params::GIParameters)::Symbol
    if params.appendix_a_momentum_sandwich
        return :appendix_a_momentum_sandwich
    elseif params.appendix_a_closed_form
        return :appendix_a_closed_form
    elseif params.appendix_a_derivative_g
        return :appendix_a_derivative_g
    elseif params.appendix_a_smearing
        return :appendix_a_3d_a7a8
    elseif params.coulomb_1d_smear
        return :coulomb_1d
    end
    return :pointwise
end

"""
    central_potential_values(params, m1, m2, r; mode = central_potential_mode(params))
    central_potential_values(params, masses::ConstituentMasses, r; kwargs...)

Diagonal central potential `V(r)` (GeV) on the mesh `r` for the given
constituent masses, evaluated through the path chosen by `mode`. Pass `mode`
explicitly to compare paths on identical inputs (as the Appendix-A audit
scripts do). Note the `:appendix_a_momentum_sandwich` mode returns the same
closed-form diagonal as `:appendix_a_closed_form`; the nonlocal sandwich part
lives in the Hamiltonian builder, not in these pointwise values.
"""
function central_potential_values(
    params::GIParameters,
    m1::Real,
    m2::Real,
    r::AbstractVector;
    mode::Symbol = central_potential_mode(params),
)
    if mode == :pointwise
        return [central_potential(ri, params) for ri in r]
    elseif mode == :appendix_a_3d_a7a8
        return smeared_central_values(params, m1, m2, r)
    elseif mode == :coulomb_1d
        return coulomb_1d_smeared_central_values(params, m1, m2, r)
    elseif mode == :appendix_a_derivative_g
        return appendix_a_derivative_central_values(params, m1, m2, r)
    elseif mode == :appendix_a_closed_form
        return appendix_a_closed_central_values(params, m1, m2, r)
    elseif mode == :appendix_a_momentum_sandwich
        return appendix_a_closed_central_values(params, m1, m2, r)
    else
        error("unknown central potential mode: $mode")
    end
end

function central_potential_values(
    params::GIParameters,
    masses::ConstituentMasses,
    r::AbstractVector;
    kwargs...,
)
    return central_potential_values(params, masses.m1_GeV, masses.m2_GeV, r; kwargs...)
end

function potential_diagonal(params::GIParameters, m1::Real, m2::Real, r::AbstractVector)
    return central_potential_values(params, m1, m2, r)
end

function potential_diagonal(params::GIParameters, masses::ConstituentMasses, r::AbstractVector)
    return potential_diagonal(params, masses.m1_GeV, masses.m2_GeV, r)
end
