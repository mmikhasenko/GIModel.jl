# Public API (exported from GIModel.jl):
#   central_potential_mode, central_potential_values

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

function potential_diagonal(params::GIParameters, m1::Real, m2::Real, r::AbstractVector)
    return central_potential_values(params, m1, m2, r)
end
