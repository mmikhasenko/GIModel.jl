#!/usr/bin/env julia
# Numerical readiness check for a reverse-mode state-energy primitive.
# This verifies the eigenvalue pullback on GIModel's actual HO Hamiltonian; it
# deliberately does not differentiate the adaptive beta/basis search.

using Pkg
Pkg.activate(@__DIR__)

using GIModel
using LinearAlgebra
using Printf

const BETA_GEV = 0.65
const NBASIS = 12
const STEP = 1.0e-5

params, mq = load_parameters_and_quark_masses(default_parameters_path())
masses = ConstituentMasses(mq["c"], mq["c"])

function hamiltonian_at_b(b)
    varied = GIParameters(
        params;
        potential = ConfinementPotential(params.potential; b = b),
    )
    return Matrix(GIModel.oscillator_central_matrix(
        varied,
        masses,
        0,
        BETA_GEV,
        NBASIS,
    ))
end

lowest_energy(b) = eigmin(Symmetric(hamiltonian_at_b(b)))

b = params.potential.b
factorization = eigen(Symmetric(hamiltonian_at_b(b)))
energy = first(factorization.values)
wave = view(factorization.vectors, :, 1)

dH_db = (
    hamiltonian_at_b(b + STEP) - hamiltonian_at_b(b - STEP)
) / (2STEP)
dE_db_finite_difference = (
    lowest_energy(b + STEP) - lowest_energy(b - STEP)
) / (2STEP)
dE_db_hellmann_feynman = dot(wave, dH_db * wave)
relative_difference = abs(
    dE_db_finite_difference - dE_db_hellmann_feynman,
) / abs(dE_db_finite_difference)

@printf "HO charmonium 1S: beta = %.2f GeV, nbasis = %d\n" BETA_GEV NBASIS
@printf "E = %.16f GeV\n" energy
@printf "dE/db finite difference = %.16f GeV^-1\n" dE_db_finite_difference
@printf "dE/db Hellmann-Feynman = %.16f GeV^-1\n" dE_db_hellmann_feynman
@printf "relative difference = %.3e\n" relative_difference

relative_difference < 1.0e-8 || error(
    "Hellmann-Feynman pullback failed the finite-difference check",
)
