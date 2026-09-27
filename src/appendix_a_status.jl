# Spin-independent central potential: which paper path the active parameters use.
# See `docs/appendix_a_from_paper.md` and `docs/formula_map.md`.
#
# Public API (exported from GIModel.jl):
#   CentralPotentialPath, central_potential_path

"""
    CentralPotentialPath(name, paper_refs, notes)

Human-readable description of the spin-independent central construction in
effect, as returned by [`central_potential_path`](@ref): the method `name`, the
paper equations it implements (`paper_refs`) and short `notes`.
"""
struct CentralPotentialPath
    name::String
    paper_refs::String
    notes::String
end

"""Return a short description of which Appendix A construction is in effect for `params`."""
central_potential_path(params::GIParameters)::CentralPotentialPath =
    central_potential_path(params.central)

central_potential_path(::AppendixAMomentumSandwich) = CentralPotentialPath(
    "appendix_a_momentum_sandwich",
    "Closed-form smeared G̃(r), S̃(r), plus central Coulomb momentum sandwich G' = A(p)G̃A(p) in the solver's native p² representation",
    "Shipped spin-independent GI central prescription. The central-method type selects exactly one construction.",
)

central_potential_path(::AppendixAClosedForm) = CentralPotentialPath(
    "appendix_a_closed_form",
    "Closed-form Gaussian-smeared G̃(r) and S̃(r) using τ_k and smeared-linear formulas; no G' momentum sandwich",
    "Diagonal bracket path. Use to isolate smearing effects before enabling `appendix_a_momentum_sandwich`.",
)

central_potential_path(::AppendixADerivativeG) = CentralPotentialPath(
    "appendix_a_derivative_g",
    "First finite-difference derivative term for Gaussian-smearing G(r): G + ∇²G/(4σ²); S(r)=br+c pointwise",
    "Derivative-expansion diagnostic, selected explicitly; not the shipped central prescription.",
)

central_potential_path(::AppendixASmearing3D) = CentralPotentialPath(
    "experimental_3d_convl_a7a8",
    "(A7)–(A8) smearing style via `smear_3d_radial` on pointwise G and S; not (A12)–(A13)",
    "Can remove small-r binding; the shipped preset uses the closed-form momentum sandwich.",
)

central_potential_path(::Coulomb1DSmearing) = CentralPotentialPath(
    "coulomb_1d_gauss_on_mesh",
    "1D Gaussian renormalization of G(r) on the radial grid; S(r) = br + c kept pointwise; same σ as contact (A9); not (A12)–(A13)",
    "Selected by `central = \"coulomb_1d_smear\"` in `[potential]`.",
)

central_potential_path(::PointwiseCentral) = CentralPotentialPath(
    "pointwise_fd",
    "Eqs. (11)–(13) orientation: V = b r - 4α_s/(3r) + c on the FD mesh",
    "Semirelativistic kinetic + this V is the raw diagnostic baseline; the shipped prescription is the closed-form Appendix-A momentum-sandwich path.",
)
