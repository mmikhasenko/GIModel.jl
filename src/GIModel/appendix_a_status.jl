# Spin-independent central potential: which paper path the active parameters use.
# See `docs/appendix_a_from_paper.md` and `docs/formula_map.md`.

struct CentralPotentialPath
    name::String
    paper_refs::String
    notes::String
end

"""Return a short description of which Appendix A construction is in effect for `params`."""
function central_potential_path(params::GIParameters)::CentralPotentialPath
    if params.appendix_a_smearing
        return CentralPotentialPath(
            "experimental_3d_convl_a7a8",
            "(A7)–(A8) smearing style via `smear_3d_radial` on pointwise G and S; not (A12)–(A13)",
            "Can remove small-r binding; default is `appendix_a_smearing = false` in `parameters.provisional.toml`.",
        )
    end
    return CentralPotentialPath(
        "pointwise_fd",
        "Eqs. (11)–(13) orientation: V = b r - 4α_s/(3r) + c on the FD mesh",
        "Semirelativistic kinetic + this V is the main diagnostic baseline until (A12) is implemented.",
    )
end
