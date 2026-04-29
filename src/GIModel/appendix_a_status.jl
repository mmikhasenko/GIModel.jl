# Spin-independent central potential: which paper path the active parameters use.
# See `docs/appendix_a_from_paper.md` and `docs/formula_map.md`.

struct CentralPotentialPath
    name::String
    paper_refs::String
    notes::String
end

"""Return a short description of which Appendix A construction is in effect for `params`."""
function central_potential_path(params::GIParameters)::CentralPotentialPath
    mode = central_potential_mode(params)
    if mode == :appendix_a_derivative_g
        return CentralPotentialPath(
            "appendix_a_derivative_g_proxy",
            "First finite-difference derivative term for Gaussian-smearing G(r): G + ∇²G/(4σ²); S(r)=br+c pointwise; A12/A13 coefficients still PDF-audit gated",
            "Comparator path for Appendix-A work. Precedence: `appendix_a_derivative_g` wins over the older 3D and 1D diagnostic toggles.",
        )
    end
    if mode == :appendix_a_3d_a7a8
        return CentralPotentialPath(
            "experimental_3d_convl_a7a8",
            "(A7)–(A8) smearing style via `smear_3d_radial` on pointwise G and S; not (A12)–(A13)",
            "Can remove small-r binding; default is `appendix_a_smearing = false` in `parameters.provisional.toml`.",
        )
    end
    if mode == :coulomb_1d
        return CentralPotentialPath(
            "coulomb_1d_gauss_on_mesh",
            "1D Gaussian renormalization of G(r) on the radial grid; S(r) = br + c kept pointwise; same σ as contact (A9); not (A12)–(A13)",
            "Precedence: `appendix_a_derivative_g` then `appendix_a_smearing`; otherwise optional `coulomb_1d_smear` in `[potential]`.",
        )
    end
    return CentralPotentialPath(
        "pointwise_fd",
        "Eqs. (11)–(13) orientation: V = b r - 4α_s/(3r) + c on the FD mesh",
        "Semirelativistic kinetic + this V is the main diagnostic baseline until (A12) is implemented.",
    )
end
