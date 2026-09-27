# Equations to code

This page maps each equation of Godfrey and Isgur (1985) to the code that
implements it. Names in `code font` without a link are internal functions;
they are listed so that you can find the implementation, not as supported API.
Paths are relative to the repository root.

## Spectrum: Section II and Appendix A

| paper | content | implementation |
|---|---|---|
| Eq. (1) | ``H = \sqrt{p^2+m_1^2} + \sqrt{p^2+m_2^2} + \tilde V`` | kinetic term in `relativistic_hamiltonian` (`src/hamiltonian.jl`) for the mesh; [`ho_p2_matrix`](@ref) and the matrix square root in `src/harmonic_oscillator_basis.jl` for the oscillator basis |
| Eqs. (2)–(3) | color-singlet ``G(r) = -4\alpha_s(r)/3r``, ``S(r) = br + c`` | `static_coulomb_G`, `static_confinement_S` (`src/running_coupling.jl`); parameters in [`ConfinementPotential`](@ref) |
| Eq. (4) | contact and tensor hyperfine terms | through Eq. (A15), below |
| Eqs. (6)–(7) | vector and Thomas spin-orbit terms | through Eqs. (A15)–(A16), below |
| Eqs. (12)–(13), Fig. 2 | running coupling ``\alpha_s(Q^2)`` and ``\alpha_s(r) = \sum_k \alpha_k\,\mathrm{erf}(\gamma_k r)`` | [`alpha_s_q`](@ref), [`alpha_s_r`](@ref); constants `ALPHA_COEFFS`, `ALPHA_GAMMAS` (`src/constants.jl`) |
| Eq. (14) and text | staged diagonalization | [`fixed_spectrum`](@ref) (stage 1), [`add_intra_meson_mixing`](@ref) (stage 2), [`add_isoscalar_annihilation`](@ref) (stage 3) |
| Eqs. (A7)–(A8) | Gaussian smearing of a potential | direct convolution `smear_3d_radial` (`src/smearing_appendix_a.jl`), used only by the diagnostic [`AppendixASmearing3D`](@ref) method |
| Eq. (A9) | smearing width ``\sigma_{12}(m_1, m_2)`` | [`contact_smearing_sigma`](@ref) |
| Eqs. (A12)–(A14) | closed-form smeared ``\tilde G``, ``\tilde S``, ``\tau_k`` | `smeared_coulomb_G_closed`, `smeared_confinement_S_closed` (`src/smearing_appendix_a.jl`) |
| after Eq. (A14) | momentum factors ``A(p)\tilde G A(p)`` and ``(m_1m_2/E_1E_2)^{1/2+\epsilon_i}`` | `appendix_a_momentum_sandwich_matrix` (`src/hamiltonian.jl`) for the central term; [`RelativisticFactors`](@ref) holds the ``\epsilon_i`` and switches |
| Eq. (A15) | contact, vector spin-orbit and tensor operators | `smeared_contact_kernel` (`src/contact_hyperfine.jl`); [`fine_structure_radial_kernels`](@ref) and the shared angular assembly in `src/spin_fine_structure.jl` |
| Eq. (A16) | scalar (Thomas) spin-orbit operator | same as (A15), `scalar_11` and `scalar_22` kernels |
| Eq. (A17) | ``\langle i\vert f(p)g(r)\vert j\rangle = \sum_n \langle i\vert f(p)\vert n\rangle\langle n\vert g(r)\vert j\rangle`` | [`ho_operator_matrix`](@ref), [`ho_p2_matrix`](@ref), [`ho_r2_matrix`](@ref); the variational ``\beta`` and basis refinement in [`OscillatorSolver`](@ref) |

Every spin-dependent operator is built once and used both to assemble the
Hamiltonian and to report its contribution, so the reported pieces always add
up to the mass. The finite-difference and oscillator solvers share the angular
and mass algebra; they differ only in how radial matrices are represented.

### Details that matter

- **The contact density** is the Laplacian of the smeared Coulomb potential,
  ``\sum_k \alpha_k\,\delta_{\tau_k}(r)``, not ``\alpha_s(r)\,\delta_\sigma(r)``.
- **The tensor operator** uses the Pauli ``S_{12}``, so the Appendix A spin
  bracket carries a factor ``1/12``.
- **Eq. (A16)** contains only the confinement derivatives of the ``11`` and
  ``22`` pairs; the Coulomb derivative enters through (A15) alone.
- **The contact term acts in every partial wave.** Smearing gives it finite
  range, and Eqs. (23)–(26) of the paper keep it in P waves.

## Mixing: Section II D and Table III

| paper | content | implementation |
|---|---|---|
| text near Eq. (14) | tensor ``L = J\pm1`` mixing | [`tensor_mixing_components`](@ref), [`TensorMixing`](@ref) |
| text near Eq. (14) | antisymmetric spin-orbit ``{}^1L_J``–``{}^3L_J`` mixing | [`spin_orbit_mixing_components`](@ref), [`AntisymmetricSpinOrbit`](@ref) |
| Eq. (16) | general annihilation matrix element | `isoscalar_general_annihilation_solution` (`src/pseudoscalar_annihilation.jl`), [`isoscalar_annihilation_block`](@ref) |
| Eq. (17) | smeared wavefunction at the origin ``S_L(\Psi)`` | `_sL_smearing_factor` (`src/pseudoscalar_annihilation.jl`), from the momentum-space wave |
| Eq. (18a) | pseudoscalar `P1` amplitude | [`PaperP1Annihilation`](@ref) |
| Eq. (18b) | pseudoscalar `P2` amplitude, solved as a mass-dependent fixed point | [`PaperP2Annihilation`](@ref) |
| Table III | amplitudes ``A({}^3S_1) = 2.5``, ``A({}^3P_2) = -0.8`` and `P1`/`P2` constants | [`AnnihilationAmplitudes`](@ref), `[annihilation]` in the parameter file |

Diagonalization of every mixing block goes through one routine,
[`diagonalize_mixing_block`](@ref), which returns a [`MixingResult`](@ref).

## Decays: Section IV and Appendices B–D

| paper | content | implementation |
|---|---|---|
| Appendix B, Eqs. (B1)–(B36) | flavor, spin and spatial wavefunction conventions | `QuarkModelTransitions/src/flavor_algebra.jl`, `spin_algebra.jl` |
| Eq. (19) | pseudoscalar emission ``g\,\boldsymbol\sigma\cdot\boldsymbol q + h\,\boldsymbol\sigma\cdot\boldsymbol p'`` | [`PseudoscalarEmission`](@ref) (`pseudoscalar_emission.jl`) |
| Appendix C, Eqs. (C1)–(C3) | helicity amplitudes, partial waves, width | `partial_wave_projection` (`transition_amplitudes.jl`), [`decay_width`](@ref) |
| Table IV | single-Gaussian reduced amplitudes | the Table V reference backend (`strong_decays.jl`), used by GIPaper |
| Appendix D | E1, M1 and M2 photon transitions | [`PhotonEmission`](@ref), `photon_transition_class` (`photon_emission.jl`, `radiative_decays.jl`) |
| Eqs. (D7)–(D9) | leptonic widths | [`LeptonicCurrent`](@ref), `leptonic_decay_factor` (`annihilation_widths.jl`) |
| Table VII (a)–(c) | leptonic, two-photon and gluonic annihilation | [`LeptonicCurrent`](@ref), [`TwoPhotonAnnihilation`](@ref), [`GluonicAnnihilation`](@ref) |
| Table VII (d) | charge radii | [`charge_radius_squared`](@ref) |

## Central-potential variants

[`CentralPotentialMethod`](@ref) selects how the spin-independent potential is
built. Only one of them is the paper's; the rest exist to measure how much each
ingredient matters:

| method | construction |
|---|---|
| [`AppendixAMomentumSandwich`](@ref) | closed-form ``\tilde G`` and ``\tilde S`` with ``A(p)\tilde G A(p)``: **the paper's prescription**, used by the shipped parameters |
| [`AppendixAClosedForm`](@ref) | closed-form ``\tilde G + \tilde S`` without the momentum factor |
| [`PointwiseCentral`](@ref) | unsmeared ``G + S`` |
| [`Coulomb1DSmearing`](@ref), [`AppendixASmearing3D`](@ref), [`AppendixADerivativeG`](@ref) | diagnostic smearing approximations (mesh only) |

[`central_potential_path`](@ref) reports which one a parameter set uses.
