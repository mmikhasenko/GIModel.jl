# Scope and limitations

This page lists what GIModel covers and what it deliberately leaves out. The
purpose is to make clear which questions the package can answer.

## Covered

- **The complete 1985 spectrum algorithm**: the fixed-sector Hamiltonian with
  all Appendix A operators, spectroscopic mixing, and isoscalar annihilation
  mixing, solved in the paper's oscillator basis with automatic convergence.
- **An independent solver**: finite differences, as a cross-check.
- **Physical states**: masses, contributions, and signed components through
  every mixing stage, with solver-native wavefunctions in position and
  momentum space.
- **All observables of the paper's Tables III–VII**, computed from the solved
  states: isoscalar mixing, strong decays in the paper's reduction, radiative
  transitions, and leptonic, two-photon and gluonic widths, and charge radii.
- **Beyond the paper**: the Eq. (19) emission operator evaluated on the solved
  wavefunctions rather than a common Gaussian, with all partial waves and
  coherent mixing (see [Strong decays beyond the paper](@ref)).

## Not covered

- **Other decay models.** There is no quark-pair-creation (``{}^3P_0``) model.
  Decays that need a new quark pair (for example ``a_1 \to K^*\bar K``),
  emission of mesons other than pseudoscalars, and hadronic or semileptonic
  weak transitions are not implemented.
- **Coupled channels.** Mass shifts and mixing induced by decay channels
  (threshold effects, pole dressing) are not included. This matters for states
  near or above open-flavor thresholds, such as the ``X(3872)``.
- **Isospin breaking.** As in the paper, ``u`` and ``d`` are degenerate.
- **Electromagnetic decays of mixed isoscalars.** Photon emission and leptonic
  widths of states from [`compute_isoscalar_spectrum`](@ref) are not available
  through the public API, because the nonstrange component is labeled `:q`
  (see [Isoscalar flavor mixing](@ref)). Two-photon widths are available.
- **Baryons and exotics.** Only ``q\bar q`` mesons.
- **Top quarks.** The 1985 parameter set has no top mass; toponium rows of
  the paper are not computed.

## Things to keep in mind

- **Mixing depends on the requested basis.** Only requested levels enter a
  mixing block ([Computing a spectrum](@ref)). Include every level that can mix
  with the state you study.
- **Convergence is certified for energies.** Wavefunction-sensitive quantities
  (widths, values at the origin, cancellation-prone overlaps) should be checked
  separately ([Solvers and convergence](@ref)).
- **Kinematics use model masses.** Decay momenta follow from the masses of the
  solved states. Use [`mass_correction_factor`](@ref) to compare at measured
  masses.
- **Couplings of the emission operator are inputs.** ``g`` and ``h`` are not
  predicted and must be calibrated.
- **Known reference discrepancies.** The mixing angles in the 1985 figure
  captions ([Mixing angles: a paper erratum](@ref)) and the amplitudes of
  transitions from radially excited ``\eta`` states
  ([Results at a glance](@ref)) differ from the paper for reasons that are
  documented, not hidden.
