# GIModel API

Every exported name of `GIModel`, grouped by topic. The Manual pages explain how the pieces fit together.

```@contents
Pages = ["gimodel.md"]
Depth = 2:2
```

## Parameter types and loading

The model inputs: see [Inputs and parameters](@ref).

```@docs
load_parameters_and_quark_masses
load_parameters
load_quark_masses
default_parameters_path
GIParameters
ConfinementPotential
RelativisticSmearing
RelativisticFactors
FineStructure
AnnihilationAmplitudes
QuarkMassTable
```

## Central potential methods

How the spin-independent potential is built.

```@docs
CentralPotentialMethod
AppendixAMomentumSandwich
AppendixAClosedForm
PointwiseCentral
Coulomb1DSmearing
AppendixASmearing3D
AppendixADerivativeG
central_potential_method
CentralPotentialPath
central_potential_path
central_potential_values
```

## Quarks and mesons

```@docs
Meson
ConstituentMasses
reduced_mass
flavor_label
is_equal_flavor
AbstractQuark
LightQuark
StrangeQuark
HeavyQuark
charge
flavor_symbol
mass_GeV
```

## Levels and quantum numbers

```@docs
BasisState
spectrum_levels
FineStructureMultiplet
orbital_angular_momentum
orbital_label
```

## Spectrum stages and states

See [Computing a spectrum](@ref).

```@docs
compute_spectrum
fixed_spectrum
add_intra_meson_mixing
central_spectrum
Spectrum
MixedSpectrum
CorrectedSpectrum
CentralSpectrum
MixedState
CorrectedState
CentralState
spectrum_state
parameters
```

## Mixing records and blocks

```@docs
StateMixing
MixingMechanism
AntisymmetricSpinOrbit
TensorMixing
IsoscalarAnnihilation
MixingBlock
MixingResult
diagonalize_mixing_block
same_j_mixing
spin_orbit_mixing_components
tensor_mixing_components
```

## Annihilation mixing

See [Isoscalar flavor mixing](@ref).

```@docs
compute_isoscalar_spectrum
add_isoscalar_annihilation
PaperP1Annihilation
PaperP2Annihilation
CalibratedP1Annihilation
MomentumIntegralSmearing
isoscalar_annihilation_block
pseudoscalar_annihilation_block
annihilation_basis_input
pseudoscalar_annihilation_basis_input
fix_annihilation_phase
fix_annihilation_phase!
isoscalar_pseudoscalar_annihilation_solution
isoscalar_general_annihilation_solution
isoscalar_general_s1_solution
```

## Solver types

See [Solvers and convergence](@ref).

```@docs
RadialSolver
OscillatorSolver
FiniteDifferenceSolver
SpinTerms
OscillatorConvergence
numerics_provenance
channel_solution
fixed_channel_solution
ChannelRadialSolution
RadialChannelKey
SectorComputation
solve_sector
resummed_channel_solution
contact_hyperfine_nonperturbative_states
```

## Wave representations and operations

See [Wavefunctions](@ref).

```@docs
RadialWave
OscillatorWave
MeshWave
radial_wave
physical_components
physical_state_amplitude
physical_transition_amplitude
sample_wave
wave_norm
radial_expect
radial_overlap
radial_derivative_overlap
wave_mean_squares
MomentumWave
OscillatorMomentumWave
MeshMomentumWave
momentum_wave
momentum_expect
momentum_overlap
momentum_functional
mock_momentum_wave
```

## Potentials and spin operators

Building blocks of the Hamiltonian, mainly for diagnostics and tests.

```@docs
alpha_s_q
alpha_s_r
contact_smearing_sigma
spin_dot
LdotS
tensor_triplet_LJ
tensor_triplet_offdiag_sameJ
fine_structure_split
fine_structure_components
fine_structure_radial_kernels
fine_structure_grid_operator
radial_cross_expect_udr
ho_p2_matrix
ho_r2_matrix
ho_operator_matrix
```

## Internal names referenced above

These functions are not exported. They appear here because exported docstrings
refer to them.

```@docs
GIModel.contact_hyperfine_shift
GIModel.fine_structure_grid_matrices
GIModel.ho_fine_structure_matrices
GIModel.orthonormalize_physical_basis
```

```@docs
convergence
```
