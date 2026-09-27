# Model inputs

The ordinary q/s/c/b spectrum uses **18 numerical physics inputs**. The reviewed
runtime preset is returned by `default_parameters_path()`; its historical filename
`data/parameters.provisional.toml` is retained for compatibility. It is a runnable
model configuration, unlike the Table II transcription in GIPaper.

| Inputs | Count | Shipped values | Authoritative runtime location |
|---|---:|---|---|
| Constituent masses q (u=d), s, c, b | 4 | 0.220, 0.419, 1.628, 4.977 GeV | preset `[masses]` (stored in MeV) |
| Confinement b, c | 2 | 0.18 GeV², −0.253 GeV | preset `[potential]` |
| Smearing sigma0, s | 2 | 1.80 GeV, 1.55 | preset `[relativistic_smearing]` |
| Relativistic exponents c, t, vector SO, scalar SO | 4 | −0.168, 0.025, −0.035, 0.055 | preset `[relativistic_factors]` |
| Gaussian coupling weights | 3 | 0.25, 0.15, 0.20 | `src/constants.jl`: `ALPHA_COEFFS` |
| Gaussian coupling scales | 3 | 0.5, √10/2, √1000/2 GeV | `src/constants.jl`: `ALPHA_GAMMAS` |

The six Gaussian values are fixed by the published approximation; the other
12 are configurable. This table documents their values rather than supplying
another runtime source. The asymptotic coupling 0.60 is the sum of the weights.
Historical `Lambda_MeV=200` belongs to the paper's motivation of that profile;
it is not evaluated by the runtime and is rejected in a runtime parameter file.

TOML loading requires all ordinary model fields, including explicit central-method
and spin switches. Unknown fields, incomplete sections and nonfinite numbers
are errors. Metadata is descriptive and does not enter the Hamiltonian.
Programmatic constructors still support deliberately simplified diagnostic models.
`load_quark_masses` reads only the masses section, so it can also read Table II
reference data; `load_parameters` cannot silently interpret that reference file
as a different Hamiltonian.

The optional `[annihilation]` section contains six additional phenomenological
values for **isoscalar mass mixing**. These are outside the ordinary 18-number
contract. Omitting the whole section uses `AnnihilationAmplitudes()` defaults;
supplying it requires all six fields. These are not decay-width operators.

Meson flavors, quantum numbers, solver choice, basis sizes, mesh boundaries,
integration tolerances and convergence criteria are separate computation choices.
An HO energy convergence certificate applies to the requested energies; it does
not certify every derived overlap or observable. Quadrature warnings likewise
need assessment for the integral concerned, even if energy checks pass.

Decay operators and their additional inputs belong to QuarkModelTransitions;
GIPaper owns reference comparisons and external assignments. See the
[rate input ledger](../GIPaper/docs/rate_input_ledger.md) and
[release editing plan](release_cleanup_plan.md).

Transition matrix elements and generic widths use the masses stored in the
resolved states. Comparisons with other masses or momenta use
`QuarkModelTransitions.mass_correction_factor` separately at fixed wavefunctions;
see the [transition guide](../QuarkModelTransitions/README.md#mass-and-momentum-comparisons).
