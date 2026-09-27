# QuarkModelTransitions.jl

Bare quark-model transition amplitudes built from solved states and wave
functions supplied by `GIModel.jl`.

This submodule owns transition operators, external-state adapters, channels,
helicity and partial-wave representations, amplitude composition, and width
normalizations. `GIModel.jl` provides the Schrödinger solver and overlap engine and includes
this module. Dependencies are managed by the repository-root `Project.toml`;
this folder has no separate package environment.

## Start with states calculated by GIModel

The primary workflow is `compute_spectrum` → `physical_state` → `matrix_element`
→ `decay_width`. The transition calculation uses the calculated masses, radial
wavefunctions, and the components/mixing retained by that spectrum. No Gaussian
replacement of the solved waves is made.

Run each self-contained block in the GIModel environment (`julia --project=.`
from the repository root). Load the transition API with
`using GIModel.QuarkModelTransitions`.
The finite-difference grids below keep the examples quick; refine the solver and
check the observable's convergence before quoting quantitative predictions.
The parameter file and selected spectrum treatment remain model assumptions.

### Photon transition with calculated charmonium waves

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = PhotonEmission(masses)
amplitude = matrix_element(final, operator, initial)
width_MeV = decay_width(amplitude)
@assert width_MeV > 0
```

### Pseudoscalar emission with calculated light-meson waves

This evaluates `a1+ → rho+ + pi0`, using GIModel radial solutions and masses.
The neutral pion's explicit flavor components are combined in the isospin limit;
this is a stated flavor assumption, not an inferred neutral-state mixing fit.
The ordered daughters are `(surviving, emitted)`.

`g` and `h` are additional direct and recoil couplings in GeV^-1. They are not
predicted by the spectrum. The values below only demonstrate the API: no coupling
fit is performed, so the resulting width is not a calibrated prediction. A
quantitative study must state how the couplings were calibrated with its chosen
wavefunctions and external masses.

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
a1, vector, pion_basis = BasisState(1, "P", 3, 1), BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
charged = compute_spectrum(params, Meson(masses, :u, :d);
    levels=[a1, vector, pion_basis], solver)
initial = physical_state(charged, a1)
rho = physical_state(charged, vector)
# Resolve pi0 = (u ubar - d dbar)/sqrt(2) in the isospin limit.
up = compute_spectrum(params, Meson(masses, :u, :u); levels=[pion_basis], solver)
down = compute_spectrum(params, Meson(masses, :d, :d); levels=[pion_basis], solver)
u, d = physical_state(up, pion_basis), physical_state(down, pion_basis)
pion = PhysicalState("pi0", (u.mass_GeV + d.mass_GeV)/2,
    [(basis=c.basis, coefficient=sign*c.coefficient/sqrt(2), wave=c.wave)
     for (state, sign) in ((u, 1), (d, -1)) for c in state.components];
    provenance=(source=:isospin_combination,))
final = TwoMesonChannel(rho, pion)
g, h = 0.7, 0.3 # Illustrative GeV^-1 inputs, not fitted GI predictions.
operator = PseudoscalarEmission(g, h, masses)
amplitude = matrix_element(final, operator, initial)
@assert partial_waves(amplitude) == [PartialWave(0, 1), PartialWave(2, 1)]
width_MeV = decay_width(amplitude)
@assert width_MeV > 0
```

The emitted pion is treated as an elementary field by `PseudoscalarEmission`.
Its mass and flavor components enter; the radial integrals involve the parent
and surviving daughter. Using calculated GI waves improves those inputs within
this emission model; it does not turn the operator into a three-meson
pair-creation calculation.

## Optional: assumed oscillator waves

For analytic checks or an approximate calculation, construct `PhysicalState`
objects with chosen waves and masses instead. Here both states use a single
oscillator function with an assumed width of 0.5 GeV; no GI spectrum is solved.

```julia
using GIModel, GIModel.QuarkModelTransitions
wave = OscillatorWave(0, 0.5, [1.0])
initial = PhysicalState("psi", 3.10, [(
    basis=BasisState(1, "S", 3, 1; flavors=(:c, :c)),
    coefficient=1.0, wave=wave,
)])
final = PhysicalState("eta_c", 2.98, [(
    basis=BasisState(1, "S", 1, 0; flavors=(:c, :c)),
    coefficient=1.0, wave=wave,
)])
operator = PhotonEmission(QuarkMassTable("c" => 1.628))
amplitude = matrix_element(final, operator, initial)
width_MeV = decay_width(amplitude)
@assert width_MeV > 0
```

`OscillatorWave(L, beta, coefficients)` stores a normalized expansion in
harmonic-oscillator functions. Here `[1.0]` selects the lowest radial function
for the chosen orbital angular momentum `L`; `beta` sets its width in GeV.

## Supported transitions and observables

Use `?matrix_element` to construct an amplitude, `?decay_width` for its physical
width in MeV, and `?mass_correction_factor` for a separate fixed-wave comparison.
Each function has one help entry covering its supported argument combinations.

**Initial states:** `PhysicalState` obtained with `physical_state(spectrum, level)`.

**Final states for matrix elements:** `PhysicalState` for photon emission,
`TwoMesonChannel(surviving, emitted)` for pseudoscalar emission, `Vacuum()` for
a meson-current matrix element, `TwoPhotonChannel()` for two photons.

**Final states for widths:** `PhysicalState`, `TwoMesonChannel`,
`LeptonNeutrinoChannel(m_lepton; ckm=...)`, `MasslessLeptonPair()`,
`TwoPhotonChannel()`, `TwoGluonChannel()`, `ThreeGluonChannel()`.

**Operators:** `PhotonEmission`, `PseudoscalarEmission`, `LeptonicCurrent`,
`TwoPhotonAnnihilation`, `GluonicAnnihilation` (widths only).

| Operator | Matrix-element final state and result | Width final state |
| --- | --- | --- |
| `PhotonEmission` | Daughter `PhysicalState`; `RadiativeAmplitude` | Same daughter |
| `PseudoscalarEmission` | Ordered `TwoMesonChannel`; `TransitionAmplitude` | Same channel |
| `LeptonicCurrent` | `Vacuum()`; dimensionless reduced current in `AnnihilationAmplitude.value` | `LeptonNeutrinoChannel` for pseudoscalars, `MasslessLeptonPair` for vectors |
| `TwoPhotonAnnihilation` | `TwoPhotonChannel()`; width amplitude in GeV^(1/2) in `AnnihilationAmplitude.value` | `TwoPhotonChannel()` |
| `GluonicAnnihilation` | No differential or helicity amplitude is exposed | `TwoGluonChannel()`, `ThreeGluonChannel()` |

The generic width functions return **MeV** and include their prescribed phase
space, spin averages, and symmetry factors. M1 matrix elements are moments in
nuclear magnetons; E1/M2 matrix elements are width amplitudes in MeV^(1/2).
The two-photon amplitude has width `1000abs2(result.value)` in MeV.
A leptonic current requires an explicit physical final channel to obtain a width.
Gluonic rates currently require a single flavor-diagonal component.
Axial tau widths remain available through the scalar `axial_tau_width` in GeV.

For reference reproduction, qualified `QMT.TableVReference` remains a compatibility
operator accepting `ReferenceState`/`TwoMesonChannel` inputs and returning a
`TransitionAmplitude`; it is separate from the native operators above. Here
`QMT` abbreviates `QuarkModelTransitions`.

### Leptonic decay from a calculated state

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = LeptonicCurrent(:electromagnetic, masses)
current = matrix_element(Vacuum(), operator, initial)
width_MeV = decay_width(MasslessLeptonPair(), operator, initial)
@assert width_MeV > 0
@assert isfinite(current.value)
```

The electromagnetic vector current derives charge and angular coefficients from
explicit flavor components. Other current kinds require `AnnihilationTerm`
coefficients. Pseudoscalar weak widths also require the CKM magnitude explicitly;
the vector width uses the massless-lepton approximation.

## Mass and momentum comparisons

`matrix_element` and `decay_width` use the masses stored in their input states.
Use `mass_correction_factor` for a comparison
at another external mass or momentum, keeping wavefunctions, mixing, constituent
masses, and operator settings fixed. This does not recompute a spectrum.

| Operator | Comparison inputs | The factor multiplies |
| --- | --- | --- |
| `LeptonicCurrent`, `TwoPhotonAnnihilation` | `target_mass` in GeV | Current/amplitude `.value` |
| `PhotonEmission` | `target_momentum` in GeV | Photon amplitude/moment `.value` |
| `PseudoscalarEmission` | `target_momentum` in GeV, `partial_wave` | Selected partial-wave amplitude |
| `GluonicAnnihilation` | `target_mass` in GeV | Integrated width |

```julia
using GIModel, GIModel.QuarkModelTransitions
import GIModel.QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = TwoPhotonAnnihilation(masses, AnnihilationTerm((:c, :c), 4/9))
result = matrix_element(TwoPhotonChannel(), operator, final)
factor = mass_correction_factor(TwoPhotonChannel(), operator, final; target_mass=2.98)
comparison_amplitude = factor * result.value
comparison_width_MeV = 1000abs2(comparison_amplitude)
@assert factor ≈ (2.98/final.mass_GeV)^1.5
```

The two-photon example includes the integrated phase space in its amplitude.
That is not universal: leptonic and M1/strong-emission widths have additional
mass or momentum factors. Do not treat every squared amplitude correction as
a full width correction. Gluonic factors already multiply the width directly.
Emission comparisons reevaluate momentum-dependent kernels. A zero reference
amplitude, or closed reference emission channel, has no multiplicative ratio.

The radiative primitives remain available as internal helpers for inspecting
individual overlap kernels. The photon example above uses the primary API.

The initial and final spectroscopic states select `DirectM1`, `HinderedM1`,
`AllowedE1`, `SpinFlipE1`, or `SpinFlipM2`; the user does not separately name
the multipole. Explicit `:u` and `:d` flavor components determine the current.
An averaged `(:q,:q)` component is intentionally rejected because it cannot
distinguish `e_u+e_d` from `e_u-e_d`. `recoil_order=0` means the leading term
for the selected class; order 2 adds the implemented relative `(qr)^2` M1
correction. The returned amplitude records the selected class and order, and
`verbose=true` reports that selection while evaluating it.

Paper-reproduction scripts may supply a package-internal resolved current for
legacy coarse flavor states. Such row-specific information belongs to the
paper audit and is not part of the ordinary transition API.

It cannot yet compute general quark-pair-creation decays (for example
`a1 -> K* Kbar` through `s sbar` creation), emission of a vector, scalar,
axial, or tensor meson, general hadronic or semileptonic weak transitions, or
continuum-induced mixing/coupled-channel pole dressing. In particular,
`PseudoscalarEmission` is an elementary one-pseudoscalar emission operator,
not a `3P0` pair-creation model.

The submodule also contains the frozen Godfrey--Isgur Table IV/V reproduction
backend. Its `StrongDecayModel`, `DecayChannel`, and scalar factorization
helpers are compatibility/research-audit concepts, not a second recommended
user workflow.

`PseudoscalarEmission(g,h,masses)` evaluates GI Eq. (19) end to end. It derives
sparse flavor and spin contractions, recouples the spectroscopic states,
evaluates one general `L_i,m_i -> L_f,m_f` kernel on analytic SHO or native
HO/FD waves, constructs the compressed helicity vector, and projects every
allowed partial wave. Physical-state components interfere coherently and
Eq. (C2) supplies the normalization-aware width in MeV.

The elementary emitted field must be a resolved `J^P=0^-` state. Components
must carry explicit ordered flavors; the averaged `(:q,:q)` sector is rejected
until the caller resolves its physical `u ubar`/`d dbar` isospin combination.
`TwoMesonChannel` is ordered. For `PseudoscalarEmission`, construct it as
`TwoMesonChannel(surviving, emitted)`; the second state must be `J^P=0^-`.
For two pseudoscalars, reverse the arguments to select the other GI Fig.-14
assignment. The API never silently sums both descriptions. The frozen Table
IV/V backend remains a reference calculation and is never selected by
`PseudoscalarEmission`.

## Radiative and annihilation observables

These APIs live in the transition submodule: M1/E1/M2 amplitudes and recoil,
gluonic and two-photon annihilation, leptonic factors and widths, mock-meson
operator overlaps, and charge radii. Existing callers now add
`using GIModel.QuarkModelTransitions`; qualified calls use `QuarkModelTransitions`.
Generic radial/momentum transforms and overlaps, physical-state composition,
and annihilation mixing of the mass spectrum remain in GIModel.

`ELECTROMAGNETIC_DEFAULTS` and `STRONG_DECAY_DEFAULTS` expose the additional
phenomenological inputs. They are separate from spectrum parameters and from
numerical controls. GIPaper owns the table-specific reference data.

## Public API and migration

The 21 supported names cover states and final channels, concrete operators,
`AnnihilationTerm` flavor coefficients, `PartialWave`, `partial_waves`,
`matrix_element`, `decay_width`, `mass_correction_factor`, and
`charge_radius_squared`. All are exported.

Constants, individual kernels, amplitude containers, normalization tags,
abstract operators, and reference-table helpers are internal. Qualified access
and explicit imports do not imply a supported API. Reference audits retain
explicit imports of the internal kernels they need.

`matrix_element` evaluates on-shell transitions without a momentum override.
Use `mass_correction_factor(...; target_momentum=...)` for momentum corrections;
strong emission also requires `partial_wave`. A ratio is undefined for a closed
or zero reference amplitude.

The [transitions guide](https://mmikhasenko.github.io/GIModel.jl/dev/manual/transitions)
and the [API reference](https://mmikhasenko.github.io/GIModel.jl/dev/api/transitions)
cover every supported name. `scripts/audit_documentation.jl` executes all
docstring and README examples and checks help links; it runs with the tests.
Use `?matrix_element` and `?decay_width` as the navigation entry points.

`decay_width` returns **MeV** for transition and reference amplitudes. Specialized
gluonic, leptonic, and `m1_radiative_width` functions return **GeV**;
`two_photon_amplitude` is in **GeV^(1/2)**. Charge radii are in **GeV^-2**;
use `QuarkModelTransitions.HBARC_FM2` to convert to fm².
