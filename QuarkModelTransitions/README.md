# QuarkModelTransitions.jl

Bare quark-model transition amplitudes built from solved states and wave
functions supplied by `GIModel.jl`.

This package owns transition operators, external-state adapters, channels,
helicity and partial-wave representations, amplitude composition, and width
normalizations. `GIModel.jl` remains the Schrödinger solver and overlap engine;
it does not depend on this package.

```julia
using GIModel
using QuarkModelTransitions
```

The primary user path is deliberately small:

```julia
initial = physical_state(spectrum, label)
final = TwoMesonChannel(first, second)
operator = PseudoscalarEmission(g, h, quark_masses)
amplitude = matrix_element(final, operator, initial; kinematics = OnShell())
partial_waves(amplitude)
decay_width(amplitude)
```

The package also contains the frozen Godfrey--Isgur Table IV/V reproduction
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

These APIs moved from GIModel into this package: M1/E1/M2 amplitudes and recoil,
gluonic and two-photon annihilation, leptonic factors and widths, mock-meson
operator overlaps, and charge radii. Existing callers now add
`using QuarkModelTransitions`; qualified calls use `QuarkModelTransitions`.
Generic radial/momentum transforms and overlaps, physical-state composition,
and annihilation mixing of the mass spectrum remain in GIModel.

`ELECTROMAGNETIC_DEFAULTS` and `STRONG_DECAY_DEFAULTS` expose the additional
phenomenological inputs. They are separate from spectrum parameters and from
numerical controls. GIPaper owns the table-specific reference data and the
[rate input ledger](../GIPaper/docs/rate_input_ledger.md).
