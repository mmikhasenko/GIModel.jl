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

## Supported matrix elements

The generic `matrix_element(final, operator, initial)` API currently knows three
operator classes:

| Operator class | Process represented | Result |
| --- | --- | --- |
| `PseudoscalarEmission` | Solver-native GI Eq. (19) strong decay `M_i -> M_f + P`, where `P` is an elementary emitted `J^P=0^-` meson | All allowed helicity and partial-wave amplitudes, coherently summed over mixed-state components; `decay_width` gives the partial width |
| `TableVReference` | One explicitly supplied Godfrey--Isgur Table V strong-decay row | The row's preselected partial-wave amplitude and width, for reproduction and audit rather than prediction |
| `PhotonEmission` | M1, E1, or M2 radiative transition using the published Appendix-D mock-meson realization of the Eq. (22) current | A coherently composed `RadiativeAmplitude`; `decay_width` gives the partial width in MeV |

The radiative primitives remain public for inspecting individual overlap
kernels. The generic photon path is:

```julia
operator = PhotonEmission(quark_masses; recoil_order = 0)
amplitude = matrix_element(final_state, operator, initial_state)
decay_width(amplitude)
```

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

The package also computes the following observables through specialized
functions; these have not yet been adapted to `matrix_element`:

- lowest-order `1S0 -> 2g`, `3S1 -> 3g`, `3P2 -> 2g`, and `3P0 -> 2g`
  annihilation amplitudes and widths;
- pseudoscalar and tensor two-photon annihilation amplitudes;
- pseudoscalar leptonic, vector dilepton, and axial-vector tau-decay factors
  and widths.

It cannot yet compute general quark-pair-creation decays (for example
`a1 -> K* Kbar` through `s sbar` creation), emission of a vector, scalar,
axial, or tensor meson, general hadronic or semileptonic weak transitions, or
continuum-induced mixing/coupled-channel pole dressing. In particular,
`PseudoscalarEmission` is an elementary one-pseudoscalar emission operator,
not a `3P0` pair-creation model.

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
