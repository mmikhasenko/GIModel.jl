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

## Linear coupling coefficients

```julia
c = coupling_coefficients(final, operator, initial;
    kinematics = CMKinematics(0.3), topology = :quark)
Cg = Dict(c.g.partial_wave_amplitudes)
Ch = Dict(c.h.partial_wave_amplitudes)
# For every allowed PartialWave pw: H[pw] = g*Cg[pw] + h*Ch[pw].
```

`topology` accepts `:both` (default), `:quark`, or `:antiquark` on both
`matrix_element` and `coupling_coefficients`. These are linear emitter-line
contributions, not separately normalized physical decay widths. The returned
objects are ordinary `TransitionAmplitude`s, retaining component provenance.

The [A-class study](../GIPaper/docs/a_class_universality/README.md) compares all
seven Table-IV A transitions with native waves. It explicitly records an
unresolved unequal-mass convention discrepancy: the existing final-wave
relative-coordinate gradient yields the **spectator** fraction divided by two
for common-beta SHO. The requested emitter-fraction extension is stored as a
separate, generally failing comparator; it has not been inserted into the
operator as an empirical correction.

For repeated surveys, both calls accept `integral_cache=Dict()`. This optional,
caller-owned cache reuses only spatial integrals with identical waves, orbital
labels, q, routing fraction, and couplings; it never combines flavor traces.
Do not mutate wave coefficients or mesh arrays while reusing that dictionary.
