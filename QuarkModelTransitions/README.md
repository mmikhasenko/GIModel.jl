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
amplitude = matrix_element(final, operator, initial; kinematics = OnShell())
partial_waves(amplitude)
decay_width(amplitude)
```

The package also contains the frozen Godfrey--Isgur Table IV/V reproduction
backend. Its `StrongDecayModel`, `DecayChannel`, and scalar factorization
helpers are compatibility/research-audit concepts, not a second recommended
user workflow.

Phase III supplies the wave-independent algebra for native GI Eq. (19): sparse
flavor contractions, constituent-spin matrix elements, `L dot S`
spectroscopic recoupling, a shared spatial-integral basis, complete compressed
helicity coefficient vectors, and Appendix-C partial-wave projection. Phase IV
now supplies `PseudoscalarEmission(g,h,masses)` and one general
`L_i,m_i -> L_f,m_f` spatial evaluator for analytic SHO and native HO/FD waves,
including real or complex off-shell momentum. End-to-end physical-state
composition and relativistic normalization remain in progress.
