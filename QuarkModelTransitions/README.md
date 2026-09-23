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

The current implementation contains the frozen Godfrey--Isgur Table IV/V
reference backend and the resolved-state/kinematics foundation for native
operators. The next implementation stage derives spin, angular, flavor, and
helicity-to-partial-wave algebra before adding numerical GI Eq. (19) matrix
elements.
