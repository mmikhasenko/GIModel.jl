# Julia AD backend research

Research checked on 2026-08-06. Links below are primary project documentation or
source. Backend behavior changes, so the eventual decision must be backed by
version-pinned probes in this repository.

## Project-specific conclusion

Do not select a backend merely from the dependency list of a machine-learning
framework. Lux supports several AD systems because it is infrastructure used by
many downstream workloads. GIModel needs a much narrower contract.

The likely architecture is:

1. expose separate prepared FD and fixed-beta HO Hamiltonian-assembly kernels;
2. give eigenvalues an explicit symmetric spectral derivative;
3. test Enzyme as the primary compiler-level backend;
4. use ForwardDiff and finite differences as small-kernel/reference checks after
   generic numeric storage exists;
5. keep a backend-neutral public API only if two implementations actually pass
   the same derivative tests.

## Does Enzyme require parametric numeric types?

No. Enzyme differentiates the compiled `Float64` primal program at LLVM level.
Its API marks immutable scalar/struct inputs as `Active` and mutable storage as
`Duplicated`. It does not inject a `Dual` or `TrackedReal` into the primal field.

Parametric numeric types remain worthwhile for:

- ForwardDiff and other operator-overloading backends;
- `Float32`/`BigFloat` numerical checks;
- correct promotion and fewer hidden conversions;
- making the numerical kernel a conventional generic Julia API.

They should be justified as numerical cleanup and multi-backend support, not as
an Enzyme prerequisite.

## Backend matrix for this codebase

| Backend/tool | Mutation | LAPACK/eigen issue | Needs generic scalar/storage types | Role here |
|---|---|---|---|---|
| Enzyme | strong support for in-place code | only part of Julia linear algebra is covered; explicit eigen rule still preferable | no for `Float64` primal | primary native spike on both routes |
| Mooncake | designed for mutation and existing Julia code | foreign calls need rules; test actual eigen path | no Dual injection | secondary reverse-mode spike |
| ForwardDiff | local mutation works only with compatible generic storage | cannot propagate through non-Julia code; custom spectral boundary needed | yes | leaf/JVP reference and small-parameter Jacobians |
| Zygote | array mutation is a major limitation | relies heavily on rules | usually generic out-of-place code | compatibility target only after rules isolate mutation |
| ReverseDiff | tracked values/tapes; mutation restrictions remain | explicit boundary needed | generally yes | low priority; compiled tapes risky with value-dependent control flow |
| Tracker | tracked values and mature ML compatibility | explicit boundary needed | generally yes | low priority |
| ChainRulesCore | not an AD engine | defines `frule`/`rrule` boundaries | no | lightweight place for spectral rules if adopted |
| DifferentiationInterface | backend-neutral operations and tests | delegates to backend | depends on backend | evaluation harness/public facade only after native experiments |
| FiniteDifferences | no AD trace | treats solver as black box | no | correctness oracle, already a test dependency |

## Enzyme

[Enzyme.jl](https://enzyme.mit.edu/julia/stable/) supports forward and reverse
mode over statically analyzable LLVM and explicitly supports mutating functions.
This fits the local matrix-assembly style in GIModel.

Important caveats from the
[Enzyme FAQ](https://enzyme.mit.edu/julia/stable/faq/):

- only some Julia linear algebra is supported because LAPACK and other libraries
  contain foreign calls;
- common BLAS operations have derivative support, but other external routines
  require rules;
- temporary storage holding active values must have a corresponding shadow;
- activity-unstable control flow can require runtime activity and extra care;
- combined reverse mode has restrictions on non-scalar returns, for which an
  in-place output or split reverse mode may be appropriate.

Consequences here:

- dense `A * G * A` products are a good target;
- `eigen(p2)`, HO quadrature eigensystems, and `eigen(H)` should not remain opaque
  calls in the AD trace;
- prepare route-specific constant eigensystems outside AD and own the
  `H -> eigenvalues` derivative explicitly;
- start with native Enzyme APIs because the project needs precise control of
  active/constant/prepared arguments.

## ForwardDiff

[ForwardDiff](https://juliadiff.org/ForwardDiff.jl/stable/) injects `Dual`
numbers into the target computation. Its
[limitations page](https://juliadiff.org/ForwardDiff.jl/stable/user/limitations.html)
states that the target and any storage used by it must be generic, and that
derivatives cannot propagate through non-Julia code.

The audit snapshot's `Float64` parameter fields and matrix allocations lost Dual
information on both routes. Parametric leaves and generic central-Hamiltonian
storage are now implemented; a ForwardDiff experiment would still need prepared
route kernels and the explicit spectral boundary described in the design notes.

For approximately 15--25 continuous model parameters and a vector of predicted
masses, ForwardDiff can still be a useful Jacobian implementation or oracle.
The expensive solve is repeated by chunk, however, so it should be benchmarked
against a reverse/VJP strategy rather than assumed optimal.

ForwardDiff does not natively consume ChainRules rules. The
[ChainRules documentation](https://juliadiff.org/ChainRulesCore.jl/stable/)
mentions `ForwardDiffChainRules.jl` as a bridge. That extra integration should be
added only if the project actually promises both interfaces.

## Mooncake

[Mooncake](https://chalk-lab.github.io/Mooncake.jl/stable/) is a source-to-source
reverse-mode system designed with first-class mutation support. Its documentation
emphasizes low-level rules for intrinsics, built-ins, and foreign calls, as well
as rule testing.

It is a credible secondary experiment because mutation in GIModel is mostly
local array construction. It is not automatically easier than Enzyme around
LAPACK. The same explicit spectral boundary remains desirable.

## Zygote

[Zygote's limitations](https://fluxml.ai/Zygote.jl/stable/limitations/) identify
array mutation and foreign calls as common blockers and recommend custom
ChainRules rules when operations cannot be made pure.

Rewriting GIModel's local assembly to satisfy Zygote everywhere would impose a
large design cost for limited project benefit. Once Hamiltonian construction and
the spectral primitive are isolated, Zygote compatibility may become cheap, but
it should not drive the first transformation.

## ReverseDiff and Tracker

These operator-overloading reverse-mode systems use tracked numeric values.
They encounter many of the same concrete-storage barriers as ForwardDiff.

The
[DifferentiationInterface backend notes](https://juliadiff.org/DifferentiationInterface.jl/DifferentiationInterface/stable/explanation/backends/)
also warn that a compiled ReverseDiff tape records the branches taken at
preparation time. Reusing it after value-dependent control flow changes can
silently produce wrong results. GIModel has state sorting, beta selection, and
some parameter-dependent branches, so a tape cannot span the orchestration
layer safely.

Neither backend offers enough project-specific advantage to justify early work.

## ChainRules and the eigensystem boundary

[ChainRulesCore](https://juliadiff.org/ChainRulesCore.jl/stable/) is a lightweight
rule-definition package, not an AD engine. It allows a package to define
forward rules (`frule`) and reverse rules (`rrule`) without depending on Zygote
or another full backend.

ChainRules already contains symmetric/Hermitian eigenvalue rules in
[`symmetric.jl`](https://github.com/JuliaDiff/ChainRules.jl/blob/main/src/rulesets/LinearAlgebra/symmetric.jl).
The eigenvector formulas divide by eigenvalue differences. Its more general
factorization rules explicitly list degenerate-matrix support as a TODO in
[`factorization.jl`](https://github.com/JuliaDiff/ChainRules.jl/blob/main/src/rulesets/LinearAlgebra/factorization.jl).
The implementation cites Mike Giles's
[matrix derivative report](https://people.maths.ox.ac.uk/gilesm/files/NA-08-01.pdf).

GIModel can use the simpler eigenvalue-only formula and define a narrower
primitive than the full `eigen` factorization. The same primitive applies to FD
and fixed-beta HO matrices; only assembly differs. This avoids promising
eigenvector derivatives before their state/gap semantics are settled.

## DifferentiationInterface

[DifferentiationInterface](https://juliadiff.org/DifferentiationInterface.jl/DifferentiationInterface/stable/)
provides common JVP, VJP, gradient, and Jacobian operations plus preparation and
testing infrastructure. Its backend documentation notes that Enzyme activity
and multiple-argument handling are not fully represented by the wrapper and
recommends trying Enzyme's native API when that matters.

Recommended use:

- native APIs during the first Enzyme/Mooncake experiments;
- `DifferentiationInterfaceTest` to compare supported operators later;
- a DI-based public entry point only after the required activity/context model
  is demonstrated, not before.

## Lux as a comparison, not a template

Lux's current
[`Project.toml`](https://github.com/LuxDL/Lux.jl/blob/main/Project.toml) uses core
AD abstractions directly and exposes several full AD engines through weak
dependencies/extensions. Its
[AD support page](https://lux.csail.mit.edu/stable/manual/autodiff) gives the
strongest support to Enzyme, Zygote, and ForwardDiff, with other backends at
lower support tiers.

That architecture is appropriate for an ML framework whose users choose their
own AD stack. GIModel should initially carry at most a light rule abstraction or
one optional backend extension. Supporting many engines is not itself a project
goal.

## Required backend experiment

For each candidate backend, test matched prepared FD and fixed-beta HO problems
and record:

1. derivative operator supported: directional derivative/JVP, VJP, full Jacobian;
2. agreement with finite differences and analytic `dE/dc = 1`;
3. agreement with the explicit Hellmann--Feynman result;
4. compile/preparation time, steady-state time, and allocations;
5. behavior for one unequal-mass sector and several orbital channels on each
   route;
6. failure behavior near an intentionally constructed small eigenvalue gap;
7. exact package and Julia versions.

A backend is not “supported” because one scalar derivative runs once, nor does
passing FD imply passing HO. It is supported for a named route when this matrix
passes with documented tolerances and failure modes.
