# Help and public API policy

Help follows the GIModel discoverability work introduced in commit `e5a2246`:
explain purpose, inputs, units, how to obtain an object, and what to do next.

Every supported binding (all are exported) has a source docstring
with a self-contained Julia example and explicit `## Related` links.
`matrix_element`, `decay_width`, and `mass_correction_factor` each have one consolidated docstring covering
all dispatches. Examples import their dependencies, construct concrete inputs,
and include checks of useful properties. They require
no prior notebook cells, external paper data, or current working directory.
Primary transition examples start with GIModel spectra and adapt their solved
states with `physical_state`. Modest example grids require convergence checks
for quantitative work, and strong couplings are explicitly illustrative inputs.
Manual oscillator waves are an optional approximation; primitive-kernel examples
may use them as small analytic demonstrations, not reproduced measurements.

Use `?matrix_element`, `?decay_width`, and `?mass_correction_factor` to enter the graph. Useful paths include:

- `physical_state` → `PhysicalState` → `matrix_element` → `decay_width`.
- `PseudoscalarEmission` → `matrix_element` → `partial_waves`.
- `PhotonEmission` → `matrix_element` → `mass_correction_factor`.
- `LeptonicCurrent` → `Vacuum` → `matrix_element` → `decay_width`.
- `TwoPhotonAnnihilation` → `TwoPhotonChannel` → `decay_width`.
- `GluonicAnnihilation` → `TwoGluonChannel` → `decay_width`.

## Namespace revision (0.2)

The 21 supported names are exported: state and final-channel constructors,
concrete operators, flavor inputs, partial-wave selection, and observable methods.
Reference-table helpers, primitive kernels, constants, result containers,
abstract extension types, normalization tags, and inferred photon tags are internal.
Explicit imports used by reference audits do not make those helpers supported API.
The [generated inventory](discoverability_graph.md) lists the supported names.

## Verification

From the repository root:

```sh
julia --project=. QuarkModelTransitions/scripts/audit_documentation.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

The first command regenerates the inventory and graph. The same audit runs in
the package tests without rewriting files. It discovers all public names from
the loaded module (not a hand-picked list), checks live help, checks every source
public docstring for an example and related links, executes each code block in a
fresh module, and validates package and GIModel link destinations. Every public
entry must reach every other entry; merely having incoming/outgoing links would
miss disconnected clusters. README examples run in fresh modules too.

The source scanner recognizes this package's ordinary triple-quoted docstrings
and declarations. If macro-generated public bindings or another documentation
syntax are introduced, extend the scanner. The checks cover help content and
example execution, not the rendering of links in Pluto or a documentation site.
