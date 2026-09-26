# Differentiation for parameter fitting

Design update, 2026-09-26. These notes supersede earlier proposals in this
directory for new public prepared-energy APIs. They record the intended next
experiment, not implemented differentiation support.

## User objective and API

The objective is to tune quark-model parameters against predictions: trial
parameters → recompute states and observables → evaluate a loss against data →
differentiate the loss → update parameters. Preserve the existing flow through
`GIParameters`, `Meson`, `fixed_spectrum` / `compute_spectrum`, and
`spectrum_state`. The initial experiment targets fixed-sector masses; it does
not promise differentiation of the full mixed spectrum.

Do not introduce `prepare_fixed_channel`, `fixed_channel_energy`, or a mandatory
parameter-vector representation merely to accommodate AD. Numerical preparation
and derivative rules should initially be internal. Reuse the same Hamiltonian
assembly for ordinary and differentiated calculations. Parameter packing, tied
quark masses, and optimization coordinates can live in the fitting example or
adapter; experimental data need not become a dependency of the model package.

## Hellmann–Feynman and the changing state

For an isolated normalized eigenstate, differentiating
`E = ψ† H ψ` gives `dE = ψ† (dH) ψ`: the two wavefunction-response terms cancel
by `Hψ = Eψ` and `d(ψ†ψ) = 0`. This is not a frozen-wavefunction approximation.
Every new trial parameter point requires a new eigenstate solve. The spectral
rule avoids explicitly computing the eigenvector derivative for this energy
gradient only.

For a mass loss `L = Σᵢ ((Mᵢ - targetᵢ)/sᵢ)^2`, propagate the mass derivatives
through the residuals and weights. Widths, radii, transition amplitudes, and
mixing matrices built from parameter-dependent radial states generally also
need wavefunction response. Applying a spectral rule to the final mixing
matrix does not remove the need to differentiate its entries through those
radial states. Never label a stopped-wave derivative as the total fitting
gradient.

The current HO solver continuously refines beta within a bracket and optimizes
the highest requested radial level. The older grid-winner/margin discussion is
insufficient for this implementation. A fixed-beta derivative and a derivative
of the reoptimized prediction are different contracts. The envelope theorem
removes the beta-response term only for the optimization objective at a smooth
stationary optimum, not automatically for other levels or observables. AD
through comparison-driven search steps is not implicit differentiation of that
optimum. Preserve the requested channel level count when comparing predictions.

## Next experiment: Enzyme with the existing design

1. Use an existing fixed-spectrum scalar mass loss with explicitly fixed beta
   and basis settings. Test Enzyme on the actual assembly and solver path before
   undertaking broad refactoring; pin the tested Julia and package versions.
2. Address demonstrated blockers internally. If the eigensolver needs a custom
   derivative, implement the spectral rule at that boundary. Keep useful
   mutation and share parameter-independent factorizations and quadrature data;
   do not freeze mass-dependent factors or cache predictions across parameters.
3. Extend checks to active spin terms, independent constituent masses, and a
   tied equal-flavor mass. The existing Hellmann–Feynman probes establish a
   numerical identity, not successful reverse AD or a working fitting loop.
4. Validate gradients against finite differences of the same loss, recomputing
   predictions at each perturbed point. Sweep finite-difference step sizes and
   independently check gradient convergence with basis and quadrature resolution.
5. Demonstrate a small mass fit before expanding the public API. Define adaptive
   solver and state-tracking semantics before claiming derivatives of reoptimized
   or mixed-spectrum predictions; add wave-dependent observables separately.

Enzyme does not require dual-number storage: ordinary `Float64` fields and
arrays are not inherently barriers to its reverse mode. This does not eliminate
activity, shadow-memory, external-call, or eigensolver support constraints.
Genericize types only for a demonstrated numerical or backend requirement.
See the [Enzyme FAQ](https://enzyme.mit.edu/julia/stable/faq/) and
[custom-rule documentation](https://enzyme.mit.edu/julia/stable/generated/custom_rule/).

Acceptance is a correct and numerically converged gradient of the actual
prediction loss, with explicit treatment of crossings and unresolved gaps.
Agreement for a derivative of an internal Hamiltonian alone is an intermediate
check. The paper caveat is recorded in `report/main.tex`, under “Derivatives”.
