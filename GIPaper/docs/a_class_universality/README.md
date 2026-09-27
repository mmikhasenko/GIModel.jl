> Report refresh, 2026-09-27: imported from `codex/a-class-universality`
> (363fc54). These archived kernels intentionally use central P/D waves and
> contact-resummed S waves; the all-L contact fix does not affect that
> prescription. `verify_a_class_artifacts.py` checks the archived traces and
> `plot_a_class.jl` renders them without overwriting their computation hashes.
> A full new `study_a_class.jl` calculation requires the study branch's
> `coupling_coefficients` and topology-selection API, not yet present in this
> checkout's QuarkModelTransitions. The report refresh uses the verified archive.

# A-class realistic-wave universality study

This deliverable tests nine reduced Table-IV transitions (seven A, one A′, and
one A″) using the existing
native Eq. (19) operator. **It is a formal common-q, unmixed-basis kernel survey,
not a list of physical open decays.** The numerical results and plots are in
[numerical_summary.md](numerical_summary.md). The publication target is the
[threshold histogram](threshold_histogram.pdf), with its [caption](threshold_caption.md).
The [finite-q histograms](histograms.pdf) and [nine transition panels](transition_panels.pdf)
are validation figures.

There is an unresolved physics convention discrepancy: the requested
`m_emitter/[2(m_emitter+m_spectator)]` SHO coefficient is not produced by the
current operator. We preserve that fact, including a failing-comparator flag,
rather than alter Eq. (19) to force the requested value. Numerical certification
below applies to the implemented kernel only. It is not certification of the
requested unequal-mass extension.

## API audit and implementation

The starting implementation already had `PhysicalState`, `TwoMesonChannel`,
`PseudoscalarEmission`, `CMKinematics`, native radial/derivative overlaps,
angular projection, and coherent component assembly. Committed tests covered
SHO equal-mass identities, generic orbital integrals, coherent state mixing,
HO/FD spatial agreement, and the frozen Table IV/V backend. There was no public
linear-coupling decomposition, topology selector, seven-class census, or
trace-to-plot ledger. GIModel needs no new domain concept or solver change.

`coupling_coefficients(final, operator, initial; kinematics, topology=:both)`
returns `(g=TransitionAmplitude, h=TransitionAmplitude)` through two unit-coupling
calls to `matrix_element`. The latter now also accepts `:quark` and `:antiquark`;
the default is unchanged. Both complex helicities and partial waves recombine
as `g*Cg+h*Ch`. The test covers arbitrary signed couplings, unequal masses,
complex q, and topology additivity. No symbolic algebra or alternate amplitude
implementation is introduced.

## Census and meaning of a trace

| ID | Parent → surviving daughter + elementary P | Relative L |
|---|---|---:|
| A1 | 1³S₁ → 1¹S₀ + P | 1 |
| A2 | 1³P₂ → 1¹S₀ + P | 2 |
| A3 | 1³P₂ → 1³S₁ + P | 2 |
| A4 | 1³P₁ → 1³S₁ + P | 2 |
| A5 | 1¹P₁ → 1³S₁ + P | 2 |
| A6 | 1³D₃ → 1¹S₀ + P | 3 |
| A7 | 1³D₃ → 1³S₁ + P | 3 |
| A8 (A′) | 1¹P₁ → 1³P₀ + P | 1 |
| A9 (A″) | 1³D₃ → 1¹P₁ + P | 2 |

The finite applicability rule is explicit in `flavor_routes()`: the parent
ordered pair is ud/du (isovector); uu/dd/ss (resolved isoscalar basis); us/ds/su/sd
(strange); or c/b with each u/d/s, in both orientations. Each light emitter may
become u, d, or s; the spectator stays unchanged and the emitted flavor is fixed
by that transfer. This gives 90 flavor/topology routes, hence 810 reduced
transition routes. Heavy constituents never emit light pseudoscalars; cc, bb,
and cb/bc are excluded. Heavy-light records retain the original strange versus
nonstrange light constituent and the final exact flavor pair separately.

Neutral uu/dd/ss emitted fields and isoscalar parents are **resolved basis
components**, not eta/eta-prime or omega/phi mixing predictions. Isovector
neutral coherent superpositions are not counted again alongside their resolved
components. All mixing coefficients in this basis survey are explicitly one;
no physical pole mixing is inferred. The API can handle coherent physical
states, but that would be a distinct follow-up study. Flavor multiplicities
are correlated bookkeeping, not statistical independence or probabilistic
weights. Charge-conjugate and isospin-degenerate entries remain separate.

S waves use the central Hamiltonian with contact hyperfine resummed independently
for singlet and triplet. P/D waves use the native central radial solution; their
spectroscopic spin labels enter the transition algebra. Fine-structure-induced
radial distortions, annihilation/flavor mixing, tensor mixing, and antisymmetric
spin-orbit mixing are not included. State masses are eigenvalues of precisely
these wave Hamiltonians, not PDG masses. The elementary field's 0.14 GeV mass is
an explicit placeholder unused by the fixed-q kernel; it must never be used to
infer which channels open.

## Unequal-mass convention audit

The paper's Eq. (19), journal p. 200 (PDF page 12), specifies a gradient on the
final wave and quark/antiquark signs. The scan was checked against the local
paper page image; the OCR line is incomplete. The equal-mass SU(6) reduction
does not by itself specify a unique beyond-paper moving unequal-mass extension.

In the current implementation set `r=r_q-r_qbar`, `alpha=m_s/(m_e+m_s)`, and
`eta=-1` for quark emission, `+1` for antiquark emission. The plane wave is
`exp(eta*i*alpha*q*z)`. For common-beta ground-state SHO,

```
F = <S|exp(eta*i*alpha*q*z)|S> = exp[-alpha²*q²/(4*beta²)]
<z exp(eta*i*alpha*q*z)> = eta*i*alpha*q/(2*beta²) F
∂z psi_S = -beta²*z psi_S
```

Together with the Eq. (19) topology signs, these give `Ch/Cg=alpha/2`.
Projection onto the maximum A-class relative wave gives the same ratio for all
seven classes; independent Gaussian tests enforce this for both emitter lines.
Thus the implementation gives

```
g + (1/2) * m_s/(m_e+m_s) * h,
```

whereas the task's proposed extrapolation gives the emitter fraction. They
coincide only at equal masses. Changing the coordinate fraction or adding a
boost term merely to enforce the emitter expression would be an unverified
operator change. This study makes neither change. `requested_emitter_hg`,
`implemented_spectator_hg`, `SHO_hg_re/im`, and `requested_formula_pass` expose the
difference row by row for the original A family. `mass_formula_applicable=false`
marks the two new families; their unequal-mass comparators are evaluated directly
through the same SHO operator, with no guessed mass rescaling.
`equal_mass_SU3_hg` stores the family-specific threshold value (+1/4, −1/4, +1/8);
no strange mass is secretly replaced.

A future resolved moving-frame derivation should update the operator, the
independent unequal-mass tests, and regenerate the study together. Until then,
heavy-light findings are conditional on the implemented convention.

## Added A′ and A″ transitions

The daughter orbital label and total J are explicit in the transition catalog.
This is essential: channel spin equals the daughter's total J (the emitted P has
J=0), not its constituent spin. A′ uses a central ³P₀ daughter, and A″ a central
¹P₁ daughter. P-wave fine-structure radial distortions are still excluded under
the same unmixed central-wave prescription as the existing parent P/D waves.

The equal-mass common-beta **threshold** ratios are −1/4 for A′ and +1/8 for A″.
The full native Eq. (19) evaluator has finite-q polynomial corrections for these
P-wave-daughter channels even with SHO waves. Tests therefore check their
extrapolated ratios, and independently compare finite-q analytic SHO overlaps
with sampled-mesh derivative/overlap integrals. They do not force the threshold
constants at finite q. This corrects the earlier discussion's suggestion that
all three families have constant SHO ratios away from threshold.

The publication histogram uses each family's own matched SHO threshold limit,
so U_hg=1 is the common reference, including for the negative A′ coefficient.
All nine transitions enter, with the five original sector colors. The A family
has seven times the transition multiplicity of either new family; this weighting
is explicit, not an equal-family average. Sector counts are 108 isovector, 162
isoscalar, 216 strange, 162 charmed and 162 bottom-flavored combinations.

## Numerics and certification

The common-beta SHO comparator is an actual `OscillatorWave(L,0.4,[1])`
evaluation through the same public API, flavor/angular projection, constituent
masses, topology, and q. There is no Gaussian removed from the native result.

The nonzero threshold sequence is 0.12, 0.06, 0.03, 0.015, 0.0075 GeV.
Both coefficients are divided by q^L; adjacent pairs give q² Richardson
extrapolations, and the last two extrapolations must agree to 0.2% separately
for g and h. All five uncollapsed samples are retained. q=0 in the derived
file means extrapolation; the operator is never evaluated at q=0 for this
reduction. The separate common-q scan uses 0.15, 0.30, 0.60 GeV.

The light pilot now covers all nine transitions before the census expands. The initial
HO default versus FD600/1200 at rmax=28 GeV⁻¹ had a worst coefficient difference
of 1.929% in A7. Increasing the HO basis to 48–128 with a 1 keV energy tolerance
and enlarging the FD box to 48 GeV⁻¹ did not remove it. The retained pilot tables
record this negative convergence finding; they are not silently overwritten.
The final FD comparison uses 1200 and 1800 interior points at rmax=48 GeV⁻¹.

Consequently, the **3% HO/FD gate is calibrated after the pilot**, not a claim of
sub-percent precision. FD refinement must also be below 0.3%. The exact measured
errors are stored for every common-q/threshold route. Failed points remain in
the comparisons and plots (crosses in the panels). A deviation comparable to the
HO/FD discrepancy is not evidence of a physical wave effect. The residual pilot
discrepancy has not been localized; interpreting percent-level effects requires
further solver/operator quadrature investigation. The study supports only
larger deviations that survive the recorded numerical envelope.

For native v=(Cg,Ch) and matched SHO s, the primary statistic is
`U_hg=(Ch/Cg)_native/(Ch/Cg)_SHO`, with reference one. `Rg` and `Rh` retain complex
phases. If a component is below 1e-8 times its vector norm, scalar ratios are
flagged unstable. The fallback compares unit coupling vectors after aligning
their overall complex phase; a zero vector is flagged as having no direction.
Real and imaginary parts are exported explicitly, never discarded implicitly.

## Artifacts and reproducibility

Run from the repository root using the existing scripts environment:

```sh
julia --project=GIPaper/scripts GIPaper/scripts/study_a_class.jl --pilot
julia --project=GIPaper/scripts GIPaper/scripts/study_a_class.jl
julia --project=GIPaper/scripts GIPaper/scripts/plot_a_class.jl
```

Tables are committed as lossless `.csv.gz` files (`gzip -dc FILE.csv.gz` to read);
the runner also writes uncompressed local copies. The plotter accepts either.

- `traces.csv`: one raw operator piece per exact route, q and backend. The readable
  `trace_id` includes sector, parent/daughter/emitted flavors, transition, both
  bases, relative L, topology, q, backend and piece. Written before reductions.
- `waves.csv`, `wave_samples.csv`: native solver settings, convergence record,
  eigenvalues, parameter hash, exact HO expansion coefficients or FD samples.
  Exact mass-degenerate solves may be shared; raw flavor traces are not merged.
- `comparisons.csv`: complex Cg/Ch, matched SHO, Rg/Rh/U_hg, vector fallback,
  formula-discrepancy flags, and all source trace ids (including every threshold
  sample). Primary plots use HO; both FD alternatives remain here for audit.
- `certificates.csv`: per-route solver, refinement and threshold errors.
- `plot_membership.csv`: every histogram-bin contribution and every panel point,
  including coincident entries, joined directly to raw trace ids. Histogram
  bins count routes once at each q; they do not average them. No failed numerical
  certificate is filtered out. Unstable scalar entries have bin zero and are
  represented by the vector-fallback artifact if any occur.
- `run.toml`: census, grids, tolerances, Julia version, parameter SHA-256 and source file hashes.
- `threshold_histogram.pdf/png`: the single-panel publication target, with
  consistently typeset LaTeX subscripts. `threshold_caption.md` records scope
  and numerical qualifications; `threshold_bins.csv` stores the exact bin edges.
- `histograms.pdf/png`, `transition_panels.pdf/png`: finite-q validation figures
  and preview rasters. Colors encode the same five sectors throughout.

Table V's on-shell inputs, reference backend, amplitudes and tests are unchanged.
The report deliberately calls no decay width routine. The study is not a new fit
of g/h, a reproduction of Table V's correction column, or a physical branching
fraction prediction.

`python3 GIPaper/scripts/verify_a_class_artifacts.py` independently reconstructs
all exported complex coefficients (including threshold extrapolations) from raw
trace ids, recomputes the ratios, and checks the census and every plot membership.
It reads the committed gzip tables using only Python's standard library.

The fixed-q operator depends on constituent masses, not the external mass labels.
HO/SHO records share external reference masses; FD records retain their own
solver eigenvalues for audit. Those small external-mass differences do not enter
any coefficient or comparator in this off-shell study.

## Nine-transition threshold result

The primary histogram contains all 810 combinations. Of these, 798 pass the
unchanged numerical gates. All 90 A′ and all 90 A″ threshold combinations pass;
the twelve flags are the original A6/A7 cases. A′ spans U_hg=0.825275–1.174724,
and A″ spans 0.880702–1.104238. No scalar ratio is unstable at threshold.
The original seven-transition results were compared by trace id against the
preceding archive: all 40,320 original operator-piece values are exactly
unchanged. The new families broaden the low side of the distribution without
altering the original calculation.

## Original seven-transition baseline and remaining work

The following numbers describe the original A-only baseline. The current
nine-transition results, including per-family ranges and flags, are generated
in [numerical_summary.md](numerical_summary.md).

The threshold range is 0.921537–1.221715 across all sectors. Equal-mass light
kernels can therefore depart appreciably from the common-beta coupling ratio
without invoking any unequal-mass correction. The A3/A4/A5 ratios coincide for
fixed flavor routing in this wave prescription because those transitions share
the same central P-wave and triplet S-wave radial inputs; the angular coefficients
change the overall amplitude, not this ratio. These entries remain separate.

Of 2,520 primary HO estimates, 2,484 pass the numerical gate. All 630 points at
q=0.60 GeV pass. Twelve routes fail at each of threshold, q=0.15, and q=0.30 GeV:
A6/A7 ending in ss, starting either in ss or a light-strange pair. The maximum
HO/FD coefficient discrepancy is 3.8748%, versus 0.0213% for FD grid refinement.
Those failures remain plotted and must not be interpreted as established
physical deviations. Locating the residual HO/FD difference is still open.
The maximum threshold extrapolation discrepancy is 1.95e-7 (relative), and no
scalar ratio requires the small-denominator fallback in this run.

The requested emitter-fraction formula fails for all 420 unequal-mass parent
routes at threshold; it agrees for the 210 equal-mass routes. This is a
convention/implementation finding independent of native-wave convergence.
Consequently this deliverable does **not** claim to have validated the task's
specified unequal-mass formula. It provides the operator audit, reproducible
countercheck, and complete conditional survey needed to resolve that issue.

## Paper integration

The publication histogram has an inset upper-right sector legend. Transition
families, the nine-transition/810-combination census, and numerical qualifications
belong in the caption. The separate report repository stores its wrapper and
caption alongside the other figure sources in `sources/a_class_threshold.jl`
and `sources/a_class_threshold_caption.tex`; the wrapper calls this study's
renderer and copies its vector PDF to `figures/a_class_threshold.pdf`. It accepts
an optional GIModel.jl checkout path, so an isolated computation worktree can
supply the paper without duplicating the numerical implementation.
