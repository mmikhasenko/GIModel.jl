# Release claims audit

The active folder-by-folder checklist is [Release editing plan](release_cleanup_plan.md).
This document preserves the investigation baseline, not the remaining task list.

Follow-up: the authorized operator migration and rate-input tracing are documented
in the rate input ledger (retired during GIPaper cleanup). The findings
below describe the pre-cleanup baseline; observable ownership and the active
Lambda entry have since been corrected.

Audit date: 2026-09-26. Source baseline: `cfada91`. This is an investigation and cleanup basis, not a change to physics or package ownership. Source paths below are relative to the repository root. Paper-equation attribution here records the implementation's documented provenance; this audit does not independently rederive every equation from the original PDF.

## Claims and verdicts

| Claim | Verdict | Evidence and qualification |
|---|---|---|
| GIModel is the central reusable package | Confirmed | Both companion projects depend on GIModel. Its spectrum pipeline accepts model parameters, constituent masses, quantum numbers and solver controls without reference catalogs. |
| GIModel has few dependencies | Confirmed for direct runtime dependencies | Six: three external packages (KrylovKit, QuadGK, SpecialFunctions), three standard libraries (LinearAlgebra, Printf, TOML). FiniteDifferences is test-only. This is not a count of transitive packages. |
| GIModel has a clear API and internal pipeline | Substantially true for spectrum computation | `compute_spectrum → fixed_spectrum → fixed_channel_solution`, followed by `add_intra_meson_mixing`; native waves and coherent physical components are reusable. The exported surface also includes diagnostics, annihilation models and observable implementations. |
| Neither companion package is used by GIModel | Confirmed in current runtime source and project declarations | No companion imports, external includes or companion-data reads in `src/`. The parameter loader reads a caller-supplied TOML; its default points to GIModel's own `data/`. Comments mentioning GIPaper are not dependencies. |
| The two companions sit independently above GIModel | Needs refinement | GIPaper also depends on QuarkModelTransitions, importing `DecayChannel` for its Table V adapter. Actual dependencies: GIPaper → QuarkModelTransitions → GIModel, plus GIPaper → GIModel. |
| GIPaper holds reference numbers and comparison scripts | Confirmed, but broader than historical paper numbers | It also holds modern experimental kinematic masses, assignments, historical masses, calibration targets and comparison policies. |
| QuarkModelTransitions is a more advanced computation package | Confirmed, with maturity qualification | It has flavor/spin algebra, helicity/partial-wave projection, coherent composition and a working frozen Table V backend. Native Eq. (19) is an incomplete public workflow, not merely an advanced complete one. |
| The model needs only a handful of numerical inputs | True only with a clearly stated scope | Ordinary q/s/c/b spectra depend on 12 configurable scalar values and six fixed coupling-profile values. Optional flavor annihilation and observables add inputs. |
| No other ad-hoc numbers enter computation | Too strong as written | No per-resonance mass patch was found in the ordinary spectrum path. But additional fixed phenomenological inputs, calibration targets and numerical controls exist. They must be inventoried rather than all called either “parameters” or “ad hoc.” |
| Additional experimental constants are carefully book-kept | Partly true | The mass registry is explicit and fails on unknown assignments. Other observable/calibration numbers remain in source defaults and scripts. There is no single complete input ledger. |

## User value and release contract

GIModel is already useful for a user who wants meson spectra, contribution breakdowns, spin mixing, native wavefunctions and parameter studies. `GIParameters`, `Meson`, `RadialSolver` and `SpinTerms` distinguish model, flavor content, discretization and Hamiltonian content. `physical_components` exposes signed mixtures; `radial_wave` refuses to silently replace a mixed physical state with one precursor. These are strong foundations.

The simplest supported workflow is to load the shipped parameter set with `default_parameters_path()`, select a meson and levels, compute a spectrum, and inspect its states/waves. The root README's relative `data/...` example instead assumes the repository working directory. Its reproduction-oriented title, five deliverables, many audit links and “provisional” configuration obscure the small package entry point.

GIPaper is useful for checking a reconstruction against a specified reference and studying residuals. Modern PDG kinematics and historical GI predictions are different inputs: a modern-mass calculation is not automatically an exact reproduction of the original kinematics. Its mass-input README makes this distinction well, including unavailable assignments and uncertainty limitations.

QuarkModelTransitions is useful today for the frozen paper backend and advanced algebra work. Its README's `matrix_element(final, operator, initial)` sketch is not a runnable native-decay example. `pseudoscalar_emission.jl` implements internal S-to-S columns, but there is no specialized public `matrix_element` for `PseudoscalarEmission`; the generic method throws. `docs/phase4_spatial_checkpoint.md` explicitly lists composition, emitted-field identity, normalization conversion and higher orbital integrals as unfinished. Release documentation should label this capability experimental/incomplete or finish it before promising it.

## Exact spectrum input count

Count numerical scalar values, not TOML keys, aliases, integer algebra coefficients or solver controls. This count covers the shipped q/s/c/b parameterization, before optional flavor annihilation.

| Group | Values in computation | Count | Source and backward trace |
|---|---|---:|---|
| Constituent masses | 0.220, 0.419, 1.628, 4.977 GeV | 4 | `[masses]` → `quark_masses_from_raw` → `Meson` / `ConstituentMasses` → kinetic energies, smearing and momentum/spin factors. u/d/q are aliases of one light value, not three independent inputs. |
| Confinement | b=0.18 GeV², c=−0.253 GeV | 2 | `[potential]` → `ConfinementPotential` → static/smeared confinement and scalar spin-orbit kernels. |
| Universal smearing | σ₀=1.80 GeV, s=1.55 | 2 | `[relativistic_smearing]` → `contact_smearing_sigma` → Gaussian widths, closed-form potential/contact/spin kernels. |
| Relativistic exponents | εc=−0.168, εt=0.025, εso,v=−0.035, εso,s=0.055 | 4 | `[relativistic_factors]` → contact, tensor and vector/scalar momentum sandwiches. |
| Running-coupling weights | 0.25, 0.15, 0.20 | 3 | `src/constants.jl: ALPHA_COEFFS` → `alpha_s_r/q` and smeared Coulomb/contact/spin sums. |
| Running-coupling scales | 0.5, √10/2, √1000/2 GeV | 3 | `src/constants.jl: ALPHA_GAMMAS` → the same coupling and smeared kernels. |
| **Total** | **12 configurable + 6 fixed** | **18** | **14 potential/operator values plus four constituent masses.** |

For a particular unequal-flavor meson only two constituent mass values are used, so its corresponding accounting is 16 values; for equal flavor it is 15. These are input-value counts, not a claim that the historical fit had 18 statistically independent fitted parameters. The fixed Gaussian profile is an additional modeling prescription even though callers cannot vary it through `GIParameters`.

`Lambda_MeV = 200` is present in the shipped TOML but not loaded into the parameter object or consulted by the Gaussian running coupling. Do not count it as an active input. `alpha_s_critical = 0.60` appears in GIPaper's clean Table II data; in runtime code 0.60 is the sum of the three weights, not another input. Changing either reference key does not change the implemented coupling.

The backward computation path is:

1. `parameters.jl` / `quark_mass_table.jl` load the active configuration and masses.
2. `running_coupling.jl`, `smearing_appendix_a.jl` and `contact_hyperfine.jl` combine those inputs with the fixed Gaussian profile.
3. `hamiltonian.jl` or `harmonic_oscillator_basis.jl` assembles the relativistic central operator; `fixed_channel_solver.jl` adds contact and diagonal fine structure before diagonalization.
4. `spin_fine_structure.jl` and `spectrum.jl` build the remaining compatible mixing blocks. Branches use angular quantum numbers, masses and declared switches, not experimental resonance labels.
5. Optional flavor-annihilation APIs act afterward and require a separately declared model/amplitude/target policy.

No fitted spin-orbit/tensor bridge strength is present in this path; the loader explicitly rejects legacy `k_spin_orbit` and `k_tensor`. This is positive evidence, not a proof that every physical formula is correct.

## Additional physics inputs outside ordinary spectra

| Scope | Additional inputs | Location / consequence |
|---|---|---|
| Isoscalar annihilation | P1: 0.50 and 0.548 GeV; P2: 0.55 and 1.17 GeV; vector A=2.5; tensor A=−0.8 | Six values in `[annihilation]` and duplicated constructor/loader defaults in `parameters.jl`. P1 and P2 are alternatives; all six do not act on every result. `compute_spectrum` does not automatically consume them. |
| Calibrated isoscalar control | Four target masses: 0.520, 0.960, 1.440, 1.630 GeV | Targets belong to `GIPaper/src/table_iii_annihilation.jl`; GIModel requires explicit targets for `CalibratedP1Annihilation`. A fitted-to-target result must not be reported as an independent prediction. |
| Mock-meson electromagnetic overlaps | Fitted exponents 0.7 and 0.5 | Defaults in `mock_meson_overlaps.jl`. Higher-level M1/E1 wrappers call these defaults without exposing them as a common observable configuration. |
| Charge radius | Fitted exponent f=0.2 | `annihilation_widths.jl`, documented as fitted to π⁺. This is not a quadrature tolerance or universal constant. |
| Photon recoil approximation | β=0.40 GeV | `radiative_decays.jl`. A physical approximation scale, distinct from variational HO β. |
| Universal / unit conventions | αEM=1/137.036; GF=1.1663787×10⁻⁵ GeV⁻²; ℏc²=0.19733²; nucleon mass=0.93827 GeV | Constants in GIModel's observable files. No new spectrum fit parameters, but still numerical inputs with precision/provenance to document. |
| Native pseudoscalar emission | g, h, explicit constituent mass table | `QuarkModelTransitions/src/pseudoscalar_emission.jl`; no hidden g/h calibration defaults. The mass table is copied into the operator. |
| Frozen strong-decay model | A, S0, β; default calibration amplitudes +12.4 and −11.0 MeV¹ᐟ²; β=0.40 GeV | `calibrate_strong_decay_model` in `strong_decays.jl`. A/S0 are derived from anchors and breakup momenta, so do not count derived values and calibration inputs twice. `reproduce_table_v.jl` relies on these default anchors. |
| External kinematics | Input-state masses, lepton masses; separate fixed-wave mass/momentum corrections | Passed to observables or transition states; paper harnesses use GIPaper's registry. Model eigenvalues, mock masses and measured masses have distinct meanings. |
| Experimental validation anchors in scripts | fπ=0.1307 GeV; ℏ=6.582119×10⁻²⁵ GeV s; pion lifetime=2.6033×10⁻⁸ s; dilepton widths 7.04, 5.55, 2.34 and 1.34 keV | Hard-coded in `GIPaper/checks/audit_table_vii.jl:476–484`, outside the mass registry. These validate/compare observables and do not enter the mass Hamiltonian. Their dataset/edition must be recorded separately; the mass-registry README explicitly says legacy widths are not updated PDG 2026 widths. |

Exact rational/color/spin factors (4/3, 32π/9, Clebsch–Gordan coefficients, charge fractions), unit conversions and normalization factors are formula constants, not empirical tuning. Likewise 0.3 and 0.75 in Table IV reduced polynomials are documented formula coefficients; their decimal spelling alone is not evidence of a fit.

Additional fixed flavor assumptions deserve names: `_named_flavor_state` encodes η/η′ mixtures with coefficients 0.5 and ±1/√2. Those are flavor-state conventions, not a new arbitrary number for each channel. The synthetic masses 1.0, 0.1 and 2.0 in `algebraic_decomposition.jl` construct auxiliary states for the projection algebra; they are not experimental kinematic inputs to a predicted width.

## Numerical controls are another ledger

Finite-difference defaults are 450 points and rmax=24 GeV⁻¹, six retained levels, relativistic kinetic energy and dense diagonalization. HO defaults start at 24 basis states, step by eight up to at least 80, use two successive 10⁻⁴ GeV energy checks, a 0.25:0.10:2.35 GeV beta bracket and 0.002 GeV beta refinement. Optimized HO β is a variational output, not a fitted model constant.

Additional controls live below `RadialSolver`: HO operator quadrature has rtol=10⁻¹⁰, initial order 64 and cap 8192; narrow spin kernels use 10⁻⁸. Momentum transforms have range/sample defaults (including the observable 60 GeV cap and 900-point default). Small-r floors, Bessel-series cutoffs, mixing tolerances and P2 iteration damping/tolerance also exist. Their legitimacy depends on convergence evidence, not on hiding them from the input count.

Two release caveats follow. First, default FD is a fast profile, not the separately documented 2400-point/rmax=32 precision profile. Second, `ho_operator_matrix` warns and returns its last matrix on quadrature exhaustion, whereas basis convergence fails on exhaustion. An energy convergence certificate is not automatically a certificate for every quadrature or cancellation-sensitive observable. State the scope of certificates explicitly.

## Package boundaries and bookkeeping issues

1. **Dependency separation succeeds; conceptual ownership is inconsistent.** Root README says GIModel knows nothing about transition operators, but it exports M1/E1 amplitudes, radiative widths, leptonic/two-photon/gluonic observables and their constants. GIPaper README instead says GIModel owns amplitudes and widths. Choose one ownership contract before moving files. Generic radial/momentum operations should remain in GIModel; operator-specific physics can then be placed consistently.
2. **One private cross-package call exists.** QuarkModelTransitions calls `GIModel.radial_derivative_overlap`, which is defined but not exported. Its boundary test only checks selected transition exports and the absence of a reverse Project dependency. It does not enforce a public-only GIModel API. GIPaper has a qualified-name export check, but that is not a general dependency/data-isolation test.
3. **The parameter files do not have interchangeable semantics.** GIModel's active file says `provisional_raw_digitization`; GIPaper's clean file says `audited_active_solver_inputs` but intentionally lacks solver-prescription switches. Loading the clean file with GIModel's permissive loader selects pointwise central behavior and disables sandwiches. A user can load apparently cleaner numbers and silently compute a different Hamiltonian.
4. **Unknown/inactive keys can look effective.** Most loader keys are read selectively with fallbacks; only selected removed legacy keys are rejected. A typo or `Lambda_MeV` edit can appear accepted without affecting computation. Missing configuration sections can also change physics through fallback defaults.
5. **Mass bookkeeping is a good pattern to extend.** GIPaper pins a mass dataset, source hash/citation, selection, context-label assignments and averages. Unknown aliases error; known unassigned entries return `nothing`; historical masses are separate. Uncertainties are recorded but not propagated. Calibration amplitudes, observable exponents and non-mass experimental anchors need comparable records.
6. **Some case-specific policy is legitimate but must remain in the adapter.** `reproduce_table_v.jl` has named parent flavor sets, spectroscopic assignments and a charmed-only heavy-fraction policy; strange rows deliberately retain the paper's equal-mass form factor. That is reproduction policy, not universal quark-mass physics. Store these assignments/conventions as explicit data and keep generic kernels independent of display labels. The frozen Table V operator checks reference labels intentionally; do not treat it as the native general operator.
7. **Documentation trails the source.** QuarkModelTransitions README says spatial integrals begin in Phase IV although the first slice exists; the root README calls native Eq. (19) a future extension. Root verification prose says “both packages” although the script tests three. These are user-visible maturity/ownership inconsistencies.
8. **Release packaging needs a separate check.** No tracked LICENSE or `.github` workflow was found. None of the three Project files declares a Julia compat entry. Companion `[sources]` entries use local sibling paths; current tests in this checkout do not establish clean downstream installation. These are release tasks, not evidence of a physics error.

## Cleanup order and acceptance criteria

| Order | Work | Concrete acceptance criterion |
|---|---|---|
| 1 | Agree package scope and maturity | One consistent ownership table; native transitions explicitly marked incomplete; runnable minimal examples for released capabilities. |
| 2 | Establish the numerical-input ledger | Every empirical/default physics input has owner, units, source, role, configurability, consumers and calibration/validation status. Distinguish 18 spectrum inputs from optional extensions and numerics. |
| 3 | Make configuration unambiguous | A named audited production preset; inactive/unknown keys rejected or explicitly classified; reference-table data cannot silently masquerade as an equivalent production preset. |
| 4 | Strengthen boundaries | Core loads and runs with companion directories absent; companion runtime calls use the documented public API; no reference-data access in core computation. |
| 5 | Consolidate reference policy | Calibration anchors and assignments in GIPaper data with provenance; no hidden resonance-specific correction in generic kernels; separate calibration channels from validation channels. |
| 6 | Release verification and documentation | Declared Julia support, license, clean-install checks, repeatable tests and scoped convergence evidence; concise user-first README. |

The recommended next folder pass starts with root metadata, README and `data/`, then `src/` and `test/`, followed by QuarkModelTransitions and GIPaper. This resolves the shared contract and input definitions before changing downstream adapters.

## Verification record

Source inspection covered all three package declarations and entry points, the complete GIModel source numeric-literal inventory, input loaders and defaults, spectrum assembly/mixing paths, transition implementation and boundary tests, GIPaper mass/calibration adapters, and selected reproduction policies. Existing generated reports were consulted as historical evidence, not counted as fresh successful gates. No reference reports were regenerated and no computational code was changed.

Fresh runtime results:

| Check | Result |
|---|---|
| GIModel `Pkg.test()` | 1,137 / 1,137 passed, about 2m57s. Emitted a coarse-grid warning for a 30 GeV constituent test and an HO quadrature-exhaustion warning at L=0, beta=0.25, nbasis=24. |
| QuarkModelTransitions test entry point | 383 / 383 passed, about 32s. |
| GIPaper test entry point | 2,638 / 2,638 passed, about 31s. |
| Isolated GIModel copy | Copied only `src/`, `data/`, Project and Manifest into a temporary directory. Loaded GIModel and computed ten charmonium states at a deliberately small 120-point FD smoke-test resolution. Neither companion was loaded or present in that copy. This checks independence, not precision. |
| Inactive Lambda probe | Changing raw `Lambda_MeV` to 99999 produced an unchanged `GIParameters` object. |
| Active versus clean configuration | Active: `AppendixAMomentumSandwich` and all three sandwich/kernel switches true. Clean: `PointwiseCentral` and all three false. |
| Transition method inspection | Only specialized `TableVReference` and generic `StrongDecayOperator` methods exist for `matrix_element`; the native operator has no specialization. Derivative overlap export check returned false. |

Companion `Pkg.test()` attempts encountered sandbox restrictions on Julia launcher/depot lockfiles. Their existing test entry points were therefore run with the installed Julia 1.11.6 binary, `--compiled-modules=no`, using the resolved GIPaper environment. Both suites passed there; this does **not** establish an isolated QuarkModelTransitions dependency-resolution pass. No new packages were installed to bypass this limitation.

Full `scripts/verify_project.sh` regeneration is outside this investigative round; a unit-test pass alone is not a release certification or a new paper-reproduction measurement.
