# Numerical inputs to the rate tables

This ledger describes the computation that the current reproduction scripts actually execute. A “paper parameter”, a measured mass, a fitted transition exponent, an angular coefficient and a quadrature tolerance are different kinds of input. None should be hidden behind the statement “zero parameters.” Reference predictions used only to score a result are not inputs to that prediction.

The operator implementations now belong to **QuarkModelTransitions**. **GIModel** supplies spectra, signed state components and generic wave operations; its isoscalar annihilation terms remain there because they change masses and mixing. **GIPaper** selects reference identities, external kinematics, calibration prescriptions and comparison targets. Its dependence on both packages is intentional.

## Where Lambda came from

`Lambda_MeV = 200` originated in the Table II extraction (commit `42b1a9f`, then copied into the first spectrum configuration in `6d71a35`). It was not invented for a particular resonance or a decay calculation.

The local paper transcription, Sec. II Eqs. (11)–(13) and the Fig. 2 caption, gives its role:

- Eq. (11) is the perturbative QCD running coupling, with Λ.
- Eq. (12) replaces it with a saturating sum of Gaussians suitable for the model.
- Fig. 2 compares that approximation against QCD with Λ=200 MeV and publishes the fitted profile `0.25 exp(−Q²) + 0.15 exp(−Q²/10) + 0.20 exp(−Q²/1000)`, Q in GeV. The caption specifies that Λ refers to the Nf=2 regime of the threshold-continuous comparison curve.

The implementation directly evaluates this published profile; it does not fit it anew from Λ. Thus Λ is **upstream historical provenance of the coupling approximation**, not an additional runtime degree of freedom. It has been removed from GIModel's active TOML, remains in GIPaper's historical/clean Table II records, and does not enter any rate table at runtime. The three weights and three scales do enter the spectrum and any runtime `alpha_s_q` call. The saturation value 0.60 is their sum, not a seventh independent number.

Evidence: `paper/vision_ocr/godfrey_isgur_1985_vision_ocr.md` lines 163–200 and 294; `GIPaper/data/raw/digitized_tables/table_ii_parameters/table_ii_parameters.csv`; `src/constants.jl`; `src/running_coupling.jl`. Paper context here was checked in the local transcription; the runtime conclusion is independently established by the code.

## Executable per-row records

Run:

```sh
julia --project=GIPaper/scripts GIPaper/scripts/trace_rate_inputs.jl
```

This executes the existing V, VI and VII pipelines, regenerates their reports, and writes:

- `input_traces/common.toml`: exact active spectrum numbers, fixed Gaussian profile, operator defaults, universal constants, numerical controls and source SHA-256 hashes.
- `input_traces/table_v.toml`: every canonical row's actual flavor-spin coefficient, reduced-amplitude class, orbital power, masses and provenance, heavy fraction, external angle, fitted A/S0 and computed results.
- `input_traces/table_vi.toml`: all 79 rows, masses, photon momentum, signed physical components **including native wave coefficients**, relevant exponents, recoil switches, additive moment and result.
- `input_traces/table_vii.toml`: row assignments, charges/flavor coefficients, masses, wave-treatment policy, computed S/Mtilde/factors/mixing coefficients where returned by the pipeline, amplitudes/radii and derived dilepton widths. All 61 canonical rows are listed, including the four unsupported top rows.

These are numerical dependency records, not an instruction-by-instruction execution dump. Derived eigenvectors/intermediates are labeled separately from independent inputs; the source hashes identify the exact implementation of algebraic constants and lower-level numerical rules. `common.toml` is a shared inventory, not a statement that every number acts on every row. The sections below specify active subsets. Missing experimental inputs stay unavailable, with no model-mass fallback.

## Shared spectrum inputs, where waves are used

The shipped mass set is q=u=d=0.220, s=0.419, c=1.628, b=4.977 GeV. Central waves use b=0.18 GeV², c=−0.253 GeV, σ₀=1.80 GeV, s=1.55 and the coupling weights (0.25,0.15,0.20) with scales (0.5,√10/2,√1000/2) GeV. Contact-distorted S waves additionally use εc=−0.168. Fixed triplet P sectors additionally use εt=0.025, εso,v=−0.035 and εso,s=0.055. A central-only wave does not consume these spin exponents.

The active central prescription is the Appendix-A momentum sandwich; contact and fine-structure sandwiches and smeared fine-structure kernels are enabled. Physical quantum numbers and the chosen set of mixing partners are discrete physics inputs. Solver basis size, beta optimization, quadrature and convergence tolerances are numerical inputs. The variational beta and physical mixing coefficients are calculated outputs, not fitted transition inputs.

## Table V: strong-decay amplitudes and their squares

Graph: **two calibration amplitudes + calibration masses + beta → A,S0; row masses → q; row coefficient/class/orbital power/heavy fraction → pure amplitude; optional external angle → physical amplitude; amplitude² → partial width**.

This frozen backend **does not use GIModel eigenvalues or solved wavefunctions**. Its GI constituent masses enter the charmed recoil fraction only. The companion “realistic factors” audit is separate.

| Input | Exact current value/source | How it enters |
|---|---|---|
| rho→pi pi calibration | +12.4 MeV¹ᐟ² | Determines A from its breakup momentum and suppressed spatial factor. |
| B→(omega pi)S calibration | −11.0 MeV¹ᐟ² | Determines S0; B here is the historical b1 label, not a bottom meson. |
| Oscillator scale | beta=0.40 GeV | q/beta and Gaussian recoil form factor. |
| Calibration masses | Registry context `V`: rho, pi, B, omega | Change q at the anchors and therefore change A/S0 globally. Exact values and mass definitions are emitted in the trace. |
| Row masses | Registry `V`; `V:1^3F_4` for the delta F-wave parent | Determine each physical channel's q. |
| Flavor-spin coefficient | Each CSV row's numeric `coefficient` | Used literally, including transcription rounding (e.g. 1.1547005); it is not recomputed from the printed symbolic formula. |
| Reduced class / orbital power | CSV `amp_class`, `qbar_power` | Select algebra and power of q/beta. |
| Convention | `LeadingS0()` | S,D,P reduce to S0; the optional Table-IV polynomial subtractions are not used in this report. A′, A″ and A0 use A. |
| Recoil fraction | 0.5 normally; c/(c+d)=1.628/(1.628+0.220) for charmed parents | Changes the Gaussian; Ac rows with the appropriate power also carry recoil. Strange rows deliberately retain 0.5 as a paper-reproduction convention. |
| External mixing angles | Strange 1P: +34°; strange 1D: +33°; charmed 1P: −41° | Rotates singlet/triplet amplitudes at the physical parent's q. These angles are prescribed by the harness, **not calculated by GIModel**. |

The exact spatial factor, including MeV conversion, is in `QuarkModelTransitions/src/strong_decays.jl`. The script's inverse “implied parent mass” diagnostic uses the reference target to explain residuals; it does not feed back into the forward amplitude. `mixing_only`/unimplemented rows remain visible with unavailable results. The fitted rho/B channels are calibration evidence, not independent validation.

The realistic-factor column is **not multiplied into this baseline amplitude**. `audit_realistic_factors.jl` separately uses q-q central and contact-distorted waves, FD defaults (450,24 GeV⁻¹), and for its momentum ratio pmax=30 GeV with 2001 samples. The reported factors are wave-overlap ratios; the paper factors are comparison targets. Do not count them as hidden corrections to the baseline Table V rate.

## Table VI: photon moments and multipole amplitudes

Graph: **GI spectrum → signed physical components; P1/vector/tensor annihilation → isoscalar components; wave kernels + charges + fitted exponents + measured q → coherent moment/amplitude; footnote additions/recoil → reported quantity**.

| Input | Value / scope |
|---|---|
| Spectrum inputs | Shared values above, for each contributing flavor/sector. |
| State assignment and mixing partners | `GIPaper/data/table_vi_states.csv` n,L,multiplicity,J,flavors,mixed; only tabulated partners are solved. No implicit D/F partner is added. |
| P1 annihilation | A_np=0.50, m_eta=0.548 GeV, plus the fixed perturbative coefficient `(2pi/3)(log(2)−1)` and Gaussian alpha evaluated at unmixed model energies. Basis: 1qq,1ss,1cc,1bb,2qq,2ss,2cc. |
| Vector/tensor annihilation | A(3S1)=2.5 and A(3P2)=−0.8, four-flavor ground radial blocks. |
| Inactive alternatives | P2=0.55/M0=1.17 and the four calibrated target masses are **not used**. |
| Magnetic kernel exponent | 0.7, historically fitted to rho→pi gamma. |
| Electric/radial-moment exponent | 0.5, historically fitted to A2→pi gamma; also used in the hindered M1 recoil moment. |
| Charges | u,c=+2/3; d,s,b=−1/3. Nonstrange neutral M1 coefficient is 1 for opposite isospin or 1/3 for equal isospin; c is 4/3, s/b −2/3. Physical flavor coefficients are then composed once. |
| External masses | Registry context `VI`, parent and daughter separately; q=(Mp²−Md²)/(2Mp). Momentum-independent M1 moments can exist without assigned masses; q-dependent results cannot. |
| Unit convention | MN=0.93827 GeV for M1 moments in nuclear magnetons. alpha=1/137.036 for E1/M2 amplitudes and conversion of M1 moments into widths. alpha does not enter the tabulated M1 moment itself. |
| Footnote c | Hindered M1 term `−coefficient*q²*E2/(24m)` in the same moment convention. No extra empirical strength. |
| Footnote g | Multiply amplitude/moment by exp(−q²/(16 beta²)), beta=0.40 GeV. |
| Footnote a | Add **+0.01 nuclear magnetons** to phi→pi gamma for the paper's supplied pi0–eta mixing contribution. This is a channel-specific physical input in GIPaper, not an emergent GIModel prediction. |

Default computation uses adaptive HO with three retained levels/channel, restricted to the explicitly selected states. `GI_TABLE_VI_SOLVER=fd` selects the separate 700-point/rmax=24 comparator. The input record stores the actual solver and wave components. P1 momentum smearing uses its declared 900-point scheme on FD; HO uses native momentum integrals. The generic FD photon transform uses pmax=30 GeV, 1501 points. These are distinct from Table VII's FD observable adapter.

E1/M2 amplitudes are in MeV¹ᐟ² and their squares are widths in MeV. M1 entries are moments; their width conversion requires alpha, q³ and the initial-spin factor. The nucleon-mass convention cancels when a consistently normalized moment is converted back to a width.

## Table VII: four different observable families

| Subtable | Computation and additional input values |
|---|---|
| (a) leptonic factors | Charge/flavor coefficient × P_P, V_V, V′_V or P′_A1. Constituent masses, wave, experimental M and derived mock mass Mtilde enter. The m/E power is the formula's fixed unity. **Neither alphaEM nor GF enters the printed decay constant.** S states use contact-distorted waves; D-vector and P-axial rows use central waves. Exact coefficients (2√3, √6, etc.), n and flavors are in the row records. The implementation here uses the declared flavor rows rather than the full Table-VI four-flavor mixing construction. |
| (b) two photons | alpha=1/137.036, constituent mass m, measured M, computed Mtilde, `(M/Mtilde)^(3/2)`, charge-weighted linear momentum integral. Pseudoscalars use contact-distorted S waves. The A2/f/f′ rows use central P waves and ideal flavor assignments. q_eff: pi-like 1/(3√2); qq isoscalar 5/(9√2); charm 4/9; strange/bottom 1/9. |
| (b), mixed eta family | P1 inputs 0.50 and 0.548 GeV and calculated mixing in the **four-state** 1qq,1ss,2qq,2ss basis. Each component uses its constituent mass and charge; the external physical mass is the registry mass of the final state. This differs deliberately from Table VI's seven-state P1 basis. No calibration-to-four-target-masses is used. |
| (c) gluons | Calculated S_L, constituent mass and alpha_s(M), where **M is the assigned experimental mass**, not the model eigenvalue. Lowest-order channel factors: S0 8pi alpha_s²/(3m²); S1 40(pi²−9)alpha_s³/(81m²); P2 32pi alpha_s²/(45m²); P0 8pi alpha_s²/(3m²). All use fixed spin-sector waves. No extra fitted transition strength. |
| (d) charge radii | Constituents and charges (pi+: +2/3,+1/3; K+: +2/3,+1/3; K0: −1/3,+1/3), distorted S wave, fitted position-smearing exponent **0.2**, and `(hbar c)²=0.19733²` for fm². No experimental meson mass enters the predicted radius. pi+ is the historical fit anchor. |

HO starts at 24 basis states and uses the solver's adaptive defaults. Ordinary families request at least **four radial levels**, which also affects the variational beta selection. Mixed pseudoscalars are constructed through their separate two-radial-state spectrum request. `NPTS=900` is passed to observable adapters; on native HO it does not create a 900-point momentum mesh. Constants `NGRID=1200` and `RMAX=24` remain in the script but are **unused by this HO calculation**.

There are 61 canonical Table VII rows. The implementation defines 57 rows across these families; four hypothetical top rows have no constituent top mass and are not computed. Within the defined rows, missing experimental assignments also leave mass-dependent predictions unavailable. An array entry is not proof of a finite prediction.

The extra D7/D8 validation/width section is separate from the printed table:

- GF=1.1663787e−5 GeV⁻² enters the pion weak-width demonstration, together with the registry muon and pion masses. The measured fpi=0.1307 GeV is a **validation input**, not a replacement for the model's predicted fpi.
- hbar=6.582119e−25 GeV s and lifetime=2.6033e−8 s give its comparison width.
- The 5.55 keV psi dilepton width is inverted to an experimental fpsi for a formula check.
- Predicted rho/psi/psi′/Upsilon dilepton widths instead consume their **model** decay constants, alpha and registry masses; measured 7.04/5.55/2.34/1.34 keV values only form comparison ratios.
- CKM/Cabibbo factors and QCD radiative corrections are not applied. The D9 tau-width API exists, but this harness does not evaluate a D9 width merely because the table contains an axial decay constant.

## What the tracing changes establish

There is no universal “rates use the same handful of numbers” claim. The baseline Table V rate is a calibrated analytic model; Table VI adds fitted electromagnetic kernels and explicit paper footnotes to calculated states; Table VII mostly uses wavefunction operators but includes a fitted charge-radius exponent and explicit wave/mixing conventions.

The cleanup makes numerical ownership visible without changing these prescriptions. Moving operators does not validate the external mixing angles, the +0.01 input, ideal-flavor assumptions, coefficient transcription precision or the use of current kinematics against 1985 targets. Those are now identifiable review items for the next folder-by-folder round.

## Verification of this migration

GIModel: 1,011 checks passed; QuarkModelTransitions: 551; GIPaper: 2,638.
The reference manifest resolves all 103 units with zero broken links. Regenerated
Table V and VII reports and the Table VI result CSV have no numerical changes;
only the Table VI report generation date changed. The core suite still emits
its existing coarse-grid and quadrature-exhaustion warnings; this migration
does not claim to resolve those numerical-limit cases.

## Policy records and freshness checks

`data/table_policy.toml` owns the paper's calibration amplitudes, rounded mixing
angles, charmed-row form-factor assignment, Table VI footnote-a addition and
isovector assignment, and comparison-only legacy validation anchors. Generic
operator defaults remain in QuarkModelTransitions. The scripts load this record;
they do not fit those values against the output residuals. The edition of the
legacy validation widths was not recorded; the data explicitly says so.

Run `julia --project=GIPaper/scripts GIPaper/scripts/trace_rate_inputs.jl` to regenerate the three reports
and traces together, then `julia --project=GIPaper/scripts GIPaper/scripts/check_input_traces.jl` for a fast,
read-only source-fingerprint and canonical-row coverage check. Generation aborts
if its tracked inputs change during computation. The full project gate runs both.
The checker detects stale sources; it is not an independent reimplementation of
the rates. Numerical quadrature literals in the common trace are an audited
snapshot of implementation defaults, so changes to those kernels require reviewing
the corresponding trace declarations as well as regenerating results.
