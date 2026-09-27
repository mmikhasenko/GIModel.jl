# Testing and documentation issues

Issues found while writing the documentation (September 2026). Each one needs
an investigation before a fix: the behavior may have a reason. Every entry
gives a reproduction, what is known, and questions to answer first.

All snippets assume:

```julia
using GIModel, GIModel.QuarkModelTransitions
params, mq = load_parameters_and_quark_masses(default_parameters_path())
```

---

## 1. The two solvers converge to different mixing angles

**Resolution (September 2026).** The oscillator path applied non-polynomial
momentum factors to the finite `p²` matrix. In general,
`f(P_N p² P_N) != P_N f(p²) P_N`; optimizing β during energy refinement could
hide this operator error. It affected both the Hamiltonian and the mixing
element, especially the cancellation of vector and Thomas terms.

The Hamiltonian now projects continuum momentum functions by independently
refined quadrature and retains an enlarged intermediate basis in sandwiches.
Wave-interface spin sandwiches separately refine their auxiliary operator
basis with the input wave held fixed. This includes diagonal expectations and
cross elements with different β or L. Failure to converge raises an error.

An independent momentum-space integral for analytic oscillator waves checks
the operator; mesh comparisons and zero-padding invariance provide additional
regressions. `test/heavy/mixing_convergence.jl` compares all three systems,
including a fixed common β, with a 0.1° tolerance for residual wave/grid errors.
A three-system check gives:

| system | FD 900 | corrected HO |
|---|---:|---:|
| u s̄ | 4.30° | 4.31° |
| c ū | −25.78° | −25.79° |
| b ū | −27.96° | −27.97° |

Published angles are not used as numerical targets: the regression compares
two independent implementations of the same Hamiltonian. Remaining differences
are governed by wave and grid convergence. The original observations and
investigation questions follow for context.

**Symptom.** The singlet–triplet mixing angle of unequal-mass P waves differs
between the finite-difference (FD) and oscillator (HO) solvers by 1–2°, and
each value is stable under refinement of its own solver.

| system | FD 450 | FD 1800 | FD 2400/32 | HO | HO, ΔE ≤ 1 keV |
|---|---:|---:|---:|---:|---:|
| u s̄ | 4.30° | 4.30° | 4.31° | 2.34° | 2.34° |
| c ū | −25.77° | −25.78° | −25.79° | −25.04° | −25.04° |
| b ū | −27.93° | −27.96° | −27.96° | −26.82° | −26.82° |

The masses agree to well below 1 MeV.

**Reproduction.**

```julia
levels = spectrum_levels(1; L_labels = ("P",))
for solver in (FiniteDifferenceSolver(ngrid = 1800), OscillatorSolver())
    spec = compute_spectrum(params, Meson(mq, :u, :s); levels, solver)
    println(spectrum_state(spec, "1^1P_1").mixings[1].components)
end
```

**Known.** The off-diagonal element is a small difference of two larger
terms (for u s̄: vector +7.8 MeV, Thomas −8.3 MeV, total −0.6 MeV), and the
diagonal gap is 13.5 MeV, so a 1% difference in either piece moves the angle
visibly. The GIPaper reports use FD.

**Questions.**
- Does the off-diagonal element itself differ between solvers? Compare
  `spin_orbit_mixing_components` for the same pair of waves, per piece.
- In HO, the ¹P₁ and ³P₁ sectors each choose their own β. The cross matrix
  element then pairs oscillator waves with different β. Is that path
  (quadrature, momentum sandwich across two bases) as accurate as the
  same-β one? Test by forcing a common β (`OscillatorSolver(beta_grid = [β])`).
- Is the FD cross-sandwich `radial_cross_expect_momentum_sandwich` exactly the
  same operator as the HO one? Check it on analytic waves.
- Which answer does an independent calculation support (e.g. the
  Godfrey–Kokoski 1991 angles: +5°, −26°, −31°)?

---

## 2. Isoscalar states cannot enter radiative or leptonic decays

**Resolution.** Isoscalar spectra now retain explicit nonstrange-isoscalar
provenance. Standard photon and electromagnetic leptonic currents resolve that
channel coherently into uū and dd̄ with coefficients 1/√2. Ambiguous ordinary
`:q` spectra still fail. Explicit effective-charge annihilation operators keep
the original coarse components, avoiding double counting. Regression tests
cover mixed ω/φ currents, φ radiative decay, and mass corrections. The Table VI
audit now uses public photon currents with explicit isospin states; its charged
A1/A2 spin-flip rows use u d̄ rather than an internal emitter coefficient.

**Symptom.** States from `compute_isoscalar_spectrum` carry the nonstrange
component as flavor `:q`. `PhotonEmission` and `LeptonicCurrent` reject `:q`,
so ω → e⁺e⁻, φ → e⁺e⁻ (through its small n n̄ admixture) and φ → ηγ all fail.

**Reproduction.**

```julia
iso = compute_isoscalar_spectrum(params, Meson(mq, :q, :q), Meson(mq, :s, :s);
    levels = spectrum_levels(2; L_labels = ("S",)),
    pseudoscalar = PaperP1Annihilation(), amplitudes = Dict(("S", 3, 1) => 2.5))
phi = physical_state(iso, BasisState(1, "S", 3, 1; flavors = (:s, :s)))
decay_width(MasslessLeptonPair(), LeptonicCurrent(:electromagnetic, mq), phi)
# ArgumentError: specify explicit u, d, s, c, b or t flavor, not q
```

**Known.** Rejecting `:q` is deliberate: in general `:q` cannot tell
(uū + dd̄)/√2 from (uū − dd̄)/√2. In an isoscalar spectrum, however, the
nonstrange channel is always the isoscalar combination (the code already marks
it `isoscalar_coherent`). GIPaper's Table VI audit works around this with the
internal `PhotonEmitter` currents. Two-photon widths work because
`AnnihilationTerm` takes the effective charge explicitly.

**Questions.**
- Should the isoscalar channel carry its own flavor tag (for example `:n` for
  (uū + dd̄)/√2) that the currents accept?
- Or should `physical_state` expand an isoscalar `:q` component into explicit
  `:u` and `:d` components with 1/√2 coefficients?
- Whichever is chosen, can GIPaper then drop its internal workaround and
  reproduce Table VI through the public API?

---

## 3. Photon emission rejects states with a D-wave component

**Resolution (requested policy).** Unsupported nonzero D-wave components
raise `ErrorException("Not implemented yet. Please submit issue if needed,
and/or PR with implementation.")`. No components are silently discarded and
no projection helper is introduced. General D→P kernels remain future physics
work, requiring a derivation and validation against the implemented S/P cases.

**Symptom.** When D waves are requested, tensor mixing gives the J/ψ a small
³D₁ component (amplitude ≈ 0.01). `PhotonEmission` then refuses the state.

**Reproduction.**

```julia
spec = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(2),
    solver = OscillatorSolver())
matrix_element(physical_state(spec, "1^3S_1"), PhotonEmission(mq),
    physical_state(spec, "1^3P_1"))
# ArgumentError: photon kernel supports S--S spin flips and S--P transitions; got L=1 -> L=2
```

**Known.** The paper's Appendix D operators cover S–S (M1) and S–P (E1, M2),
and Table VI uses states without D admixture. The same restriction makes
radiative decays of D-wave states, such as ψ(3770) → χ_cJ γ (a D → P E1
transition), impossible.

**Questions.**
- Is the restriction a limitation of the Appendix D formulas, or only of the
  implemented angular kernels?
- Would a general E1 kernel (any L → L ± 1) and M1 kernel (same L) cover the
  missing cases, and can they be checked against the S–P results?
- Until then, should small components be allowed to drop with a warning, or
  must the user choose the basis explicitly (current behavior)?

---

## 4. Different default solvers

**Resolution.** Both general spectrum entry points now default to FD.
GIPaper comparisons instead default to HO, and their production spectrum and
mixing-angle audits select `OscillatorSolver()` explicitly, following the
original investigation. FD remains available as an explicit cross-check. The
general FD default preserves
`compute_spectrum` behavior and supports diagnostic central-potential variants
that native HO deliberately rejects. The paper's oscillator method remains an
explicit `solver = OscillatorSolver()` choice. Combining already-computed
spectra still retains their solver, with incompatible solvers rejected. This
preserves existing single-meson defaults; GIPaper's distinct HO default is
documented explicitly, and historical FD reports retain their solver labels.
Warm timings on the development machine (one BLAS thread, `spectrum_levels(2)`):
u–s, FD 5.13 s / HO 13.99 s; c–c, FD 2.84 s / HO 6.20 s.
These timings are illustrative, not performance guarantees.

**Symptom.** `compute_spectrum` defaults to `FiniteDifferenceSolver()`, while
`compute_isoscalar_spectrum` defaults to `OscillatorSolver()`. A user who
combines `compute_spectrum` results in `add_isoscalar_annihilation` gets FD,
and one who calls `compute_isoscalar_spectrum` gets HO.

**Known.** HO is the paper's method and certifies its energies; FD with
default settings is faster in some cases and was the historical default for
the GIPaper reports. Results differ by less than 1 MeV in masses, but see
issue 1 for angles.

**Questions.**
- What are the run times of the two defaults for typical requests (one
  meson, `spectrum_levels(2)`)?
- Would switching `compute_spectrum` to HO change any recorded GIPaper report
  beyond its tolerance?
- Is there a reason to keep FD as the default (speed, robustness for light
  or high-L states)?

---

## 5. The keyword constructor of `GIParameters` builds a different model

**Resolution.** Construction from scratch requires all six keyword fields.
The loader and copy constructor already supply them; repository tests and
research probes use those paths and need no model changes. Diagnostic
`RelativisticFactors()` defaults remain available when deliberately requested.
Missing physics fields now fail immediately instead of silently creating a
pointwise, unrelativized model.

**Symptom.** `GIParameters(; potential, smearing)` fills the other fields with
`PointwiseCentral()` and `RelativisticFactors()` with all ε = 0 and all
momentum sandwiches off. That is not the Godfrey–Isgur model, and nothing
warns the user.

**Reproduction.**

```julia
p = GIParameters(; potential = params.potential, smearing = params.smearing)
p.central, p.factors
```

**Known.** The docs recommend the copy constructor
`GIParameters(params; potential = ...)`, which keeps everything else. The
simplified defaults are used by diagnostic tests.

**Questions.**
- Which tests and scripts rely on the simplified defaults?
- Should the keyword defaults be the paper's model, with the diagnostic
  variants requested explicitly?
- Or should the keyword constructor require every field?

---

## 6. Signs of spectroscopic mixing eigenvectors are arbitrary

**Resolution.** All mixing blocks now use the annihilation rule by default:
positive overlap with the assigned precursor, in ascending unmixed-mass order.
A zero overlap uses the largest component. The common helper is shared with
annihilation; `phase_anchor=1` remains an explicit compatibility option. Only
overall state phases change, not eigenvalues, relative internal signs or widths.

**Symptom.** Annihilation eigenvectors are phase-fixed (positive overlap with
the assigned basis state), but spin-orbit and tensor eigenvectors keep the sign
returned by the diagonalization.

**Reproduction.**

```julia
spec = compute_spectrum(params, Meson(mq, :u, :s); levels = spectrum_levels(1; L_labels = ("P",)))
[(c.basis.label, c.coefficient) for c in physical_components(spec, "1^3P_1")]
# the assigned ³P₁ component is ≈ −0.997
```

**Known.** Widths do not depend on an overall state phase, and relative signs
inside a state are correct. The sign still matters when a mixed state is
composed further (for example, spin mixing followed by annihilation) or when
amplitudes are compared with published signs. A similar problem with
oscillator-basis column signs was found and fixed earlier.

**Questions.**
- Can the annihilation phase rule be applied in `diagonalize_mixing_block`
  for every mechanism?
- Does any report or test depend on the current signs?

---

## 7. Missing conveniences

**Resolution.** `convergence(spec, state)` returns the source-channel energy
certificates for every component (FD has no automatic certificate). Native
oscillator momentum waves are callable as `phi(p)`, with no mesh transform.
`superpose(states, coefficients; mass_GeV, label)` constructs coherent,
optionally normalized combinations, including interference in the norm.
The mass must be supplied rather than inferred from a superposition. The
`AnnihilationTerm` documentation now explains the color factor √3 and the
Appendix-D factor 2 in the pseudoscalar coefficient, and distinguishes CKM,
state mixing, electric charge, and current normalization.

Each is small, but each forces users into internals:

- **Convergence certificate.** It is reachable only through
  `spec.computation.channel_cache`. A public accessor (for example
  `convergence(spec, state)`) would let users check it directly.
- **Momentum waves.** An `OscillatorMomentumWave` cannot be evaluated at a
  momentum through the public API; plotting goes through `sample_wave` and a
  numerical transform.
- **Isospin states.** Building π⁰ or ω from `:u` and `:d` spectra takes a
  hand-written component comprehension (see the strong-decay tutorial). A
  helper for linear combinations of `PhysicalState`s would remove it.
- **Leptonic current coefficients.** The `AnnihilationTerm` coefficient for a
  weak pseudoscalar decay is `2√3` in the docstring example; its meaning
  (color and normalization factors) is not explained, so users cannot
  construct other cases confidently.

**Question.** Which of these belong in the public API, and under what names?
