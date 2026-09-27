# Open issues

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
