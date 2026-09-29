# Comparing GIModel with Godfrey–Isgur (1985)

S. Godfrey and N. Isgur, Phys. Rev. D 32, 189 (1985), doi:10.1103/PhysRevD.32.189.
The model uses the paper's Table II parameters and literal Appendix A operators.
**Nothing is tuned toward the paper**; differences are reported, not fitted.
Full account: https://mmikhasenko.github.io/GIModel.jl/dev/ → "The 1985 paper".

## Setup for any comparison

- Use `OscillatorSolver()` explicitly (the paper's method). The package default
  is finite differences.
- The comparison layer is the separate `GIPaper` package in the repository
  (`GIPaper/`), which depends on GIModel:
  ```julia
  using GIModel, GIPaper
  params, mq = load_parameters_and_quark_masses(GIPaper.model_parameters_path())
  ref = load_reference_spectrum(reference_spectrum_path("charmonium"))
  rows = compare_reference(params, mq, ref; solver = OscillatorSolver())
  ```
  Paper reference values live only in GIPaper; GIModel never sees them.
- GIPaper decay kinematics use pinned PDG 2026 masses with no fallback to model
  masses. `reference_meson` has no fallback masses either: don't invent one.

## What agrees, and how well

| paper | GIModel | agreement |
|---|---|---|
| Figs. 4–9 spectra (209 states, 7 sectors) | `compute_spectrum` | mean \|Δ\| 2.6–6.1 MeV per sector, max 24.8 MeV (charmed) |
| Table III isoscalar mixing | `compute_isoscalar_spectrum` | P2: component RMS 0.015, masses ~11 MeV; P1: η good, η′ and radials worse (RMS 0.08, ~56 MeV); vector/tensor blocks to printed precision |
| Tables IV–V strong decays | frozen single-Gaussian backend | all amplitudes, most within a few % (refit A = 1.644, S₀ = 3.291 vs paper 1.67, 3.27) |
| Table VI radiative | `PhotonEmission` on solved waves | median model/paper 1.00 (M1), 1.01 (E1); sign 58/62 |
| Table VII gg, γγ, e⁺e⁻, f, radii | annihilation operators | medians 1.02–1.06; mixed-isoscalar γγ 0.78 |

The paper's spectra are read off printed level diagrams, so they carry a few
MeV of digitization uncertainty, which is comparable to the mean residuals.

## Known discrepancies: don't "fix" these

1. **Same-J mixing angles in the captions of Figs. 4, 7 and 9.** GIModel
   reproduces the masses but not most caption angles (u s̄: +4.3° vs +34°;
   b c̄: +69.2° vs −53°). The same authors' later papers with the same
   parameters agree with GIModel within a few degrees: Godfrey–Kokoski,
   PRD 43, 1679 (1991); Godfrey, PRD 70, 054017 (2004); Godfrey–Moats,
   PRD 93, 034035 (2016). The captions appear to be an erratum. Compare
   against both, and state which one you used.
2. **Transitions from radially excited η states** (Table VI). These involve
   strong cancellations and depend on the isoscalar basis, so large ratios are
   expected.
3. **Table V is a simplification.** The paper used one Gaussian with
   β = 0.40 GeV for every meson, and the leading constant S₀ instead of the
   full q-dependent polynomial of Table IV. The native `PseudoscalarEmission`
   (solved waves, nodes, mixing) goes beyond this and gives different numbers
   by design. The two paths are never mixed.
4. **³S₁–³D₁ admixtures near ψ(3770)**: still open.
5. **Toponium rows are not computed.** Table II has no top mass.

## Common mistakes seen in comparisons

- Comparing with FD results and blaming the physics for sub-MeV differences.
- Requesting too few levels, which changes the mixing blocks, then comparing
  compositions.
- Using the other quark order (c ū vs u c̄), which flips the mixing angle.
- Reading the angle from masses instead of from `physical_components`.
- Comparing η′ with a two-state pseudoscalar basis; the paper uses four states.
