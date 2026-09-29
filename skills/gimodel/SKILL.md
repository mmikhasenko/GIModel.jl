---
name: gimodel
description: Compute meson masses, wavefunctions, mixing and decays with GIModel.jl, the Julia implementation of the Godfrey–Isgur (1985) relativized quark model. Use when the user asks for quark-model meson spectra (charmonium, bottomonium, D, D_s, B, B_c, K, light and isoscalar mesons), singlet–triplet or η–η′ mixing, radial wavefunctions, radiative (E1/M1), strong (A→Bπ), leptonic, two-photon or gluonic widths, or compares results with Godfrey–Isgur, Phys. Rev. D 32, 189 (1985).
---

# GIModel.jl

Julia ≥ 1.11. Not registered yet:

```julia
using Pkg; Pkg.add(url = "https://github.com/mmikhasenko/GIModel.jl")
using GIModel                          # spectra, states, wavefunctions
using GIModel.QuarkModelTransitions    # decays; a submodule, nothing extra to install
```

Docs: https://mmikhasenko.github.io/GIModel.jl/dev/ (manual → conventions,
spectra, solvers, isoscalar, transitions; tutorials; "The 1985 paper").
Every exported function has a docstring. Read `?name` before guessing arguments.

## Core recipe

```julia
params, mq = load_parameters_and_quark_masses(default_parameters_path())   # Table II
spec = compute_spectrum(params, Meson(mq, :c, :c);
    levels = spectrum_levels(2),            # all n ≤ 2 levels in S, P, D
    solver = OscillatorSolver())            # the paper's method, certified energies
spec                                        # prints a table with the contributions
s = spectrum_state(spec, "1^3P_1")          # or spectrum_state(spec, 1, "P", 3, 1)
s.mass_GeV, s.central_GeV, s.contact_shift_GeV, s.fine_structure_shift_GeV
physical_components(spec, "1^3S_1")         # signed (basis, coefficient, wave) list
```

To get a quick table without writing code, run
`julia --project=<env> scripts/spectrum.jl c c 2 [--L SPD] [--solver ho|fd]`
from this skill's directory. It prints CSV plus a numerics-provenance line.

Smoke values (HO, default parameters, `spectrum_levels(2; L_labels=("S","P"))`):
η_c 2.9671, J/ψ 3.0914, χ_c2 3.5482 GeV; the `u d̄` 1¹S₀ is 0.149 GeV. If
these come out wrong, the setup is wrong.

## Rules that prevent silent errors

**Labels and ordering**
- Levels are `"n^(2S+1)L_J"`, e.g. `"1^3P_2"`. `n` starts at 1 in *every*
  partial wave, so the lowest P wave is `1P`, not `2P` (the spectroscopic
  convention used by the PDG and older papers).
- `Meson(mq, f1, f2)` is f1 (quark) and f2-bar (antiquark). The masses are
  symmetric, but the singlet–triplet mixing angle flips sign when you swap the
  quark order. `Meson(mq, :c, :u)` is a D⁰.
- Flavors: `:u :d :s :c :b`, plus `:q` for the averaged light quark. There is no
  top quark and no isospin breaking.
- After mixing, a state keeps the label of the unmixed level it is assigned to,
  by ascending mass. The label identifies the state, not its composition.

**Mixing depends on the levels you request**
- A mixing block contains only the levels you asked for. `spectrum_levels(1;
  L_labels=("S","P"))` has no D waves, so the J/ψ shows up as *unmixed*. Always
  request every level that can mix with the state you study (same J^PC, and
  n ± 1). The one exception is radiative transitions (see Transitions).
- Unequal-mass mesons mix ¹L_L and ³L_L (for example D₁ and D₁′).
  Self-conjugate mesons do not. The lower state is
  `cos θ |¹L_L⟩ + sin θ |³L_L⟩`. Read the angle from
  `physical_components`; don't compute it from masses.
- `radial_wave(spec, label)` *throws* for a mixed state. That is intended. Use
  `physical_components` or `physical_state`. Never pick one component's wave as
  a stand-in for the physical state.

**Solvers: numerics versus physics**
- The default is `FiniteDifferenceSolver()`: fast and accurate to <1 MeV for
  most states, but it has no convergence check. Use `OscillatorSolver()` for
  numbers to quote and for any comparison with the paper.
- HO certifies **energies only**. Widths, values at the origin, radii and
  cancellation-prone overlaps have to be checked separately: tighten
  `energy_tolerance_GeV` or raise `nbasis`, or refine the FD `ngrid`/`rmax`,
  and compare. Gluonic, two-photon and leptonic widths are the most sensitive.
- If HO and FD disagree beyond about 1 MeV, suspect convergence, not physics.
- Solver settings must not change results. `SpinTerms(contact_hyperfine=false,
  fine_structure=..., same_j_spin_orbit=..., tensor=...)` changes the physics on
  purpose, to isolate a term.
- Report `numerics_provenance(solver)` next to any quoted number.

**Parameters**
- Parameter objects are immutable. To vary one, copy it:
  `GIParameters(params; potential = ConfinementPotential(params.potential; b = 0.20))`.
- Quark masses are not in `GIParameters`. For a continuous mass, use
  `Meson(:c, :c, ConstituentMasses(m, m))` or
  `Meson(HeavyQuark{:up}(m, :Q), LightQuark(mq["q"]))`.
- Units: GeV for masses and momenta, GeV⁻¹ for r (×0.1973 → fm), and **MeV for
  widths** from `decay_width`. The TOML file stores masses in MeV; the loader
  converts them to GeV.

**Isoscalars (η, η′, ω, φ, f₂, f₂′)**
- `compute_spectrum` never applies annihilation mixing. Use
  `compute_isoscalar_spectrum(params, Meson(mq,:q,:q), Meson(mq,:s,:s); levels,
  pseudoscalar = PaperP1Annihilation() | PaperP2Annihilation(),
  amplitudes = Dict(("S",3,1) => params.annihilation.s1_A, ...))`. Channels
  you don't list stay ideally mixed.
- For pseudoscalars, the levels must include both `1^1S_0` and `2^1S_0`: the
  paper uses a four-state basis, and the η′ moves by about 300 MeV without it.
  P2 matches Table III compositions better than P1.
- Labels repeat across channels, so select with
  `BasisState(1, "S", 1, 0; flavors = (:s, :s))`.
- To add c c̄ or b b̄, solve each channel with the same params and solver, then
  call `add_isoscalar_annihilation`.

## Transitions

Always follow the same pattern: `physical_state` → operator →
`matrix_element(final, op, initial)` (note the order: ⟨f|O|i⟩) →
`decay_width` (MeV).

```julia
spec = compute_spectrum(params, Meson(mq, :c, :c);
    levels = spectrum_levels(2; L_labels = ("S", "P")), solver = OscillatorSolver())
psi, eta_c = physical_state(spec, "1^3S_1"), physical_state(spec, "1^1S_0")
decay_width(matrix_element(eta_c, PhotonEmission(mq), psi))            # M1, 2.4 keV
decay_width(MasslessLeptonPair(), LeptonicCurrent(:electromagnetic, mq), psi)
decay_width(TwoPhotonChannel(), TwoPhotonAnnihilation(mq, AnnihilationTerm((:c,:c), 4/9)), eta_c)
decay_width(TwoGluonChannel(), GluonicAnnihilation(mq, M -> alpha_s_q(M)), eta_c)
```

- `PhotonEmission` picks the multipole itself (`verbose = true` shows which).
  It **rejects `:q`**: build states from `:u`/`:d`, or take them from an
  isoscalar spectrum, which resolves u and d on its own.
- **`PhotonEmission` versus D waves.** It rejects any state with a nonzero
  D-wave admixture ("Not implemented yet"). With `spectrum_levels(2)`, the J/ψ
  picks up a ³D₁ amplitude of about −0.013 through the tensor force, so J/ψ → η_c γ fails. The
  docs' radiative examples request only S and P waves for this reason. That is
  an approximation: the ³S₁–³D₁ mixing is neglected. Say so when you use it,
  and never drop components from a mixed state to get past the error.
  Leptonic, two-photon and gluonic operators accept mixed states. For example,
  Γ(J/ψ → ee) changes from 10.51 to 10.46 keV when D waves are included.
- `PseudoscalarEmission(g, h, mq)`: the couplings `g` and `h` (GeV⁻¹) are
  **inputs, not predictions**, and have to be calibrated against known widths
  (see the strong-decays tutorial). The channel is
  `TwoMesonChannel(surviving, emitted)`, and the emitted meson must be 0⁻.
  Build neutral states such as π⁰ with
  `superpose([uu, dd], [1, -1]; label, mass_GeV)`.
- Kinematics use **model masses**. To evaluate at a measured mass, multiply the
  *amplitude* by `mass_correction_factor(...; target_mass | target_momentum)`.
- Not implemented: ³P₀ pair creation, emission of non-pseudoscalar mesons, weak
  hadronic or semileptonic decays, coupled channels and threshold effects
  (X(3872)-like states), baryons, exotics. Say so; don't approximate.

## Comparing with the 1985 paper

Read [references/paper-comparison.md](references/paper-comparison.md) before
claiming agreement or disagreement with Godfrey–Isgur (1985). In short:
- Masses agree to a few MeV on average across all seven sectors.
- **The mixing angles in the captions of Figs. 4, 7 and 9 are not reproduced
  on purpose.** They appear to be a paper erratum; compare with
  Godfrey–Kokoski 1991 instead. Never tune the model toward the caption angles.
- Table V strong decays use a frozen single-Gaussian (β = 0.40 GeV)
  reproduction. The native `PseudoscalarEmission` path gives different numbers,
  and that is expected.

## Before reporting a result

1. Were all levels that can mix requested?
2. Is the solver HO (or FD with a refinement study), and did you check
   convergence of the specific observable if it is a width or overlap?
3. Are units stated (GeV masses, MeV widths)? Is the quark order stated for
   mixing angles?
4. Does it rely on model masses for kinematics? Say so.
