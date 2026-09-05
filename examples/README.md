# Examples

Worked examples of **GIModel** used as a physics tool rather than as a paper
reproduction. Each one is self-contained; the curated examples run against the
model's public API only — no reaching into internals, no reference data.

The curated examples share one environment (`examples/Project.toml`): GIModel
(dev-path to the repository root) plus CairoMakie, PlutoUI and LaTeXStrings.
The exception is `played_with_model.jl`, the scratch notebook, which wants the
comparison layer and so activates `GIPaper/` instead. Either way each file
activates its own environment, so there is no setup step.

One caveat on "public API only": `chi_c_annihilation_widths.jl` calls
`GIModel.radial_grid`, which is not exported even though the exported
`fine_structure_grid_operator` needs a mesh from it. That is an API gap, not an
example reaching past the front door.

## `adaptive_ho_refinement.jl` — how an HO answer earns its certificate

A Pluto notebook that follows one charmonium `¹S₀` channel from an initial
24-function harmonic-oscillator basis to a certified result. It explains the
common-`β` variational search, replays each visited basis size as a convergence
table and plot, identifies every loud failure condition, and finishes with a
five-question interactive quiz.

```bash
julia -e 'using Pluto; Pluto.run(notebook="examples/adaptive_ho_refinement.jl")'
```

The production call returns the ordinary `ChannelRadialSolution` used by later
physics stages plus its compact `OscillatorConvergence` certificate. The visible
history is deliberately reconstructed with public, fixed-size
`OscillatorSolver(converge=false)` calls: those probes label themselves
`:unchecked` and are never confused with the certified result. No spatial mesh
or finite-difference fallback appears anywhere in the notebook.

The notebook also states the certificate's boundary explicitly: matrix-element
quadrature has its own convergence warning, which remains visible and is not
misrepresented as covered by basis refinement.

The notebook can also be executed non-interactively as a smoke test:

```bash
julia --project=examples examples/adaptive_ho_refinement.jl
```

## `heavy_quark_transition.jl` — charmonium → bottomonium

A Pluto notebook that turns the heavy-quark mass into a slider and watches
charmonium walk over into bottomonium.

```bash
julia -e 'using Pluto; Pluto.run(notebook="examples/heavy_quark_transition.jl")'
```

The premise is a property of the model: the GI medium (`b`, `c`, `σ₀`, `s`,
`α_s(Q²)`) is flavour-blind, so `m_Q` is a free external dial — nothing has to
be retuned when it moves. See the [`HeavyQuark`](../src/quark.jl) docstring.

**Controls**: `m_Q` from 1.2 to 5.4 GeV; a flavour switch between `Q Q̄`
(quarkonium — both masses move), `Q q̄` and `Q s̄` (heavy-light — only the heavy
quark moves); an axis lock; and a toggle to subtract `m₁ + m₂` and look at the
binding part alone.

**Two views**, both live on the slider:

| view | what it shows |
|---|---|
| level scheme | one column per `ⁿ²ˢ⁺¹L_J`, plus each level's trajectory in `m_Q` and the collapsing `1S` hyperfine splitting |
| `J^P` shape | each `L`-block collapsed onto the quantum numbers that survive, with the frame rescaled to the current spectrum so only the *pattern* moves |

**What the model says**, visible on the way from `m_c` to `m_b`: levels climb by
`≈ 2Δm_Q` (rest mass); the `1S` hyperfine splitting falls from ~135 MeV to
~75 MeV, the `1/m_Q` of the spin-spin contact term; the excitation pattern is
nearly invariant once the quark masses are subtracted, because the confining
string is the same in both systems; and for heavy-light the binding curves
flatten onto `Λ̄ = M − m_Q`, which is heavy-quark symmetry appearing on screen.

Recorded sweeps are in `figures/` (`charm_to_bottom_*.gif`); the notebook's last
section re-renders them, one model solve per frame.

### The API this exercises

The notebook is a fair sample of how the package is meant to be driven:

```julia
params, quark_masses = load_parameters_and_quark_masses(path)   # parameters.provisional.toml

Q = HeavyQuark{:up}(mQ, :Q)                # the dial: mass is a free parameter
meson = Meson(Q, LightQuark(quark_masses["q"]))   # or Meson(Q, Q) for quarkonium

spec = compute_spectrum(params, meson;
    levels = spectrum_levels(2; L_labels = ("S", "P", "D")),
    solver = RadialSolver(; ngrid = 360, rmax = 20.0, nlevels_per_channel = 3))

spectrum_state(spec, "1^3S_1").mass_GeV    # one level by its own label
```

Points worth copying:

* **Quarks are objects, mesons are built from them.** `Meson(q1, q2)` takes the
  masses and the flavour symbols; equal symbols mark a self-conjugate channel,
  which is what switches the same-`J` antisymmetric spin-orbit mixing off for
  `Q Q̄` and on for `Q q̄` — no flag needed at the call site.
* **`RadialSolver` is numerics, `SpinTerms` is physics.** Changing the solver
  must not move a mass; the notebook uses a coarser mesh than the default
  (`ngrid = 360, rmax = 20.0` against `450 / 24.0`) purely for speed, having
  checked the masses agree to well under an MeV across the whole mass range.
  Changing a `SpinTerms` switch is *meant* to move a mass.
* **The spectrum is staged.** `compute_spectrum` composes the complete
  fixed-sector solve and intra-meson mixing; call `fixed_spectrum` or
  `add_intra_meson_mixing` directly to inspect a production stage.
  `central_spectrum` is an independent central-only diagnostic. Read the per-state breakdown
  (`central_GeV`, `contact_shift_GeV`, `fine_structure_shift_GeV`, `mixings`)
  from the states themselves.
* **A `Spectrum` prints its own physics.** Displaying one gives a table of
  levels and contributions, which is why the notebook just shows `spec`.

## `chi_c_annihilation_widths.jl` — why χ_c0 → gg is six times χ_c2 → gg

A script that computes the two-gluon and two-photon widths of the two charmonium
P-wave states and then takes the ratio apart.

```bash
julia examples/chi_c_annihilation_widths.jl
```

| | M (GeV) | Γ(gg) | Γ(γγ) |
|---|---:|---:|---:|
| χ_c0 (³P₀) | 3.458 | 5.60 MeV | 3.25 keV |
| χ_c2 (³P₂) | 3.533 | 0.910 MeV | 0.623 keV |
| ratio 0/2 | | 6.15 | 5.21 |

The gg amplitudes land at 0.95× and 1.08× the paper's Table VII(c) entries, so
this is the audited machinery driven forward rather than re-fit.

**The point of the example.** Lowest-order spin algebra fixes
Γ(³P₀→2g)/Γ(³P₂→2g) = 15/4 at a *common* wavefunction: two back-to-back gluons
are transverse, the 0⁺⁺ couples to the natural `ε₁·ε₂` structure at full
strength, while the 2⁺⁺ has to build J = 2 out of the same two polarizations and
loses 4/15 after averaging over its five substates. (The J = 1 member has no
entry at all — Landau–Yang forbids it, which is why χ_c1 is the narrow one.)

The computed ratio is 6.15, not 3.75, because annihilation is a contact process
and the two states do **not** share a wavefunction. Spin–orbit and tensor are
attractive for J = 0 and repulsive for J = 2, so the script also solves the same
channel with `V_spin = 0` to isolate the effect:

| wave | M (GeV) | r_rms (GeV⁻¹) | S₁ |
|---|---:|---:|---:|
| central (`V_spin = 0`) | 3.523 | 2.667 | 0.1432 |
| ³P₀ | 3.458 | 2.475 | 0.1733 |
| ³P₁ | 3.526 | 2.660 | 0.1452 |
| ³P₂ | 3.533 | 2.708 | 0.1366 |

A P-wave has `R(0) = 0`, so the rate rides on the slope at the origin. The χ_c0
is pulled 65 MeV below the spin-averaged level and squeezed ~7% in radius, which
lifts its smeared wavefunction-at-origin `S₁` by 21%; the χ_c2 is pushed out and
loses 5%. That is 1.61 in |S₁|², and 3.75 × 1.61 = 6.03 — the last 2% is
`α_s(M)` evaluated at two different meson masses. Fed the central wave for both,
the script reproduces the bare 3.75 exactly.

### The API this exercises

```julia
solution = fixed_channel_solution(
    params,
    ConstituentMasses(mc, mc),
    FineStructureMultiplet("P", 3, J);
    solver = OscillatorSolver(),
    nlevels = 4,
)
wave = radial_wave(solution, 1)
M = solution.eigenvalues_GeV[1]
S1  = wavefunction_origin_smearing(wave, mc; L = 1)          # Eq. (17) smeared ψ(0)
amp = gluonic_annihilation_amplitude(:P0_2g, S1, alpha_s_q(M), mc)   # amp² = Γ
A   = two_photon_amplitude(:P2, wave, mc, M, 4/9)            # Table VII(b)
```

Points worth copying:

* **Each `³P_J` is its own solve.** The J-dependence enters through
  `fine_structure_grid_operator`, and the resulting wave — not just a spin
  factor — is what the annihilation observable sees. Using one central P-wave
  for all three J would silently lock every ratio to its spin algebra.
* **Amplitudes, not widths, are the natural output.**
  `gluonic_annihilation_amplitude` returns √Γ carrying the sign of `S₁`, which
  is what Table VII tabulates and what flips once per radial node.
* **`two_photon_amplitude` has no `:P0` kind.** The paper tabulates only light
  ³P₂ two-photon rows, so the example runs the `:P2` path on the ³P₀ wave and
  rescales by the coefficient ratio `√3/√(4/5) = √15/2` — same kernel, same
  wave, different front factor. Worth promoting to a real `:P0` kind if χ_c
  γγ rows are ever added to the audit.

## `played_with_model.jl` — the scratch pad

The oldest and only uncurated example, kept because it is what
poking at this model actually looks like: display `params` and the quark-mass
table, build `Meson(all_masses, :b, :b)`, put `central_spectrum` next to the
full `compute_spectrum`, then hand the same parameters to GIPaper's
`compare_reference` and read the bottomonium residuals against the digitized
paper spectrum.

```bash
julia -e 'using Pluto; Pluto.run(notebook="examples/played_with_model.jl")'
```

It runs on the **GIPaper** environment (`Pkg.activate(joinpath(@__DIR__, "..",
"GIPaper"))`), not the examples one, because it needs the comparison layer —
which also makes it the one file here that touches reference data.

No narrative, cells in the order they were typed. Treat it as a starting point
for your own poking, not as a document.
