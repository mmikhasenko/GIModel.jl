# Examples

Worked examples of **GIModel** used as a physics tool rather than as a paper
reproduction. Each one is self-contained; the curated examples run against the
model's public API only — no reaching into internals, no reference data.

The curated examples share one environment (`examples/Project.toml`): GIModel
(dev-path to the repository root) plus CairoMakie, GLMakie, PlutoUI and
LaTeXStrings. The density support and paper renderer use GLMakie for lit
geometry CairoMakie cannot render; the other examples draw with CairoMakie.
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

## `density_candidates.jl` — saved density states

Precompute the $\rho$, $\psi(2S)$, and $\Upsilon(3S)$ states used by the paper,
and save their signed radial amplitudes and provenance without coupling the
calculation to a plotting style:

```bash
julia examples/density_candidates.jl precompute
```

The command writes portable TOML caches to `data/density_candidates/`. Each
cache records the radial grid, coefficient-weighted signed amplitudes grouped
by `(L,S)`, channel weights, model mass, rms separation, angular-normalization
and isotropy checks, solver settings, parameter contents and hash, Julia
version, and hashes of the model sources. There is no
generation timestamp, so rerunning an unchanged calculation produces identical
files and leaves the repository clean.

The $\rho$ uses the model's common light-quark constituent mass. The three
states are predominantly $1^3S_1$, $2^3S_1$, and $3^3S_1$; their allowed
$D$-wave mixing channels are included before export. The cache stores numerical
wave-function data rather than Julia serialization, so the paper renderer can
load it without solving the spectra again.

The paper-specific layout, palette, cutaway, and caption live in the separate
report repository under `report/sources/`. See `report/README.md` for the exact
figure-production command.

## `density_3d.jl` — what a meson looks like

A script that draws the quark-antiquark **separation density** as a Monte-Carlo
point cloud, one dot per draw from `|psi(r)|^2`, in the style of the dot-cloud
pictures of atomic orbitals.

```bash
julia examples/density_3d.jl
```

**What it is not.** `r` is the separation of the quark and the antiquark, not a
position in the lab, and the cloud is not a charge distribution.

**The one decision that makes the figure.** Fix `mJ`. Summing over `mJ` gives
`sum_m |Y_Lm|^2 = (2L+1)/4pi`, an exact sphere for every state in the spectrum —
so an `mJ`-averaged picture of a meson is a ball, always. At fixed `mJ` the
panels separate: `1^1P_1` puts its lobes along `z` at `mJ = 0` and becomes a
torus at `mJ = 1`: the *same level* in two different magnetic substates.

**Reading the density, not the radial probability.** The cloud is `|psi|^2`, the
density per unit volume; `u(r)^2` is the radial probability per unit `r`. They
look nothing alike — an `r^2`-weighted 3D cloud would hollow out every S-wave at
the origin. The `1/r^2` in `psi = (u/r) Y_Lm` cancels against the `r^2` of the
volume element, which is why the radial marginal here is exactly `u(r)^2`.

**Reading the cutaway.** Two vertical half-planes meet on the `z` axis and
remove a 120° wedge facing the camera, like open doors. The retained samples
stay at their calculated positions. In the 2S comparison, the opened cloud is
an exact subset of the full cloud: removing the wedge does not refill or
resample the retained region. The six-panel figure uses a common camera,
particle radius, and accepted sample count. Scale bars lie in the image plane
through the origin; perspective sizes vary with depth.

**Lighting and depth.** GLMakie draws opaque spheres with a short-range point
light at the camera, a dim fill, and screen-space ambient occlusion. The light
uses intensity 8 and range `6 × half_width`; the previous intensity 10 and range
`14 × half_width` saturated the particles. The camera uses a fixed 55°
perspective view, which makes the spheres large enough to read as texture.
Rendering uses an invisible OpenGL window and still requires a working graphics
context.

A headlamp alone does not give the cavity enough contrast: points behind a
node remain visible through it. The cutaway therefore also applies an
**illustrative geometric depth tint**, `0.24 + 0.76 exp(-d / 0.16 fm)`, where
`d` is distance behind the nearest exposed half-plane. This is a rendering cue,
not a calculated shadow or a density colour scale. It depends only on the cut
geometry, never on radial node locations or wave-function values. No separate
slab or background shell is added. Particle concentration carries the
probability information; brightness helps distinguish the exposed faces from
the interior. Finite sphere size and occlusion limit quantitative readings
from the cloud itself.

**Checking the visible node against the calculation.** The 2S figure includes
the full mixed state's radial marginal `p(r) = sum_(L,S) U_LS(r)^2`, normalized
per fm. A dashed line marks the interpolated zero of the dominant S-channel
amplitude. The D-wave admixture can leave a small nonzero probability there,
so the shaded ring is described as radial depletion. The profile uses the
solution directly, independently of the Monte-Carlo draw and its shading.
Parameters are the repository's provisional GI model parameters.

### Composite states

`physical_components` returns the flattened `(n, L, S)` composition of a mixed
state with signed coefficients, and components that share `(L, S)` and differ
only in `n` must be **added as amplitudes before squaring** — they carry the same
angular function, so only their `r`-integral is orthogonal, not their value at
each `r`. Squaring them separately is a silent factor-level error; the script
groups by `(L, S)` first.

Interference between different `L` at the same `S` — the `3S1`/`3D1` tensor
mixing — survives the sum over the unobserved quark spin, because every `L` in
one spin group carries the same `m_L = mJ - m_s`.

### Two checks that need no reference data

Both are exact identities. The angular integral is evaluated numerically, so
its residual includes the finite quadrature error:

| check | identity |
|---|---|
| `check_angular_normalization` | `2 pi int f(x) dx == p(r)` at every radius, from orthonormality of the `Y_Lm` |
| `check_mj_sum_isotropy` | the `mJ`-summed density is spherical, since `sum_mJ \|J mJ><J mJ\|` commutes with rotations |

The first constrains the Clebsch-Gordan coefficients only through
`sum_ms CG^2 = 1`; the second pins down each one individually. The script prints
both per state (`~1e-7`, the angular quadrature floor, and `~1e-15`) and
asserts tolerances of `1e-6` and `1e-12`, respectively.

### The API this exercises

```julia
spec = compute_spectrum(params, ccbar; levels, solver = OscillatorSolver())

for piece in physical_components(spec, "1^3D_1")
    piece.basis.L_label, piece.basis.multiplicity, piece.coefficient
    u = sample_wave(piece.wave, rgrid).u     # OscillatorWave -> plotting grid
end
```

`sample_wave` is the sanctioned way to put an `OscillatorWave` on a mesh: the
oscillator solver keeps no spatial grid, and the sampled result is normalized on
the plotting grid, which is what a density needs. The `MeshWave` branch in
`u_samples` interpolates instead, so the script gives the same picture under
either solver.

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

The annihilation-width example also uses `QuarkModelTransitions`, which owns
the observable kernels; spectrum-only examples need only GIModel.
