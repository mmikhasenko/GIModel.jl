# Review: transition matrix-element API plan

Status: review of
[matrix_element_api_implementation_plan.md](matrix_element_api_implementation_plan.md)
and [transition_matrix_element_api.md](transition_matrix_element_api.md),
2026-09-20.

Reviewed against `src/strong_decays.jl`, `src/spectrum.jl`,
`src/mock_meson_overlaps.jl`, `src/model_objects.jl`,
`GIPaper/src/table_v.jl`, the canonical Table V CSV, and the paper's
Appendix B/C (`paper/vision_ocr/pages/page-012.md`, `page-039.md`,
`page-040.md`).

## Summary

The plan's *policies* are right and should be kept verbatim: the invariant
list, the "what should remain hardcoded" section, the refusal to freeze Layer D
before two operator prototypes exist, the symbolic-generation policy, and the
risk register. Those are the parts that are hard to get right and they are
already good.

The plan's *structure* has one load-bearing problem and several smaller ones.
The load-bearing problem is that it makes the partial wave an **input** to
`matrix_element`, when in both the paper's own formulation and the modern
literature the computational primitive is the helicity amplitude and the
partial waves are a **projection of it** — several at once, sharing one spatial
integral. Fixing that changes Layers C/D, the `TwoMesonChannel` type, and the
return type, so it is worth settling before Phase 1 rather than after.

Sections 1–3 below are the structural changes. Section 4 is smaller cleanups.
Section 5 answers the plan's own review questions where the code or the paper
already settles them.

---

## 1. The primitive is the helicity amplitude, not `M^{LS}`

### What the paper does

Eq. (19) is introduced with its own evaluation recipe (page-image 12):

> The calculations are most readily performed by taking `q = q ẑ` thereby
> calculating helicity amplitudes `H_m` where `m = s' = s`, and then
> transforming to the usual partial-wave basis; details of this process are
> given in Appendix C.

Appendix C then gives

- (C1) `H_m = A[M*_{j*m} → M_{jm} + P(q ẑ)]`;
- (C2) `Γ = 1/(2j*+1) · (q/2π) · Σ_m |H̃_m|²`, with
  `H̃_m = (2π)^{9/2} (2M* 2E)^{-1/2} H_m`, and the explicit non-relativistic
  choice `E/M* → 1`;
- (C3) `h_m = √2 H_m` for `m>0`, `h_0 = H_0`;
- Table XI: the helicity → partial-wave conversion for `M*_j → V + P`, keyed by
  `j*` and parity. For a spin-zero daughter, `A_L = h_0`.

Table XI is the decisive detail. One parent gives **two** partial waves from
the **same two** helicity amplitudes:

```
j* = 1, positive parity:   A_S = √(2/3) h_1 + √(1/3) h_0
                           A_D = √(1/3) h_1 - √(2/3) h_0
```

That is not incidental — it is why Table V lists
`A1 → [(ππ)_ρ π]_S` and `A1 → [(ππ)_ρ π]_D` as separate rows of one decay, and
why the `S` and `D` reduced classes exist as a pair in Table IV.

### What the plan proposes

```julia
final = TwoMesonChannel(d1, d2; relative_L = 0, channel_spin = 0)
amp = matrix_element(final, operator, parent; kinematics = OnShell())
```

`relative_L` and `channel_spin` are user inputs, and one call returns one
scalar `total`. Consequences:

1. **The user supplies something derivable.** For `M* → V + P` the channel spin
   is `S = 1`, always; for `P + P` it is `S = 0` and `L = J`. The free
   `channel_spin` argument invites combinations that are identically zero or
   forbidden, which the plan then has to catch with a selection-rule error. The
   allowed `(L,S)` set follows from `J^P` of the three states; it should be
   enumerated, not accepted.
2. **The expensive part is recomputed per partial wave.** `A_S` and `A_D` above
   are two linear combinations of the same `h_0, h_1`, which are two evaluations
   of the same spatial integrals at the same `k`. Calling `matrix_element` twice
   does that work twice, and — worse — makes the two results look independent
   when they are rigidly correlated.
3. **There is no layer for the projection.** Layer C lists `topologies`,
   `selection_rule`, `spin_factor`, `flavor_factor`, `color_factor`,
   `spatial_problem`. None of those is Table XI. The plan mentions
   "Jacob–Wick partial-wave projection" once, under Layer 7's algebra strategy,
   but never places it in the dispatch architecture.
4. **It contradicts the manifest.** `docs/paper_manifest/appendix_cd.toml`
   currently marks C1 `folded` and C3 `context`, with the note "Table IV already
   tabulates the partial-wave amplitudes, so we work directly in the
   partial-wave basis … rather than projecting from helicity amplitudes." That
   is correct *for the legacy backend* and becomes false the moment Phase 4
   evaluates Eq. (19) natively. Phase 4 un-folds C1, C3, and Table XI, and the
   plan should say so and budget for it.

### Recommendation

Make the channel the pair of daughters only, and return the whole partial-wave
set:

```julia
final = TwoMesonChannel(d1, d2)          # canonicalized, symmetrized if identical

amp = matrix_element(final, operator, parent; kinematics = CMKinematics(k))

amp.helicity          # h_0, h_1, ... the primitive
amp[PartialWave(0,1)] # A_S, a projection of it
partial_waves(amp)    # the allowed (L,S) list, derived
decay_width(amp)      # sums over the set exactly once
```

with an explicit, testable projection layer:

```julia
allowed_partial_waves(final, initial)               # from J^P + identical-particle symmetry
partial_wave_projection(final, initial)             # the Table XI / Jacob-Wick matrix
```

This is strictly more general (it is the input a Dalitz-style analysis or a
self-energy sum wants anyway), it removes two user arguments, it makes the
`1/√2` identical-particle factor derivable where the plan wants it, and it puts
Appendix C where it belongs instead of inside an operator.

The `L,S`-per-call form can survive as a convenience wrapper — `matrix_element(final, op, initial; partial_wave = PartialWave(0,1))` — but it must be the wrapper, not the primitive.

---

## 2. Factorize as (pure coefficients) × (basis of spatial integrals)

### The plan's seam

```julia
spatial_kernel(op, topology, final_basis, initial_basis, kinematics)
spatial_matrix_element(kernel, wa, wb, wc)
```

Both signatures carry `final_basis` and `initial_basis` into the spatial layer.
That means the spatial result is keyed by the full spectroscopic identity of
the external states, so nothing is shared between two rows that differ only in
`J`, or between the `S` and `D` projections of one decay, or between the
components of a mixed state that share a radial wave.

### How the literature factorizes it

Roberts and Silvestre-Brac (Few-Body Syst. 11, 171 (1992)) and, in its
cleanest modern form, Burns (arXiv:1403.7538 — already cited in
`quark_model_resonances.md`) show that for both `³P₀` and elementary-emission
operators the spin-space algebra factorizes into coefficients that are
**independent of the radial wave functions**, multiplying a small basis of
orbital matrix elements. Schematically

```
M = Σ_topologies Σ_μ  ξ_μ(quantum numbers) · I_μ(k; w_A, w_B, w_C)
```

where `μ` runs over a handful of orbital/multipole labels, `ξ_μ` are pure
numbers (rationals and square roots of rationals), and `I_μ` are the integrals.
The same structure is visible in Table IV of the paper: the classes
`A, A′, A″, A₀, S, D, P` *are* the small basis, and the per-row `coefficient`
column *is* `ξ`.

### Recommendation

Split Layer C/D at that line instead:

```julia
# pure numbers; no wave functions; exactly testable against Table V's coefficient column
orbital_decomposition(op, topology, final_basis, initial_basis)
    -> Vector{Pair{OrbitalLabel, Float64}}

# expensive; keyed only by (waves, label, k); cacheable and shared
spatial_integral(op, label::OrbitalLabel, wA, wB, wC, k)
```

Payoffs, in order of importance:

1. The generated algebra becomes a pure function of small integers, so the
   Phase 3 gate ("exact signs and magnitudes match the paper convention") is a
   test on rational arithmetic, with no numerics in the way.
2. One spatial integral serves every `(L,S)`, every `J` in a multiplet, and
   every mixed-state component pair that shares radial waves. For the `³P₀`
   convolution in Phase 5 this is the difference between a usable and an
   unusable vertex grid for WP2's self-energy integrals.
3. The Table IV classes stop being a legacy-only concept: they become the
   `OrbitalLabel` set for the elementary-emission operator, which makes the
   legacy backend a genuine special case of the general one rather than a
   parallel implementation.

---

## 3. Replace `StateView` with a resolved `PhysicalState`

### The problem `StateView` solves, and a better solution

`StateView{S,T}` exists because `physical_components` needs the owning
`Spectrum` — waves live in `Spectrum.computation.channel_cache`, deliberately,
not on the state (`docs/code_architecture.md`). So "a state" is genuinely a
`(spectrum, state)` pair, and the plan is right to bundle them.

But a *lazy view* keeps the `Spectrum` alive through the whole calculation, and
`matrix_element` then has to re-enter `physical_components` and re-derive the
mass. Resolving once is simpler and gives more:

```julia
struct StateComponent{B,C,W}
    basis::B
    coefficient::C
    wave::W
end

struct PhysicalState{C}
    label::String
    mass_GeV::Float64
    components::Vector{StateComponent{...}}
    provenance::...
end

physical_state(spec, "1^3S_1") -> PhysicalState
```

Then `matrix_element(final, op, initial)` takes three `PhysicalState`s and
**never sees a `Spectrum`**. That:

- answers the "three separately owned spectra" requirement structurally rather
  than by convention — there is nothing left to disagree about;
- makes the invariant "one physical mass per external state" a *field* rather
  than a rule someone must remember not to break (the plan currently states it
  three times in prose, which is a sign it has no structural home);
- makes every layer above the spectrum testable with synthetic states — the
  Phase 2 exit gate asks for "synthetic tests covering … simultaneous
  three-state mixing", which is painful to construct through a `Spectrum` and
  trivial with a hand-built `PhysicalState`;
- gives `explain` everything it needs from the amplitude alone.

This subsumes the plan's Layer A entirely and most of Layer B. `physical_components`
keeps working as the lower-level accessor; `physical_state` is
`physical_components` + mass + identity, which is a few lines.

It also answers review question 2: `StateView` should not be public, because it
should not exist.

### On `physical_decay_amplitude` (question 3)

Keep it internal. Its signature
`physical_decay_amplitude(kernel, parent, daughter1, daughter2)` is already
specialized to "one ket, two bras", which is the only case in scope; a public
name promises a generality it does not have. If a third case appears later, the
right public object is the composer over `PhysicalState`s, not a function named
after strong decay.

---

## 4. Smaller points

### 4.1 The deprecation ceremony is oversized for this package

`matrix_element(a::StrongDecayAmplitude)` has exactly two call sites outside its
own definition: `test/strong_decays.jl:107` and one line of
`docs/strong_decays_tutorial.qmd`. The package is at `0.2.0`, where SemVer
permits breaking changes on a minor bump.

The plan's five-step sequence (add alias → warn → add 3-arg method → update docs
→ "remove only at a breaking release boundary") buys nothing and costs a
permanently ambiguous name during the exact period when the whole point of the
work is to make the name unambiguous. Recommend: rename to
`reduced_matrix_element` in one commit, no alias, bump to `0.3.0`, and delete
step 5 from Phase 7.

The *rest* of Phase 0 — the units audit and the factorization tests — is real
and should stay. The units are indeed misdocumented: as the code factorizes it,
`coefficient` and `reduced` are both dimensionless (the fit in
`calibrate_strong_decay_model` divides `MeV^{1/2}` by `MeV^{1/2}`), only
`spatial_overlap` and `total` carry `MeV^{1/2}`, and the `StrongDecayAmplitude`
docstring currently labels all four `MeV^{1/2}`.

### 4.2 `convention::Symbol` should become a dispatch singleton

`:table_iv` / `:leading` is threaded as a `Symbol` keyword through
`reduced_decay_amplitude`, `strong_decay_amplitude`, `calibrate_strong_decay_model`,
and `decay_amplitude`, validated by hand at the top of one of them. The plan
says these must be "explicitly named approximations, not hidden keywords" but
does not list the change.

The repo already has the idiom and the precedent: `CentralPotentialMethod`
singletons replaced five booleans and a TOML string. Do the same:

```julia
abstract type ReducedAmplitudeConvention end
struct TableIVPolynomial <: ReducedAmplitudeConvention end   # the printed formula
struct LeadingS0         <: ReducedAmplitudeConvention end   # the paper's numeric column
```

This is cheap, is confined to `strong_decays.jl`, removes the hand-rolled
validation, and makes the approximation visible in every amplitude's type.

### 4.3 Concrete fields in `TransitionAmplitude`

The two documents disagree — `{T,O,I,F,K}` in the design, `{T,O,I,F,K,N}` in the
plan — and both leave `terms`, `normalization`, and `provenance` untyped. Every
parameter is anchored to a field, which satisfies the repo's rule, but the
untyped fields undo the benefit. Type them:
`terms::Vector{TransitionTerm{T}}`, `provenance::TransitionProvenance`.

With `PhysicalState` from §3, `I` and `F` are concrete already, and
`kinematics`/`normalization` are small singletons — so the parameter list can
shrink rather than grow.

Separately, the plan's width dispatch:

```julia
decay_width(::TransitionAmplitude{<:Any,<:Any,<:Any,<:Any,<:Any,GITableVNormalization})
```

is five wildcards to reach the sixth parameter. Either put the normalization
first, or — better — dispatch on the normalization object itself:

```julia
decay_width(amp) = partial_width(amp.normalization, amp)
partial_width(::GITableVNormalization, amp) = ...
```

Same behaviour, one readable signature per convention, and an unsupported
normalization gives a clean `MethodError` naming the type.

### 4.4 `explain` duplicates `show`

`StrongDecayAmplitude` already has a three-argument `show` that prints the
factor table with `[PAPER]` / `[DERIVED]` tags, and it is the best piece of
discoverability in the current decay code. Adding `explain(amp)` as a new
exported verb means a user who types `amp` at the REPL gets the short form and
never learns the long one exists.

Recommend: make `show(io, MIME"text/plain", amp)` the full decomposition, as it
is today, and add `explain` only if there is genuinely more than a REPL display
can hold (per-topology breakdown for a mixed state with many terms, say) — in
which case name what it adds, e.g. `term_table(amp)`.

The `operators(:strong_decay)` registry is fine and worth keeping; it has no
`show` equivalent.

### 4.5 Export budget

`length(names(GIModel))` is currently **169**. Direction 1 of the July plan was
*fewer, cleaner exports* (64 at the time). The API plan adds roughly eighteen
more names.

Recommend exporting: `matrix_element`, `reduced_matrix_element`, `decay_width`,
`TwoMesonChannel`, `PartialWave`, `physical_state`, the operator types, the two
kinematics types, the normalization types, `TransitionAmplitude`. Keep
internal until a second operator exists: `topologies`, `selection_rule`,
`spin_factor`, `flavor_factor`, `color_factor`, `orbital_decomposition`,
`spatial_integral`, `TransitionTerm`, `StateComponent`. The plan already applies
this reasoning to Layer D; extend it to Layer C.

### 4.6 Below-threshold behaviour must be decided in Phase 1, not Phase 5

`decay_momentum` returns `0.0` below threshold and `spatial_overlap` returns
`0.0` for `q ≤ 0`. That is correct for a width table and directly contradicts
the plan's own invariant that nothing may silently return a documented zero
where an off-shell vertex is wanted: a self-energy integral (WP2, and
"definition of done" item 10) needs the vertex on the closed-channel side,
where `k` is not real and the answer is not zero.

`CMKinematics(k)` is introduced in Phase 1 but the closed-channel convention is
never specified. Decide there: either `CMKinematics` takes complex `k` and the
zero guard moves into `OnShell` only, or a separate `SubThreshold` kinematics
type carries the analytic continuation. Deferring it to Phase 5 means Phase 1's
type is wrong for its stated purpose.

### 4.7 Quasi-two-body daughters are in Table V, not out of scope

The plan puts "three-body final states and unstable-daughter spectral
functions" out of scope, which is right for the *physics*. But the canonical CSV
has rows with `ε`, `δ₂`, and `κ` daughters — that is what the `A′`, `A″`, `A₀`
classes and the `[A0/A]` realistic-factor bracket exist for
(`A1 → (ππ)_ε π`, `Q1 → (Kπ)_κ π`, `D → (ηπ)_{δ₂} π`, …). The legacy backend
must keep producing them, so the typed `TwoMesonChannel` has to be able to hold
a daughter that is a *label with a mass* and no radial wave.

This is a required Phase 1 decision (does `TwoMesonChannel` accept a
non-`PhysicalState` daughter, and does the type system stop you from asking a
native operator to evaluate one?), not a Phase 6 concern.

### 4.8 The Phase 3 target set should include a `mixing_only` row

The proposed initial set is good. Add one row the CSV marks `mixing_only` —
`f′ → ππ` (`+1.1` observed, coefficient exactly `0` in the pure basis) or
`E → δ₂π`. These are the sharpest possible test of the coherent composer,
because the pure-component amplitude is exactly zero and the physical amplitude
is not. Nothing else in the target set distinguishes a correct coherent sum from
an incoherent one as cleanly.

### 4.9 Name the existing gate instead of inventing a new one

Phase 0 asks to "snapshot current Table V outputs". The project already has
that gate: `bash scripts/verify_project.sh` runs both suites plus
`GIPaper/scripts/reproduce_table_v.jl`, and the standing acceptance criterion
for refactors in this repo is byte-identical residual reports. Phase gates
should cite it rather than describe a parallel mechanism.

---

## 5. Answers to the plan's review questions

Where the code or the paper already settles a question, it should be closed in
the plan rather than left open.

**Q1 — argument order.** `matrix_element(final, operator, initial)` is right and
matches the repo: `radial_overlap(wx, wy, f)`,
`physical_transition_amplitude(kernel, spec, left, right)`,
`mock_meson_overlap(mwx, mwy, …)` all put the bra first. Keep it.

**Q2 — should `StateView` be public.** It should not exist; see §3.

**Q3 — should `physical_decay_amplitude` be public.** Internal; see §3.

**Q4 — canonical normalization.** The relativistic two-body normalization, with
`GITableVNormalization` as a reference backend carrying a documented conversion.
The paper forces this: (C2) sets `E/M* → 1` explicitly ("consistent with our
nonrelativistic calculation"), so the Table V convention cannot be used under a
loop integral, which is the stated endpoint. The conversion factor is written
down in C2 — `H̃_m = (2π)^{9/2}(2M* 2E)^{-1/2} H_m` — so this is a derivation,
not a fit.

**Q5 — first native operator.** GI Eq. (19), not `³P₀`, and not both. Four
reasons: it is the only one with a 220-row numeric oracle inside this repo, so
every layer above the spatial integral can be validated before new physics
enters; its conventions are fully written down (Appendix B flavor operators
B22–B25, spin functions B26–B29, spatial phases B31–B37, Appendix C projection);
it is a one-body operator between two waves, so it reuses `radial_overlap`,
`momentum_functional`, and `momentum_overlap` with no new numerics; and it
leaves the risky shifted convolution for after the algebra, composition, and
normalization layers are proven. Phase 4/5 ordering already implies this — close
the question.

**Q6 — minimal defensible derived set.** The plan's six, plus one `mixing_only`
row (§4.8).

**Q7 — angular-algebra dependency.** Checked the registry against the plan's
own criteria:

| package | 3j/CG | 6j | 9j | deps | note |
|---|---|---|---|---|---|
| `PartialWaveFunctions.jl` | yes | no | no | **none** | your own; Condon–Shortley; `_doublearg` integer API; julia ≥ 1.2 |
| `WignerSymbols.jl` | yes (exact) | yes | no | Primes, LRUCache, HalfIntegers, RationalRoots | exact rational-root arithmetic |
| `Wigxjpf.jl` | yes | yes | yes | `Wigxjpf_jll` binary | the only one with 9j |

No package supplies all three with a small footprint. Recommendation:
`PartialWaveFunctions.jl` for CG and Wigner-d (zero dependencies, and its
doubled-integer convention matches the "combine in the `L·S` order with Wigner's
conventional Clebsch-Gordan coefficients" statement of B37), and generate the
6j/9j set offline under the plan's own deterministic-generation policy. For
`L ≤ 4` and `S ≤ 1` that table is small, and generating it removes the
dependency question rather than answering it.

Note that no library settles the *phase* question, which is the one that
actually matters here. Three of the ledger entries the plan asks for already
exist and only need collecting: `fix_outer_phase` (outer lobe positive),
`fix_annihilation_phase!` (positive overlap with the ascending-mass precursor),
and the oscillator-basis phase fix. A fourth needs checking rather than
asserting: B31–B36 give `Ψ_{0LL}` an alternating sign following `sgn Y_{LL}`
while fixing radial excitations positive at large `r`. The radial half agrees
with `fix_outer_phase`; the `(-1)^L` lives in the angular factor and should be
supplied by Condon–Shortley `Y_{LM}` in the generated algebra. Confirm that it
is applied exactly once before writing it into the ledger as a convention.

**Q8 — identical mixed final states.** Canonicalize in the `TwoMesonChannel`
constructor, and derive the exchange phase from `(-1)^{L+S}` once the allowed
`(L,S)` set is enumerated (§1). The manual `1/√2` row factors then become a
consequence rather than data. This is only possible with the §1 change.

**Q9 — HO/FD tolerance.** Cannot be answered from the desk; it must come from
the same kind of convergence certificate the repo already produces
(`audit_ho_convergence.jl`, `audit_fd_convergence.jl`). State that as the
method and set the number from a node-sensitive channel's certificate.

**Q10 — which `A,S0` behaviour is reproduction-only.** `A`, `S0`, `β = 0.40`,
the `:leading` convention, and the Table IV class polynomials are all
reproduction-only: they are a two-parameter fit to two rows under a single-`β`
SU(6) reduction. The `g`, `h` of Eq. (19) are the operator parameterization.
The relation between them is a calibration study (`S₀ = 3hβ` is the paper's
own identification, and the `S` class carries `-(1/2)(g + h/4)q²/β²`), so the
two couplings are in principle recoverable from the two fitted numbers — but
only within the SU(6)/single-`β` reduction that Phase 4 is replacing. Treat any
agreement as evidence, not identity, exactly as the risk register says.

---

## 6. What to change in the plan document

Minimal edit list, in dependency order:

1. §5 — `TwoMesonChannel` loses `relative_L`/`channel_spin`;
   `matrix_element` returns a partial-wave set; add
   `allowed_partial_waves` and `partial_wave_projection` to the public API.
2. §6 — add "Layer C½: helicity → partial-wave projection (Appendix C,
   Table XI)"; split Layer C/D at the `orbital_decomposition` /
   `spatial_integral` line of §2.
3. §6 Layer A/B — replace `StateView` with `physical_state` → `PhysicalState`;
   delete the Layer A/B split.
4. §8 Phase 0 — replace the five-step deprecation with a single rename at
   `0.3.0`; cite `scripts/verify_project.sh` as the gate; add the
   `convention::Symbol` → singleton change.
5. §8 Phase 1 — add the below-threshold/complex-`k` decision and the
   quasi-two-body-daughter decision to the exit gate.
6. §8 Phase 4 — state that it un-folds C1, C3, and Table XI in
   `docs/paper_manifest/appendix_cd.toml`.
7. §13 — close questions 1–7 and 10 with the answers above; leave 8 and 9 open
   until the §1 change lands and a convergence certificate exists.
8. Throughout — replace `explain(amp)` with the existing `show` idiom, and add
   an export budget line to each phase gate.
