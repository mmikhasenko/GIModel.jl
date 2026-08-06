# Engineering Work Plan

> **Scope.** This is the *code* stream, not the physics one. Per-unit paper
> status is the manifest (`docs/paper_manifest/*.toml`); the physics rationale
> behind `partial`/`missing` is [`paper_gap_ledger.md`](paper_gap_ledger.md);
> the reproduction process record is [`work_plan.md`](work_plan.md). This file
> is what a newcomer needs to continue the architecture work: what is settled,
> what is next, and which mistakes are already paid for.

## Settled invariants — do not silently break these

These were each found the hard way. Every one is pinned by tests.

1. **One normalization.** Every solve returns `u` with `∫u² dr = 1`, via the
   single definition `physically_normalized_waves`. Finite-difference used to
   return Euclidean eigenvectors (`Σu² = 1`) while the oscillator path returned
   physical ones — they differ by `√h`, and anything quadratic in `u` given the
   wrong convention is off by `h ≈ 0.05`.

2. **The wave interface is the only way to get a `u`.** Seven operations on
   `RadialWave` (`wave_norm`, `radial_expect`, `radial_overlap`,
   `momentum_wave`, `momentum_expect`, `momentum_functional`, and `smear` when
   it lands). Operations normalize internally, so consumers must not.
   *Consumers* go through the interface; *implementations* of those operations
   legitimately touch `.u`, `.r`, `.h`.

3. **Helpers that are quadratic in a wave must normalize themselves.** An
   unstated "callers must pre-normalize" requirement is invisible at the call
   site and breaks the moment one caller stops. Caught by a scale-invariance
   probe: any wave consumer must give the same answer for `u` and `13.7·u`.

4. **The oscillator basis has a fixed phase.** LAPACK's QR assigns arbitrary
   column signs — invisible numerically (a sign flip is unitary) but fatal to
   any closed-form matrix element. `orthonormalize_physical_basis` now keeps
   each column's original sign. Same lesson as `fix_annihilation_phase!`.

5. **An unimplemented path fails, it does not fall back.** Returning an empty
   result that a caller reads as "use something else" produced a 3× wrong pion
   (0.2842 vs 0.0950 GeV) with no warning.

6. **Eq. (A17) is atomic.** An exact momentum side with a mesh-projected
   potential is not the Hamiltonian of any single problem, is not variational,
   and pushed Table VII gluonic ratios out of band. Both sides move together or
   neither does.

## Verification standard

Every stage: `bash scripts/verify_project.sh` exit 0 **and** all residual
reports byte-identical apart from the timestamp lines. Check `VERIFY_EXIT` before
trusting a report diff — `set -e` aborts early and leaves stale reports that
compare clean against themselves.

Prefer **invariants over recorded numbers**: closed forms the code must
reproduce (`⟨r⟩ = 3/2` for `u = re^{-r}`), Parseval, level orthogonality,
scale invariance, virial theorem, convergence under refinement. These test the
physics; a pinned number only tests that nothing changed.

Validate self-contained mathematics **in isolation first**. Oscillator matrix
elements are a well-posed problem with no reference data needed — `p²` and `r²`
reconstructing the oscillator Hamiltonian to 1e-14 is worth more than watching a
Table VII ratio drift.

## Next stages

### 1. `ho_operator_matrix` performance (do first)

`oscillator_hamiltonian_for_beta` calls it twice per `(L, β)` — once for `G̃`,
once for `S̃` — and each runs its **own** adaptive doubling loop over the same
quadrature nodes. Across a 22-β scan that duplicates every eigendecomposition.

*Fix:* compute the decomposition once per `(L, β, nq)` and evaluate both
operators from it; consider reusing smaller-`nq` results rather than restarting
the doubling. Roughly 2×, no numerical change.

*Why first:* the full gate went from ~6 to ~20+ minutes with A17, which makes
every later stage uncomfortable to iterate on. Correctness is unaffected.

*Gate:* zero report drift.

### 2. B3 — solver types

Split `RadialSolver` into `FiniteDifferenceSolver` and `OscillatorSolver`. Now
justified rather than cosmetic: since A17 landed, the oscillator path has **no
operator mesh**, so `OscillatorSolver` genuinely has no `ngrid`/`rmax` for
operators (only a reconstruction mesh for reporting) and does have `nbasis` and
`beta_grid` — which are currently module constants that the β-railing warning
tells users to go edit in source.

Also drop the phantom `Basis` type parameter from `GIParameters` so the physics
struct stops carrying a marker about how it will be discretized. ~19 `with_basis`
call sites plus 7 dispatch sites.

*Watch:* `contact_hyperfine_nonperturbative_states` is now basis-generic, so the
dispatch it used to rely on is gone — check nothing else depends on it.

*Gate:* zero report drift.

### 3. C2 — report provenance

No residual report records the grid or solver that produced it. A report saying
`mean_abs = 6.0 MeV` cannot tell you whether that was 450 or 900 points, which
is a real gap for a reproduction project.

*Cost:* a one-time header change to **all 24 reports** — the only sanctioned
non-zero-drift stage outside A17. Review the baseline diff once, accept, done.

### 4. C3 — retire the deprecated keywords

`ngrid`/`rmax`/`kinetic`/`eigensolver`/`nlevels_per_channel` and the four spin
switches still work and still override the `solver`/`terms` objects, warning
once each. ~180 call sites, but only **two files in `GIPaper/src`**
(`comparison.jl`, `residual_report.jl`) — start there, it silences the warning
for the main comparison path.

Mechanical, no deadline, driven by the warnings themselves.

## Known open items, deliberately not fixed

- **`src/` is not fully flavor-free.** `quark_mass_table.jl:16-17` maps `"c"`/
  `"b"` to TOML keys, and `is_charm_class` gates the `DecayChannel` constructor.
  The dynamics are flavor-blind; the vocabulary is ~90% there.
- **Strange-row form factor.** The paper applies its unequal-mass Gaussian to
  charmed rows only, though strange rows have `r = 0.66`. Recorded as a
  candidate discrepancy rather than "fixed" — unlike the `S_c` coefficient, no
  printed row says the paper is wrong there.
- **Mixing blocks still run on FD waves** (annihilation, tensor/spin-orbit
  off-diagonal) — the remaining bounded gap in A17's manifest note.
- **Comparator central methods keep the mesh** deliberately: several smear
  numerically on it and are not closed-form functions of `r`, so making them
  hybrid would violate invariant 6.

## Diagnostics that do not work

Recorded because they cost real time:

- **Comparing the QR output `U` against the raw basis `B`** to test a phase
  convention. Both come from the same routine, so a flip *inside* the basis is
  invisible — it reports a uniform `-1` that cancels in a matrix
  (`S·ana·S = ana`) and makes the problem look absent. Compare element by
  element instead.
- **Watching a scored ratio drift out of band** as the primary signal that
  something is wrong. It says only "something, somewhere, downstream of
  everything". Validate the isolated mathematics first.
- **Forming oscillator matrix elements as (weight) × (polynomial).** At large
  quadrature order the two factors are ~1e-95 and ~1e19; the digits are gone
  before summation starts, so compensated summation cannot recover them. Use
  the Golub–Welsch/DVR form, where the eigenvector entries *are* the product.
