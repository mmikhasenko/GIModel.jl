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

7. **The solver carries the method; the parameters carry the model.**
   `GIParameters` has no basis marker, and `RadialSolver`
   (`FiniteDifferenceSolver` / `OscillatorSolver`) decides how the radial problem
   is solved. A stage-1 solver is recorded in `SectorComputation` so stages 2-3
   resolve their own eigenproblems the same way — otherwise one spectrum would be
   two calculations. Settings a method does not have (`kinetic`, `eigensolver` on
   the oscillator path) are absent and throw, never silently ignored.

## Verification standard

Every stage: `bash scripts/verify_project.sh` exit 0 **and** all residual
reports byte-identical apart from the timestamp lines. Check `VERIFY_EXIT` before
trusting a report diff — `set -e` aborts early and leaves stale reports that
compare clean against themselves.

Two ways that gate has been quietly wrong, both now fixed:

  - **Ordering.** `data_checks.py score-annihilation` does not compute its
    numbers, it *parses* them out of `table_iii_mixing_audit.md` and
    `isoscalar_residuals.md`. It used to run first, so the scorecard published
    the previous run's values and every change surfaced as drift one gate late,
    blamed on whatever was in flight then. Anything that reads a report must run
    after the script that writes it.
  - **Coverage.** A script the gate does not run will rot and nobody will know.
    `audit_nonmixing_contact.jl` was broken outright by the GIModel/GIPaper split
    and sat dead for a month with its reports frozen at that commit. Both FD/HO
    audits are in the gate now. If you add a report, add its script.

Prefer **invariants over recorded numbers**: closed forms the code must
reproduce (`⟨r⟩ = 3/2` for `u = re^{-r}`), Parseval, level orthogonality,
scale invariance, virial theorem, convergence under refinement. These test the
physics; a pinned number only tests that nothing changed.

Validate self-contained mathematics **in isolation first**. Oscillator matrix
elements are a well-posed problem with no reference data needed — `p²` and `r²`
reconstructing the oscillator Hamiltonian to 1e-14 is worth more than watching a
Table VII ratio drift.

## Next stages

### 1. C2 — report provenance

No residual report records the grid or solver that produced it. A report saying
`mean_abs = 6.0 MeV` cannot tell you whether that was 450 or 900 points, which
is a real gap for a reproduction project.

*Cost:* a one-time header change to **all 24 reports** — the only sanctioned
non-zero-drift stage outside A17. Review the baseline diff once, accept, done.

### 2. C3 — retire the deprecated keywords

`ngrid`/`rmax`/`kinetic`/`eigensolver`/`nlevels_per_channel` and the four spin
switches still work and still override the `solver`/`terms` objects, warning
once each. ~180 call sites remain. `GIPaper/src/comparison.jl` is already done —
`compare_reference` takes a `solver` and **throws** if given both it and the
loose mesh keywords, since one of the two would have to be discarded. That is the
pattern to copy; `residual_report.jl` still prints a hardcoded
`"finite-difference"` in its header, which C2 should fix at the same time.

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

- **Diagnosing a numerical failure by its most dramatic-looking quantity.** See
  the weight×polynomial entry below: 1e-95 and 1e19 are arresting numbers and
  they were not the mechanism. Confirm the proposed mechanism reproduces the
  observed failure before acting on it.
- **Comparing the QR output `U` against the raw basis `B`** to test a phase
  convention. Both come from the same routine, so a flip *inside* the basis is
  invisible — it reports a uniform `-1` that cancels in a matrix
  (`S·ana·S = ana`) and makes the problem look absent. Compare element by
  element instead.
- **Watching a scored ratio drift out of band** as the primary signal that
  something is wrong. It says only "something, somewhere, downstream of
  everything". Validate the isolated mathematics first.
- **Forming oscillator matrix elements as (weight) × (polynomial).** Correct
  conclusion, wrong reason on the first pass, and the wrong reason is the
  instructive part. It is *not* that ~1e-95 times ~1e19 loses digits —
  floating-point multiplication preserves relative accuracy, the summands are all
  comparable in size, and there is no cancellation to lose them to. The real
  failure is narrower: Golub–Welsch hands back the weight as the square of an
  eigenvector's first component, and once that component drops below the
  eigensolver's noise floor LAPACK returns it as **exactly zero**, deleting a
  whole node from the rule (4 of 60 nodes, 88 of 200). Use the DVR form, where
  the eigenvector entries *are* the product and that number is never requested.
  Any future rule that takes `√wᵢ` from an eigensolver inherits the trap.
