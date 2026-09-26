# Phase IV spatial checkpoint

Status: general orbital evaluator and public Eq. (19) assembly complete,
2026-09-26. Phase IV validation remains active.

## What is implemented

- `PseudoscalarEmission(g, h, quark_masses)` is the public GI Eq. (19)
  operator. Constituent masses are explicit operator data; there is no hidden
  default parameter set.
- GIModel owns the representation-specific derivative overlap
  `integral (du_final/dr) u_initial f(r) dr` for both `OscillatorWave` and
  `MeshWave`. The transition package does not inspect native wave fields.
- The equal-mass S-wave spatial columns are evaluated as

  ```text
  direct:  g q <j0(alpha q r)>
  recoil:  h s <[d/dr - 1/r] j1(alpha q r)>
  ```

  where `alpha` is the spectator constituent-mass fraction and `s` is the
  quark/antiquark plane-wave sign.
- Real and complex off-shell momenta are supported by the overlap layer.
- The Phase-III shared coefficient matrix is multiplied by these numerical
  columns without introducing a per-channel formula.
- Arbitrary `L_i,m_i -> L_f,m_f` columns use one finite plane-wave multipole
  sum. The recoil gradient has exactly two generic branches,
  `L_f -> L_f+1` and `L_f -> L_f-1`; both reduce to combinations of
  `radial_overlap` and `radial_derivative_overlap`. There are no decay-class
  or channel-specific spatial methods.

## Independent checks now passing

For normalized equal-beta oscillator ground states,

```text
<j0(q r / 2)> = exp[-q^2/(16 beta^2)]
<[d/dr - 1/r] j1(q r / 2)> = -(q/4) exp[-q^2/(16 beta^2)].
```

Combining the derived `rho -> pi pi` coefficient row `[-1,-1,-1,+1]`
therefore gives

```text
-2 q exp[-q^2/(16 beta^2)] (g + h/4),
```

which recovers the Table-IV relation `A = (g+h/4) beta` up to the separately
tracked external normalization and overall phase. A mesh evaluation of the
same oscillator wave checks the independent centered-difference route.

The first solved-wave convergence certificate uses the nonstrange central
`1S`/`2S` pair at `q=0.45 GeV`, with the adaptive native HO solution compared
to FD grids on `r_max=28 GeV^-1`:

| column | levels final/initial | HO | FD 1200 | FD 1800 | HO/FD-1200 rel. |
|---|---:|---:|---:|---:|---:|
| direct | 1/1 | 0.2974798 | 0.2975040 | 0.2975035 | 0.0081% |
| recoil | 1/1 | 0.03187284 | 0.03187425 | 0.03187485 | 0.0044% |
| direct | 1/2 | -0.01890433 | -0.01886033 | -0.01886067 | 0.233% |
| recoil | 1/2 | 0.01699235 | 0.01700901 | 0.01700881 | 0.098% |
| recoil | 2/1 | -0.02104328 | -0.02104829 | -0.02104940 | 0.024% |
| direct | 2/2 | 0.2568387 | 0.2568111 | 0.2568090 | 0.011% |
| recoil | 2/2 | 0.02751844 | 0.02751378 | 0.02751449 | 0.017% |

The cancellation-prone `1S`--`2S` direct column sets the initial measured
regression threshold at `0.3%`; it is the only entry above `0.1%`. FD movement
from 1200 to 1800 points is already much smaller than the HO/FD difference, so
the comparison is not being limited by that last FD refinement. This
certificate covers the S-wave slice only and is not the final Phase-IV
orbital-sector certificate.

The general-orbital checks add four independent constraints:

- The analytic equal-beta Cartesian result
  `<P_0|exp(s i k z)|S> = s i k/(sqrt(2) beta) exp[-k^2/(4 beta^2)]` is
  reproduced for both plane-wave signs.
- At zero momentum the two gradient directions give the conjugate oscillator
  matrix elements `+/- i h beta/sqrt(2)`.
- All three spherical gradient components have equal magnitude, while invalid
  magnetic projections are exact zeros.
- The derived `A1 -> rho pi` S/D ratio agrees to `2e-12` with the independent
  Table-IV expressions for `A=(g+h/4) beta` and
  `S=[3h-(g+h/4)q^2/(2beta^2)] beta`. Both partial waves come from one
  helicity vector and the same shared integral columns.

For actual solved nonstrange `1P -> 1S` waves at `q=0.45 GeV`, native HO and
FD (`ngrid=1200`, `rmax=28 GeV^-1`) agree by `0.063%` for the direct column,
`0.039%` for the longitudinal recoil column, and `0.035%` for each transverse
recoil column. The regression gate is `0.1%` for this non-node-sensitive
orbital slice.

## Still required for the Phase-IV gate

- Native solved HO/FD convergence for additional orbital and node-sensitive
  sectors beyond the completed S-wave and `1P -> 1S` slices.
- Separate calibration and validation channel sets.

The public method now uses the ordered
`TwoMesonChannel(surviving, emitted)` contract and verifies that the second
state is `J^P=0^-`. It never infers the role from a display label.
The Appendix-C normalization and coherent physical-state composition are part
of the same `matrix_element` route; the Table-V backend is not a fallback.
