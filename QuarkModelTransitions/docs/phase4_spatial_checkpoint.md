# Phase IV spatial checkpoint

Status: first vertical slice complete, 2026-09-23. Phase IV remains active.

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

## Still required for the Phase-IV gate

- General orbital-gradient integrals beyond S-to-S.
- End-to-end `matrix_element` composition, including an explicit and tested
  rule for which final pseudoscalar is the elementary emitted field.
- Appendix-C to relativistic-normalization conversion.
- Native solved HO/FD convergence and node-sensitive certificates beyond the
  completed S-wave slice.
- Separate calibration and validation channel sets.

The emitted-field role is intentionally not guessed from a display label. It
must be made invariant under `TwoMesonChannel` canonicalization before the
public end-to-end method is enabled.
