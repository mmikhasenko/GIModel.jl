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

## Still required for the Phase-IV gate

- General orbital-gradient integrals beyond S-to-S.
- End-to-end `matrix_element` composition, including an explicit and tested
  rule for which final pseudoscalar is the elementary emitted field.
- Appendix-C to relativistic-normalization conversion.
- Native solved HO/FD convergence and node-sensitive certificates.
- Separate calibration and validation channel sets.

The emitted-field role is intentionally not guessed from a display label. It
must be made invariant under `TwoMesonChannel` canonicalization before the
public end-to-end method is enabled.
