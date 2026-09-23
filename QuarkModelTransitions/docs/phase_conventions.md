# Transition phase conventions

This ledger is executable policy for `QuarkModelTransitions`. A convention is
not considered established until a sign-sensitive test names its source and
checks an independent identity or paper fixture.

## Established

- **Radial waves.** GIModel fixes the outermost radial lobe positive. Radial
  excitations therefore follow the large-radius sign convention stated below
  GI Eqs. (B35)--(B36).
- **Orbital states.** Spherical harmonics use the Condon--Shortley convention.
  The alternating sign of GI Eqs. (B31)--(B36) belongs to `Y_LM`; it is not an
  additional radial-wave phase.
- **Spectroscopic coupling.** States couple orbital then spin angular momentum,
  `CG(L,mL,S,mS,J,M)`, matching the `L dot S` order stated before Eq. (B37).
- **Constituent spin.** The singlet/triplet kets are GI Eqs. (B26)--(B29).
  Spherical Pauli components use `sigma_0 = sigma_z`,
  `sigma_+1 = -(sigma_x+i sigma_y)/sqrt(2)`, and
  `sigma_-1 = +(sigma_x-i sigma_y)/sqrt(2)`. Consequently
  `<chi_00|sigma_q,0|chi_10> = +1` while the antiquark value is `-1`.
- **Flavor states and transfers.** The sparse `q qbar` kets use GI
  Eqs. (B1)--(B15), including the charged-meson and heavy-light signs. Quark
  emission contracts the emitted-meson ket as
  `sqrt(2) conj(P[initial,final])`; antiquark emission uses
  `-sqrt(2) conj(P[final,initial])`, the anti-fundamental sign of Eq. (B25).
  Spectator flavor is matched exactly, so OZI-forbidden contractions are exact
  zeros rather than numerical cancellations.
- **Appendix-C helicities.** For meson emission along `+z`, `h_0 = H_0` and
  `h_m = sqrt(2) H_m` for positive `m`, as in Eq. (C3). The vector-pseudoscalar
  projection is

  `A_L = sqrt((2L+1)/(2J+1)) sum_m CG(L,0,1,m,J,m) H_m`,

  with the parity-related negative-helicity term absorbed into `h_m`. Tests
  reproduce every entry and sign of GI Table XI for parent `J = 0:5`.
- **Physical mixing.** A parent coefficient is used directly and both daughter
  coefficients are conjugated. One overall phase per external state changes
  the amplitude covariantly and leaves its width invariant.
- **Identical daughters.** Exchange contributes `(-1)^(L+S)`; forbidden odd
  exchange waves are removed and the allowed amplitude receives one
  `1/sqrt(2)` normalization.

## Pending before native Eq. (19)

- quark versus antiquark emission phases and momentum routing;
- topology ordering and the sign of each shared orbital-integral label;
- conversion from the GI Appendix-C normalization to the package's canonical
  relativistic two-body normalization.
