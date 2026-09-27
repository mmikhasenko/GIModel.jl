# Conventions and units

This page collects the conventions used throughout the package. Signs of
amplitudes and mixing angles depend on them, so check here before comparing a
number with another calculation.

## Units

| quantity | unit |
|---|---|
| masses, energies, momenta | GeV |
| distances | ``\mathrm{GeV}^{-1}`` (``1\ \mathrm{GeV}^{-1} = 0.1973`` fm) |
| ``\langle r^2\rangle`` | ``\mathrm{GeV}^{-2}`` |
| ``\langle p^2\rangle`` | ``\mathrm{GeV}^2`` |
| decay widths from [`decay_width`](@ref) | MeV |
| M1 transition moments | nuclear magnetons |
| E1 and M2 photon amplitudes | ``\mathrm{MeV}^{1/2}`` |
| two-photon amplitudes | ``\mathrm{GeV}^{1/2}`` |

Field names carry their unit when it is not obvious: `mass_GeV`,
`central_GeV`, `contact_shift_GeV`, `beta_GeV`. The parameter file stores
quark masses and the constant ``c`` in MeV (`m_c_MeV`, `c_MeV`); the loader
converts them to GeV.

## State labels

A level is written ``n\,{}^{2S+1}L_J`` and labeled in code as
`"n^(2S+1)L_J"`, for example `"1^3P_2"`:

- `n ≥ 1` counts radial excitations, starting at 1 in every partial wave (so
  the lowest P wave is `1P`, not `2P`);
- the multiplicity is 1 for spin singlets and 3 for spin triplets;
- `L` is the letter `S`, `P`, `D`, `F`, `G` for ``L = 0, \dots, 4``.

Parity and charge conjugation follow from the quantum numbers:
``P = (-1)^{L+1}`` and, for self-conjugate mesons, ``C = (-1)^{L+S}``.

After mixing, a physical state keeps the label of the unmixed level it is
assigned to, by ascending mass. The label identifies the state; its
composition is given by [`physical_components`](@ref).

## Quark order

`Meson(mq, f1, f2)` is the ``f_1\bar f_2`` meson: `f1` is the quark and `f2`
the antiquark. `Meson(mq, :c, :u)` is ``c\bar u`` (a ``D^0``) and
`Meson(mq, :u, :s)` is ``u\bar s`` (a ``K^+``).

The masses are symmetric, but one sign is not: the antisymmetric spin-orbit
operator contains ``(1/m_1^2 - 1/m_2^2)`` and changes sign when quark and
antiquark are exchanged. Mixing angles of ``c\bar u`` and ``u\bar c`` are
therefore opposite.

## Mixing angles

For unequal-mass mesons, the two ``J = L`` states mix. GIModel writes the
lower physical state as

```math
|\text{low}\rangle = \cos\theta\,|n\,{}^1L_L\rangle + \sin\theta\,|n\,{}^3L_L\rangle ,
```

with the spin-orbit angular matrix element
``\langle {}^1L_L|\,\boldsymbol L\cdot(\boldsymbol S_1 - \boldsymbol S_2)\,|{}^3L_L\rangle = +\sqrt{L(L+1)}``,
where ``\boldsymbol S_1`` is the spin of the quark (`f1`). This matches the
convention in the captions of Figs. 4, 7 and 9 of the paper. The
coefficients are read directly from [`physical_components`](@ref) or from
`state.mixings[1].components`.

With this convention the model gives positive angles for ``u\bar s`` (light
quark, heavier antiquark) and negative angles for ``c\bar u``, ``c\bar s``,
``b\bar u`` and ``b\bar s`` (heavy quark, lighter antiquark), as in the paper.

## Phases of wavefunctions

- **Radial waves.** Every radial wave is normalized with its outermost lobe
  positive. A radially excited wave therefore changes sign at the origin once
  per node. The large-radius convention matches the one stated after
  Eqs. (B35)–(B36) of the paper.
- **Orbital states.** Spherical harmonics follow the Condon–Shortley
  convention.
- **Angular-momentum coupling.** States couple orbital angular momentum first,
  then spin: ``\langle L m_L; S m_S | J M\rangle``.
- **Spin and flavor wavefunctions.** Singlet and triplet spin states follow
  Eqs. (B26)–(B29), and flavor states follow Eqs. (B1)–(B15), including the
  signs of charged and heavy-light mesons.
- **Mixed states.** An eigenstate of an annihilation block is chosen to have
  a positive overlap with the basis state it is assigned to, so a small
  admixture with a negative sign does not flip the whole state. Eigenstates of
  spectroscopic (spin-orbit and tensor) blocks keep the sign returned by the
  diagonalization; only relative signs within a state are meaningful there.

An overall phase of an external state changes the sign of an amplitude but
never a width. Relative signs between interfering amplitudes are fixed by the
conventions above.

## Flavor states

The light quarks are degenerate, so the model does not distinguish ``u`` from
``d`` dynamically. The flavor label matters in two places:

- **Charges.** Photon emission and leptonic currents need explicit `:u` or
  `:d`; the averaged `:q` is rejected there.
- **Isospin states.** A neutral pion is ``(u\bar u - d\bar d)/\sqrt2``. Build it
  from `:u` and `:d` spectra as a [`PhysicalState`](@ref) with the two
  components (see [Strong decays beyond the paper](@ref)).

For isoscalar annihilation mixing, the nonstrange channel is `Meson(mq, :q, :q)`
and stands for ``(u\bar u + d\bar d)/\sqrt2``.

## Pseudoscalar emission

The operator of Eq. (19) uses the phase conventions of Appendix C of the paper:
helicity amplitudes with the emitted meson along ``+z``, converted to partial
waves as ``A_L = \sqrt{(2L+1)/(2J+1)} \sum_m \langle L\,0; 1\,m|J\,m\rangle H_m``.
For identical daughters, exchange symmetry ``(-1)^{L+S}`` removes forbidden
waves and the allowed ones carry a factor ``1/\sqrt2``. The width is

```math
\Gamma = \frac{q}{2\pi(2J_i + 1)} \sum_{L,S} |A_{LS}|^2
```

in GeV, returned in MeV by [`decay_width`](@ref).

All mass-mixing eigenvectors use positive overlap with their assigned unmixed
state, ordered by unmixed mass. If that overlap is zero, the largest component
sets the phase. This convention is shared by spin-orbit, tensor, and annihilation
mixing. `diagonalize_mixing_block(...; phase_anchor=1)` explicitly requests a
common first-row anchor when comparison with another convention requires it.
