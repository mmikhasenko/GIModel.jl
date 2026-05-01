# Research note: GI spin-independent central operator

Date: 2026-05-01

## Question

The current solver is stuck at the spin-independent part of the Godfrey-Isgur
model. The heavy-quarkonium diagnostics show that most of the remaining error is
a common mass offset, while radial/orbital spacings are already fairly close.
This points at the central operator and its Appendix-A relativization, not at
spin-orbit/tensor tuning.

## Local state

- Active default path: pointwise `G(r) + S(r)` in `src/GIModel.jl`.
- Experimental paths:
  - `appendix_a_smearing`: direct numerical 3D convolution of both `G` and `S`.
  - `coulomb_1d`: 1D Gaussian on `G` only.
  - `appendix_a_derivative_g`: first derivative-expansion proxy,
    `G + laplacian(G)/(4 sigma^2)`, with pointwise `S`.
- The repo already records that these are diagnostic comparators, not the full
  GI Appendix-A operator.

The largest local clue is `docs/residual_reports/heavy_quarkonium_diagnostics.md`:

- charmonium mean residual: about +76 MeV; after removing common offset, mean
  absolute residual about 16 MeV.
- bottomonium mean residual: about +38 MeV; after removing common offset, mean
  absolute residual about 9 MeV.

That means the next profitable target is a paper-faithful spin-independent
central path before changing fine-structure scale factors.

## Literature facts

### Original GI model

Godfrey and Isgur, Phys. Rev. D 32, 189 (1985), is the primary authority. OSTI's
record identifies the key ingredients as a universal one-gluon-exchange plus
linear-confinement potential and says relativistic effects are crucial.

In the paper text already extracted locally, Appendix A says the raw on-shell
potentials are not used directly. The model introduces:

1. smearing/nonlocality;
2. momentum-dependent effective interactions;
3. harmonic-oscillator basis diagonalization for operators that mix p and r.

The paper also says the confinement term is assumed unmodified by relativistic
corrections, motivated by the one-dimensional QED analogy. In practice this
means:

- smear `S(r)` geometrically to get `tilde S(r)`;
- do not multiply the central `S` by additional momentum-dependent factors;
- central `G` gets both smearing and the momentum-dependent sandwich.

### Later GI-family papers make the smeared central formulas explicit

Open-access modified-GI papers quote the standard GI smearing formulas in a
cleaner form than the local OCR. The important equations are:

```text
alpha_s(r) = sum_k (2 alpha_k / sqrt(pi)) int_0^(gamma_k r) exp(-x^2) dx

G(r) = -4 alpha_s(r) / (3 r)
S(r) = b r + c

sigma = sqrt(
    s^2 * (2 m1 m2 / (m1 + m2))^2
    + sigma0^2 * (1/2 + 1/2 * (4 m1 m2 / (m1 + m2)^2)^4)
)

tilde G(r) =
    - sum_k [8 alpha_k / (3 sqrt(pi) r)] int_0^(tau_k r) exp(-x^2) dx
  = - sum_k [4 alpha_k / (3 r)] erf(tau_k r)

tau_k = 1 / sqrt(1/sigma^2 + 1/gamma_k^2)

tilde S(r) =
    b r * [
        exp(-sigma^2 r^2)/(sqrt(pi) sigma r)
        + (1 + 1/(2 sigma^2 r^2)) erf(sigma r)
    ] + c
```

These formulas match the parameters and alpha-grid already in this repo:
`alpha_k = (0.25, 0.15, 0.20)` and
`gamma_k = (0.5, sqrt(10)/2, sqrt(1000)/2)`.

The later papers then describe the central confinement part as

```text
tilde H_conf = G'(r) + tilde S(r)

G'(r) =
    (1 + p^2/(E1 E2))^(1/2) tilde G(r)
    (1 + p^2/(E1 E2))^(1/2)

Ei = sqrt(mi^2 + p^2)
```

The spin-dependent pieces get analogous `(m_i m_j / E_i E_j)^(1/2 + epsilon_i)`
sandwiches, but that is downstream of fixing the central operator.

## Expanded web research beyond the 1985 paper

The web search changes the plan in one important way: we should treat later
GI/MGI application papers as implementation documentation, not just citations.
They repeatedly reproduce the same relativization recipe in cleaner notation
than the original scan.

## User-supplied deep research synthesis

The external report at
`/Users/mikhailmikhasenko/Downloads/deep-research-report.md` adds a broader
post-1985 development map. Its main value for the present stuck point is
prioritization:

- The most reproducible and highest-value target is still the **bare
  single-channel GI Hamiltonian**: masses, wavefunctions, mixing angles, and
  radiative matrix elements.
- No official public 1985 GI reference implementation was found in the surveyed
  public sources, so a transparent reimplementation with explicit conventions is
  the right path.
- The likely convention traps are exactly the ones this repo is encountering:
  operator ordering, smearing convention, basis truncation, phase conventions,
  and decay/loop form factors.
- Later screened potentials and explicit coupled-channel/unquenched models are
  real post-GI developments, but they should not be mixed into the current
  reproduction until the bare Hamiltonian is stable.
- A practical reproduction sequence is:
  original GI 1985 single-channel spectra; then Godfrey-Kokoski 1991
  heavy-light P-wave/mixing tests; then Barnes-Godfrey-Swanson 2005 higher
  charmonium observables; then screened-potential variants; then explicit
  coupled channels.

This reinforces the immediate implementation target:

1. finish the original spin-independent relativized central operator;
2. verify heavy quarkonia first, where threshold dressing is least intrusive;
3. add unequal-mass mixing only after the equal-mass central/fine-structure base
   is trustworthy;
4. treat screening and coupled channels as later, separately flagged physics.

The report also gives a useful expectation scale: later Godfrey-style
applications caution that uncoupled relativized predictions should not be
expected to beat roughly the 10-20 MeV level in general. That is a sensible
success criterion for the low-lying bare Hamiltonian once the common central
offset is removed.

### Most useful references

1. **Pang et al. 2017, strange mesons, EPJ C 77, 861**
   - Best implementation reference found.
   - Gives the nonrelativistic pieces, running alpha, smearing width, closed-form
     `tilde G`, closed-form `tilde S`, momentum-dependent sandwiches, total
     Hamiltonian, and the SHO-basis solution strategy.
   - Crucially states the spin-independent relativized confinement Hamiltonian:
     `tilde H_12^conf = G'(r) + tilde S(r)`.

2. **Wang et al. 2018, higher bottomonium, EPJ C 78, 915**
   - Confirms the same high-level GI structure in a heavy-quarkonium setting:
     relativistic kinetic energy, `G(r)`, `S(r)`, Gaussian smearing, and the
     central Coulomb momentum sandwich.
   - Useful because our current diagnostics are `ccbar`/`bbbar`.

3. **Yang et al. 2023, charmed-strange mesons with coupled-channel effects,
   EPJ C 83, 1098**
   - Confirms the two-step relativization:
     first smear `G` and `S`; then apply momentum-dependent factors.
   - Helpful independent check that the central Coulomb factor is
     `sqrt(1 + p^2/(E1 E2))` on both sides of `tilde G`.

4. **Godfrey/Moats/Godfrey-Swanson follow-up spectroscopy papers**
   - These papers use the relativized quark model as production machinery and
     give spectra/wave functions, but usually point back to GI for details.
   - They are useful for validation targets and confidence that the original
     parameter set is still the baseline, but less useful than the MGI papers for
     reconstructing the operator.

5. **Review/status papers**
   - Reviews describe the GI model as relativistic kinetic energy, running
     coupling, flavor-dependent smearing, and factors replacing quark masses by
     quark energies.
   - They support the interpretation but do not add implementation details.

### Cross-source consensus

Across the accessible web sources, the central spin-independent recipe is:

```text
H = sqrt(p^2 + m1^2) + sqrt(p^2 + m2^2) + tilde H_conf + ...

tilde H_conf = G'(r) + tilde S(r)

G'(r) = A(p) tilde G(r) A(p)
A(p) = sqrt(1 + p^2/(E1 E2))
Ei = sqrt(mi^2 + p^2)
```

and the smearing is:

```text
tilde f(r) = integral d^3r' rho(r-r') f(r')
rho = sigma^3 / pi^(3/2) exp[-sigma^2 (r-r')^2]
```

For the GI running Coulomb ansatz, that convolution has the closed form
`gamma_k -> tau_k`. For the linear potential, it has the closed form `tilde S`
listed above.

### Important caution from the web sources

Most MGI papers later replace the linear potential by a screened potential for
high excitations:

```text
br -> b(1 - exp(-mu r))/mu
```

That is **not** part of the original GI reproduction target. It is useful later
if we study high excitations/coupled-channel effective behavior, but it should
not be mixed into the current stuck point. For now, keep `S(r)=br+c` and fix the
original GI relativization first.

### What external references decide

The expanded web evidence strongly favors this sequence:

1. implement closed-form `tilde G` and `tilde S`;
2. implement central `G' = A tilde G A`;
3. use the existing FD `p^2` eigenbasis as a practical substitute for the SHO
   machinery;
4. compare with MGI/SHO papers only after the operator is in place.

This is better grounded than trying to infer coefficients from the original
paper scan alone.

## What this means for this codebase

The current `appendix_a_smearing` branch numerically convolves the pointwise
`G(r)` and `S(r)` over a finite radial mesh. That is useful as a check, but it is
not the best implementation target. For `G`, the closed form above is safer:
smearing the running-coupling Coulomb term simply replaces each `gamma_k` by
`tau_k`. For `S`, the closed form avoids the finite-domain tail and small-r
instability of the direct convolution.

The real missing piece is not merely smearing. It is the central momentum
sandwich on `tilde G`:

```text
A(p) tilde G(r) A(p),  A(p) = (1 + p^2/(E1 E2))^(1/2)
```

This operator is nonlocal in coordinate space. GI used a harmonic-oscillator
basis because this lets them evaluate `f(p) g(r) f(p)` cleanly: `f(p)` in
momentum space and `g(r)` in coordinate space.

## Recommended implementation path

### Step 1: add closed-form smeared central functions

Add a new named central mode, e.g. `:appendix_a_closed_form`, with:

- `smeared_coulomb_G_closed(params, m1, m2, r)`;
- `smeared_confinement_S_closed(params, m1, m2, r)`;
- careful `r -> 0` limits.

Small-r limits:

```text
tilde G(0) = - sum_k 4 alpha_k tau_k * 2/(3 sqrt(pi))
tilde S(0) = 2 b/(sqrt(pi) sigma) + c
```

This should replace the naive finite-grid convolution as the first serious
spin-independent milestone.

### Step 2: compare closed form against existing direct convolution

For `G` and `S` separately:

- sample the same `ccbar` grid used in `scripts/compare_central_paths.jl`;
- compare closed form vs `smear_3d_radial` with a much larger tail;
- ignore the first couple finite-difference points if necessary, but require
  agreement away from the origin.

This verifies normalization, `sigma` units, and the `tau_k` convention.

### Step 3: test spectra with only `tilde G + tilde S`

Run the heavy-quarkonium diagnostics using the closed-form smeared central
potential but without the momentum sandwich. Expected outcome:

- short-distance Coulomb singularity softens;
- S-wave centers should move most;
- P/D/F states should move less;
- if the common offset worsens, that confirms the momentum sandwich is not
  optional.

This is a bracket, not the final model.

### Step 4: implement `A(p) tilde G(r) A(p)` on the existing FD basis

The current FD Hamiltonian already diagonalizes `p^2` in `sqrt_kinetic_matrix`.
Use the same eigenbasis to build:

```text
A = U diag(sqrt(1 + lambda/(E1(lambda) E2(lambda)))) U'
V_G_eff = A * Diagonal(tilde_G_values) * A
H = kinetic + V_G_eff + Diagonal(tilde_S_values)
```

where `lambda` are eigenvalues of the radial `p^2` operator for the same `L`.

This is not exactly the paper's HO-basis implementation, but it is the closest
operator-level analogue inside the existing solver. It directly targets the
central missing physics while keeping the current FD architecture.

### Step 5: only then revisit HO basis

If the FD momentum sandwich gives unstable or basis-sensitive results, port the
HO basis. But do not make that the first move: the closed-form `tilde G`,
`tilde S`, and FD sandwich can be implemented and tested faster, and they will
teach whether the stuck point is really the central Appendix-A operator.

## Sources

- Godfrey and Isgur, "Mesons in a relativized quark model with chromodynamics,"
  Phys. Rev. D 32, 189 (1985), DOI: https://doi.org/10.1103/PhysRevD.32.189
- OSTI record for the original paper:
  https://www.osti.gov/biblio/5736888
- Pang et al., "A systematic study of mass spectra and strong decay of strange
  mesons," Eur. Phys. J. C 77, 861 (2017):
  https://link.springer.com/article/10.1140/epjc/s10052-017-5434-0
- Wang et al., "Higher bottomonium zoo," Eur. Phys. J. C 78, 915 (2018):
  https://link.springer.com/article/10.1140/epjc/s10052-018-6372-1
- Yang et al., "The mass spectrum and strong decay properties of the
  charmed-strange mesons within Godfrey-Isgur model considering the
  coupled-channel effects," Eur. Phys. J. C 83, 1098 (2023):
  https://link.springer.com/article/10.1140/epjc/s10052-023-12275-3
- Godfrey and Moats, "Bottomonium mesons and strategies for their observation,"
  Phys. Rev. D 92, 054034 (2015):
  https://doi.org/10.1103/PhysRevD.92.054034
- Godfrey and Moats, "Properties of excited charm and charm-strange mesons,"
  Phys. Rev. D 93, 034035 (2016):
  https://doi.org/10.1103/PhysRevD.93.034035
- Barnes, Godfrey, and Swanson, "Higher charmonia," Phys. Rev. D 72, 054026
  (2005): https://doi.org/10.1103/PhysRevD.72.054026
