# Quark-model strong decays: from width tables to resonance poles

Status: literature and project-scope note, 2026-09-19

## Question

How far have constituent-quark-model calculations been pushed when one starts
from model wave functions, supplies a strong-interaction transition operator,
and tries to understand branching fractions, total widths, and eventually
resonances whose widths or threshold couplings are too important to treat as a
small afterthought?

## Short answer

Very far, but in two methodologically different senses.

1. **Decay calculations are mature as spectroscopy tools.**  With a
   pseudoscalar-emission, flux-tube-breaking, microscopic, or especially a
   \({}^3P_0\) pair-creation operator, quark-model wave functions have been used
   to compute hundreds of two-body partial-wave amplitudes, total widths, and
   branching patterns across light mesons, strange mesons, heavy-light mesons,
   quarkonia, hybrids, and baryons.  These calculations are often useful for
   assignments, selection rules, mixing angles, and promising discovery
   channels.

2. **A width table is not yet a theory of a broad resonance.**  The more
   advanced projects reuse the same quark-model vertices inside self energies,
   coupled-channel Hamiltonians, or unitary scattering equations.  They solve
   for complex poles, residues, line shapes, continuum probabilities, and mass
   shifts.  This has been done convincingly in selected sectors—notably
   charmonium/open charm, positive-parity \(D_s\), and light scalars—but not as
   one universal, globally predictive extension of the quark model.

The practical frontier is therefore not simply “can the model predict a 300
MeV width?”  It can produce such a number.  The frontier is whether the same
calculation preserves unitarity and analyticity, treats all nearby channels and
bare levels coherently, and predicts the observable pole or amplitude rather
than assigning a Breit--Wigner width to an undressed eigenstate.

## A ladder of calculations

### 0. A bare spectrum

The Hamiltonian gives real eigenvalues \(M_a^{(0)}\) and normalizable valence
wave functions \(|a\rangle\).  This is the present GIModel strength.  Above a
strong threshold these are best called *bare basis states*, not automatically
physical resonances.

### 1. Born-level decay amplitudes

For a transition operator \(T\), one evaluates

\[
  \mathcal M^{LS}_{a\to i}(k)
  = \langle B C;k,LS|T|a\rangle
\]

and converts the on-shell value into \(\Gamma_{a\to i}\) using a chosen
phase-space and state-normalization convention.  Then

\[
  \Gamma_a=\sum_i\Gamma_{a\to i},\qquad
  \mathcal B_i=\Gamma_{a\to i}/\Gamma_a.
\]

This level can retain exact radial nodes, state-dependent length scales,
spin/flavor mixing, and channel-dependent partial waves.  It is the natural
next step beyond the current single-\(\beta\) strong-decay implementation.

It is most trustworthy when the parent is isolated, the amplitude and phase
space vary slowly across its width, the daughters are reasonably stable, and
final-state rescattering is weak.  Ratios, zeros, selection rules, and dominant
versus suppressed channels are usually more robust than absolute widths.

### 2. Dressed masses and perturbative unquenching

The same vertices generate an energy-dependent self-energy matrix,

\[
  \Sigma_{ab}(E)=\sum_i\int dk\,k^2
  \frac{\mathcal M_{a\to i}(k)\mathcal M^*_{b\to i}(k)}
       {E-E_{B_i}(k)-E_{C_i}(k)+i0}.
\]

Its real part shifts and mixes the bare levels; its imaginary part supplies
open-channel widths.  Solving a real mass equation while evaluating a width on
shell is useful, but it is still an approximation when \(\Sigma(E)\) changes
rapidly or several states mix.

Loop shifts that are almost common to a complete multiplet can be absorbed in
parameters already fitted in a quenched model.  Adding them explicitly without
refitting risks double counting.  This is a central lesson of the Geiger--Isgur
and Barnes--Swanson loop programs.

### 3. Coupled-channel poles

For broad, threshold-adjacent, or overlapping states, one should solve

\[
  \det[E\mathbf 1-M^{(0)}-\Sigma(E)]=0
\]

on the appropriate complex Riemann sheet, or equivalently solve a unitary
multichannel scattering equation.  The output is a pole
\(E_p=M_p-i\Gamma_p/2\), its residues into each channel, and the real-axis
line shapes.  At this level:

- the peak, Breit--Wigner mass, and pole mass need not agree;
- a threshold cusp can resemble a state;
- one bare state can generate more than one pole;
- a pole can be mostly continuum even if it descends from a quark-model seed;
- branching “fractions” for overlapping or very broad structures require an
  amplitude-level definition and are not always positive, process-independent
  probabilities.

### 4. Scattering and finite-volume validation

The strongest modern version supplements quark-model bare states and
pair-creation vertices with direct hadron--hadron interactions, fits or predicts
scattering amplitudes, and maps the Hamiltonian into finite volume for
comparison with lattice-QCD energy levels.  The positive-parity \(D_s\) work of
Yang et al. is an important example: quark-model states, pair creation, and
\(D^{(*)}K\) interactions are put into one Hamiltonian framework and confronted
with lattice levels.

## What counts as “the width is significant”?

There is no universal cutoff in \(\Gamma/M\).  A 100 MeV state can be well
described by a Born width if it is isolated and far from thresholds, while a
sub-MeV state at an S-wave threshold can require a nonperturbative treatment.
Escalate from a width table to a coupled-channel calculation when any of the
following is true:

- \(\Gamma\) is comparable to the separation between bare levels of the same
  quantum numbers;
- a relevant threshold lies within roughly the resonance region and the loop
  function changes rapidly there;
- an S-wave channel opens nearby;
- the predicted mass shift or off-diagonal loop mixing is comparable to the
  bare splitting;
- two resonances overlap and interfere in the same measured channel;
- a nominal two-body daughter is itself broad, making the physical process
  three-body;
- the result sought is a line shape, phase motion, pole residue, or
  compositeness rather than a rough discovery channel.

The diagnostic to compute is not just \(\Gamma/M\), but also
\(\Gamma/\Delta M\), distance to each threshold, partial-wave threshold power,
and the energy derivative of the self energy.

## Representative research programs

| Program | Sector and inputs | Decay dynamics | How far it goes | Main lesson |
|---|---|---|---|---|
| Godfrey--Isgur (1985) | Whole meson spectrum; analytic SU(6) oscillator states for the strong-decay table | Elementary pseudoscalar emission, \(g\,\boldsymbol\sigma\!\cdot\!\mathbf q+h\,\boldsymbol\sigma\!\cdot\!\mathbf p'\) | Partial-wave amplitudes and width estimates, including discussion of overlapping 300--500 MeV vectors and very broad scalars | Already recognized that broad states demand interference and coupled-channel amplitudes, but did not solve that dynamical problem |
| Kokoski--Isgur (1987) | Relativized meson-model wave functions | Flux-tube breaking with one elementary strength | Broad survey of ordinary-meson amplitudes | Directly demonstrates the “QM waves + one strong operator” research program; explains much of the empirical success of \({}^3P_0\) |
| Ackleh--Barnes--Swanson (1996) | Meson wave functions | Pair creation from scalar confinement plus one-gluon exchange | Microscopic comparison with \({}^3P_0\), including amplitude ratios | Scalar confinement gives roughly the needed strength; OGE is usually subdominant but can matter in selected channels |
| Capstick--Roberts | Relativized three-quark baryon waves | \({}^3P_0\) pair creation | \(N\pi\), quasi-two-body \(\Delta\pi,N\rho,N\eta,N\omega,\ldots\) widths and missing-baryon search channels | The program extends to baryons, but broad isobars and multibody final states make the jump to a resonance amplitude harder |
| Barnes--Black--Page (2003) | SHO strangeonium and kaonium waves | \({}^3P_0\) | 43 resonances, 525 modes, 891 amplitudes through about 2.2 GeV | Width tables can be exhaustive and useful, but the authors explicitly invoke a narrow-resonance approximation |
| Barnes--Godfrey--Swanson (2005) | GI and nonrelativistic charmonium spectra; SHO surrogates for open-charm decays | \({}^3P_0\) | All open-charm amplitudes for 40 states through the 4S region | An influential discovery/assignment atlas; also an example of “GI-based” not meaning the exact GI radial waves were used in the strong amplitudes |
| Godfrey--Moats (2015) and related sector surveys | Relativized heavy-light spectra and wave functions | \({}^3P_0\) | Masses, radiative widths, strong widths, branching patterns, and assignments over charm/charm-strange sectors | This is close to the desired wave-function-to-branching-fraction pipeline, still mostly at the isolated-resonance level |
| Segovia--Entem--Fernández (2013) | One constituent model for charmonium spectroscopy and transitions | Pair creation derived from the same scalar/vector interquark interactions used in the Hamiltonian | Open-charm widths compared with \({}^3P_0\) and data | A useful test of operator consistency; “microscopic” does not automatically outperform a fitted effective vertex |
| Cornell coupled-channel program | Bare charmonium plus open-charm channels | Cornell interaction creates/couples two-meson states | Mass shifts, mixing, widths, and eventually complex pole positions | A canonical example of crossing from decay calculation to resonance theory |
| Ortega--Segovia--Entem--Fernández | Constituent \(c\bar c\) or \(c\bar s\) cores plus meson--meson channels | \({}^3P_0\) couples Fock sectors; quark-model interactions act between hadrons | Bound states, resonance properties, component probabilities, and decay observables for \(X(3872)\) and positive-parity \(D_s\) | Near-threshold anomalies become mixed core--continuum states rather than failed mass predictions |
| Hao--Lu--Zou (2022), Ni--Wu--Zhong (2023) | Numerical potential-model waves in heavy-light sectors | \({}^3P_0\) or chiral-quark transition amplitudes inside coupled channels | Systematic spectra and widths, not just a single anomalous state | Heavy-light mesons are currently one of the most developed test beds for unified spectrum-plus-width calculations |
| Törnqvist / van Beveren--Rupp / related unitarized quark models | Light scalar or vector confinement seeds plus many two-meson channels | Flavor-related strong vertices with form factors, iterated to all orders | S-matrix poles, phase shifts, line shapes, pole doubling, dynamically generated companions | The framework can describe genuinely broad resonances, but conclusions become sensitive to channel completeness, chiral constraints, and parameterization |
| Yang--Wang--Wu--Oka--Zhu (2022) | Quark-model positive-parity \(D_s\) states | Pair creation plus direct \(D^{(*)}K\) interactions in Hamiltonian EFT | Infinite-volume structure and finite-volume lattice spectrum | A model can now be tested against lattice levels, not only masses and PDG widths |

## Sector-by-sector assessment

### Heavy quarkonium

This is the cleanest place to test decay vertices because the valence wave
functions are comparatively controlled.  Below open-flavor threshold the
potential model is excellent.  Above threshold, \({}^3P_0\) atlases predict
useful dominant modes and angular-momentum patterns, while Cornell-style and
other coupled-channel calculations can produce mass shifts and poles.  The
hard cases are S-wave thresholds and the \(XYZ\) region, where a compact core,
molecular channel, and production mechanism can all matter.

### Heavy-light mesons

This is arguably the best development sector for GIModel.  It offers:

- calculable initial and daughter wave functions;
- heavy-quark-symmetry relations and mixing-angle checks;
- narrow control states such as \(D_{s2}^*(2573)\);
- threshold-sensitive \(0^+\) and \(1^+\) states;
- broad nonstrange partners;
- direct links to lattice finite-volume spectra.

The \(D_{s0}^*(2317)\), \(D_{s1}(2460)\), \(D_{s1}(2536)\), and
\(D_{s2}^*(2573)\) form an unusually compact laboratory in which one can test
bare masses, heavy-quark mixing, S- versus D-wave couplings, threshold mass
shifts, continuum admixture, and strong widths together.

### Light and strange mesons

These sectors provide the most severe test and the largest historical decay
tables.  They also expose every limitation at once: relativistic constituents,
chiral dynamics, many open modes, broad daughter resonances, overlapping
radial/orbital levels, and important final-state interactions.  The excited
\(1^{--}\) \(\rho/\omega/\phi/K^*\) family discussed already by GI is a good
eventual stress test.  It should not be the first implementation target.

The light \(0^{++}\) sector shows the endpoint of unitarization: pole doubling,
large continuum components, Adler zeros, and line shapes that are not captured
by a sum of independent Breit--Wigners.  Here the phrase “the quark-model
state's branching fractions” may itself cease to be the right observable.

### Baryons

Quark-model decay operators have been used extensively to predict pion and
quasi-two-body branching fractions and guide missing-resonance searches.  The
extra internal coordinate, dense spectrum, spin-flavor recoupling, unstable
isobars, and dominant three-body final states make full unitarization more
expensive.  It is an important comparison literature, not the natural first
extension of this meson codebase.

## What the literature says about predictive power

The successful content is mostly *structural*:

- allowed and forbidden channels;
- relative partial waves and amplitude ratios;
- node-induced suppressions;
- heavy-quark spin multiplet relations;
- mixing-angle sensitivity;
- identification of channels in which a missing state should be visible;
- threshold-driven mass shifts and core/continuum mixing once channels are
  iterated.

Absolute widths have less controlled uncertainties.  The main sources are the
pair-creation strength and its flavor/scale dependence, wave-function shape,
relativistic normalization and phase space, physical versus model masses,
missing channels, unstable daughters, and final-state interactions.  A global
fit showing that the nominal \({}^3P_0\) strength must run with the reduced mass
is evidence that a single universal phenomenological constant is too rigid,
not evidence that the entire framework is useless.

The most defensible hierarchy of claims is:

1. selection rule or amplitude zero;
2. partial-wave and branching-ratio pattern;
3. order of magnitude of an isolated state's total width;
4. precise absolute width;
5. pole and line shape of a broad/threshold state.

Each step down the list requires more dynamics and more experimental input.

## Recommended project for GIModel

### Flagship sector: positive-parity charmed-strange mesons

Build a reusable strong-decay vertex on the solver-native GI wave functions,
validate it on isolated heavy-light states, and then embed exactly the same
vertex in a minimal coupled-channel calculation of the \(D_s\) \(0^+\) and
\(1^+\) sectors.

Why this sector:

- It starts from capabilities already present: unequal-mass GI wave functions,
  physical spin mixing, and strong-decay kinematics.
- It has both safe Born-level benchmarks and nonperturbative threshold cases.
- The first channels are small in number: \(DK\) for \(0^+\), and
  \(D^*K\) (with relevant higher channels added in a controlled order) for
  \(1^+\).
- It admits comparison with established coupled-channel quark-model work and
  with finite-volume lattice results.
- It creates infrastructure reusable for charmonium/open charm and eventually
  broad light vectors.

### Work packages

#### WP1 — extend the existing matrix elements to exact-wave Born amplitudes

- Audit the current decomposition first. `matrix_element(a)` is the existing
  Table IV/V reduced interaction factor \(cX(\bar q)\), while
  `decay_amplitude` completes it with the analytic single-\(\beta\) SHO spatial
  overlap. Correct its dimensional documentation and make that scope explicit.
- Preserve the existing implementation as the analytic SHO reference. Choose
  one solver-native extension first: apply the GI Eq. (19)
  pseudoscalar-emission operator to calculated waves, or introduce a general
  \({}^3P_0\) pair-creation vertex, retaining both as independent comparisons
  where their channel coverage overlaps.
- Define normalization, recoil momentum routing, mock/relativistic phase space,
  flavor conventions, and partial-wave projection in a written convention
  ledger before fitting anything.
- Reuse `radial_overlap`, `momentum_overlap`, and
  `physical_transition_amplitude` to evaluate the spatial integrals using the
  actual initial and daughter radial waves. Treat mixed states coherently at
  amplitude level.
- Reproduce analytic SHO limits as unit tests, then compare exact GI waves to
  equal-\(\beta\) and rms-matched SHO surrogates.
- Fit the minimal vertex strength only to a declared calibration set; reserve
  other modes for validation.

Deliverable: partial-wave amplitudes, widths, branching ratios, and a sensitivity
budget for a compact set of established \(D\) and \(D_s\) decays.

#### WP2 — self energies

- Reuse the momentum-dependent WP1 vertices under the loop integral.
- Compute diagonal shifts, off-diagonal \({}^1P_1\)--\({}^3P_1\) mixing, wave
  function renormalization, and continuum probabilities.
- Demonstrate regulator and basis convergence, and separate genuine threshold
  structure from a form-factor artifact.
- Refit the small affected sector, or explicitly define subtractions, rather
  than adding loop shifts to parameters that already absorb them.

Deliverable: dressed masses and compositions for the four lowest P-wave
\(D_s\) states, with controlled switches for each channel.

#### WP3 — poles and observables

- Continue the self energy to the required sheets and solve for complex poles.
- Construct the unitary real-axis \(DK\) and \(D^*K\) amplitudes.
- Report poles and residues first; quote Breit--Wigner parameters only when the
  line shape supports that language.
- Map the Hamiltonian to finite volume and compare its level pattern with the
  published lattice spectra.

Deliverable: a small but complete spectrum-to-pole calculation, rather than a
large table of uniterated widths.

#### WP4 — broad-resonance stress test

After the heavy-light machinery is stable, apply it to the excited isovector
\(1^{--}\) sector.  Include at least the nearby \(2S\) and \(1D\) bare seeds
coherently and the dominant \(\pi\pi\), \(\omega\pi\), and relevant multipion
effective channels.  The research question is whether the GI expectation of
overlapping broad \(\rho_S\) and \(\rho_D\) structures survives as a pole and
phase-shift statement.

## Acceptance criteria

A serious result should satisfy all of the following:

- SHO analytic limits and angular recoupling identities pass automatically.
- Widths are invariant under equivalent normalization conventions after the
  documented conversion.
- Closed channels give real self energies; open channels reproduce the optical
  theorem.
- Pole positions are stable against numerical contour/grid choices.
- Results include variation under wave-function route (HO versus FD), masses,
  pair-creation scale, and channel truncation.
- Calibration and validation channels are separated.
- Bare mass, dressed real-axis mass, Breit--Wigner parameter, and pole position
  are never used as synonyms.
- Comparisons use amplitudes or pole information where data analyses provide
  them, not only PDG summary widths.

## Selected primary literature

### Operators and decay surveys

- S. Godfrey and N. Isgur, *Mesons in a relativized quark model with
  chromodynamics*, Phys. Rev. D 32, 189 (1985),
  [DOI](https://doi.org/10.1103/PhysRevD.32.189).
- R. Kokoski and N. Isgur, *Meson decays by flux-tube breaking*, Phys. Rev. D
  35, 907 (1987), [DOI](https://doi.org/10.1103/PhysRevD.35.907).
- E. S. Ackleh, T. Barnes, and E. S. Swanson, *On the mechanism of open-flavor
  strong decays*, Phys. Rev. D 54, 6811 (1996),
  [arXiv](https://arxiv.org/abs/hep-ph/9604355).
- T. J. Burns, *Angular momentum coefficients for meson strong decay and
  unquenched quark models* (2014),
  [arXiv](https://arxiv.org/abs/1403.7538).
- T. Barnes, N. Black, and P. R. Page, *Strong decays of strange quarkonia*,
  Phys. Rev. D 68, 054014 (2003),
  [arXiv](https://arxiv.org/abs/nucl-th/0208072).
- T. Barnes, S. Godfrey, and E. S. Swanson, *Higher charmonia*, Phys. Rev. D
  72, 054026 (2005), [arXiv](https://arxiv.org/abs/hep-ph/0505002).
- S. Godfrey and K. Moats, *Properties of excited charm and charm-strange
  mesons*, Phys. Rev. D 93, 034035 (2016),
  [arXiv](https://arxiv.org/abs/1510.08305).
- J. Segovia, D. R. Entem, and F. Fernández, *Strong charmonium decays in a
  microscopic model*, Nucl. Phys. A 915, 125 (2013),
  [arXiv](https://arxiv.org/abs/1301.2592).
- J. Segovia, D. R. Entem, and F. Fernández, *Scaling of the \({}^3P_0\)
  strength in heavy meson strong decays*, Phys. Lett. B 715, 322 (2012),
  [arXiv](https://arxiv.org/abs/1205.2215).
- S. Capstick and W. Roberts, *Quark models of baryon masses and decays*, Prog.
  Part. Nucl. Phys. 45, S241 (2000),
  [arXiv](https://arxiv.org/abs/nucl-th/0008028).

### Unquenching and coupled channels

- P. Geiger and N. Isgur, *When can hadronic loops scuttle the
  Okubo--Zweig--Iizuka rule?*, Phys. Rev. D 47, 5050 (1993),
  [DOI](https://doi.org/10.1103/PhysRevD.47.5050).
- T. Barnes and E. S. Swanson, *Hadron loops: General theorems and application
  to charmonium*, Phys. Rev. C 77, 055206 (2008),
  [arXiv](https://arxiv.org/abs/0711.2080).
- E. J. Eichten, K. Lane, and C. Quigg, *Charmonium levels near threshold and
  the narrow state X(3872)*, Phys. Rev. D 69, 094019 (2004),
  [arXiv](https://arxiv.org/abs/hep-ph/0401210).
- E. J. Eichten, K. Lane, and C. Quigg, *New states above charm threshold*,
  Phys. Rev. D 73, 014014 (2006),
  [arXiv](https://arxiv.org/abs/hep-ph/0511179).
- P. G. Ortega, J. Segovia, D. R. Entem, and F. Fernández, *Coupled channel
  approach to the structure of the X(3872)*, Phys. Rev. D 81, 054023 (2010),
  [arXiv](https://arxiv.org/abs/0907.3997).
- P. G. Ortega et al., *Molecular components in P-wave charmed-strange mesons*,
  Phys. Rev. D 94, 074037 (2016),
  [arXiv](https://arxiv.org/abs/1603.07000).
- W. Hao, Y. Lu, and B.-S. Zou, *Coupled channel effects for the charmed-strange
  mesons* (2022), [arXiv](https://arxiv.org/abs/2208.10915).
- R.-H. Ni, J.-J. Wu, and X.-H. Zhong, *Unified unquenched quark model for
  heavy-light mesons with chiral dynamics* (2023),
  [arXiv](https://arxiv.org/abs/2312.04765).
- Z. Yang, G.-J. Wang, J.-J. Wu, M. Oka, and S.-L. Zhu, *Novel coupled channel
  framework connecting the quark model and lattice QCD for the near-threshold
  \(D_s\) states*, Phys. Rev. Lett. 128, 112001 (2022),
  [DOI](https://doi.org/10.1103/PhysRevLett.128.112001).

### Broad light resonances

- N. A. Törnqvist, *Understanding the scalar meson \(q\bar q\) nonet*, Z. Phys.
  C 68, 647 (1995), [arXiv](https://arxiv.org/abs/hep-ph/9504372).
- E. van Beveren and G. Rupp, *On the dominance of non-exotic meson-meson
  scattering by s-channel \(q\bar q\) confinement states and the classification
  of the scalar mesons* (2001),
  [arXiv](https://arxiv.org/abs/hep-ph/0110156).
- Z.-Y. Zhou and Z. Xiao, *The origin of light \(0^+\) scalar resonances*
  (2010), [arXiv](https://arxiv.org/abs/1007.2072).
- G. Rupp, S. Coito, and E. van Beveren, *Unquenching the meson spectrum: a
  model study of excited \(\rho\) resonances* (2016),
  [arXiv](https://arxiv.org/abs/1605.04260).

## Bottom line

The literature supports a staged program, not a single leap.  Exact GI wave
functions plus a general strong-decay operator would already be a meaningful
advance over the repository's current single-oscillator decay table.  But the
moment that result is used to discuss threshold states, large mass shifts, or
overlapping broad levels, the vertex must be iterated in a unitary
coupled-channel framework.  The positive-parity \(D_s\) sector is the natural
place to build and validate that bridge; excited light vectors are where to
find out how far it really reaches.
