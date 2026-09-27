# Reading the paper

A map of S. Godfrey and N. Isgur, *Phys. Rev. D* **32**, 189 (1985),
[doi:10.1103/PhysRevD.32.189](https://doi.org/10.1103/PhysRevD.32.189), for
readers who want to check the code against the article. The article is not
distributed with this repository. Page numbers refer to the pages of the PDF.

## Where things are

| PDF pages | content |
|---|---|
| 1–2 | introduction and motivation |
| 2–4 | the Hamiltonian: relativized kinetic term, potentials, smearing, running coupling, numerical procedure |
| 4 | the three-stage diagonalization; mixing by tensor and antisymmetric spin-orbit terms |
| 5 | Table II: model parameters |
| 5–6 | light nonstrange and strange mesons (Fig. 4, with strange mixing angles) |
| 7 | isoscalar mesons, Fig. 5 |
| 8 | charmonium and charmed mesons (Figs. 6–7, with charmed mixing angles) |
| 9 | bottomonium and ``b``-flavored mesons (Figs. 8–9, with ``b``-flavored mixing angles) |
| 11 | Table III: isoscalar mixing |
| Sec. IV | strong decays (Tables IV–V), radiative transitions (Table VI), annihilation and electromagnetic properties (Table VII) |
| 35–38 | Appendix A: smearing and the effective potentials |
| 37–38 onward | Appendices B–D: wavefunction conventions, helicity amplitudes, electromagnetic operators |

## Appendix A in brief

Appendix A is where the relativization is made concrete. Its equations fall
into three groups:

1. **Smearing** (A7)–(A9): a Gaussian smearing function, the smeared
   potential as a convolution, and the mass-dependent width ``\sigma_{12}``.
2. **Smeared central potentials** (A10)–(A14): the pointwise ``G`` and ``S``,
   and closed forms for their smeared versions ``\tilde G`` and ``\tilde S``,
   with ``1/\tau_k^2 = 1/\gamma_k^2 + 1/\sigma_{12}^2``. The momentum factors
   that follow (A14) complete the central potential.
3. **Effective spin-dependent operators and the numerical method**
   (A15)–(A17): contact, tensor and spin-orbit operators built from derivatives
   of ``\tilde G`` and ``\tilde S`` with their momentum factors, and the
   factorization of matrix elements in an oscillator basis.

Equations (A1)–(A6) derive these operators from a scattering amplitude. They
motivate the construction but are not used directly as formulas.
[Equations to code](@ref) lists the implementation of each equation.

## Reading with care

- The article is typeset densely and some equations break across lines.
  Transcriptions of equations and table values in this repository were checked
  against the page images, not against text extraction.
- Spectrum values in Figs. 4–9 are read from level diagrams. The digitized
  values in `GIPaper/data/` carry an uncertainty of a few MeV.
- The paper works in the isospin limit: ``u`` and ``d`` have the same mass,
  and ``u``/``d`` splittings of heavy-light mesons are outside its scope.
- The mixing angles quoted in the captions of Figs. 4, 7 and 9 are not
  consistent with the rest of the calculation; see
  [Mixing angles: a paper erratum](@ref).
