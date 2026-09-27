# Mixing angles: a paper erratum

The captions of Figs. 4, 7 and 9 of the 1985 paper quote thirteen
singlet–triplet mixing angles for strange, charmed and ``b``-flavored mesons.
GIModel reproduces the masses in those figures but not most of the angles. This
page explains why we conclude that the published angles, not the
implementation, are at fault.

## The comparison

In the paper's convention (lower state
``\cos\theta\,|{}^1L_L\rangle + \sin\theta\,|{}^3L_L\rangle``; see
[Conventions and units](@ref)), the ground-state P-wave angles are:

| system | GIModel (FD) | GI 1985 caption | Godfrey–Kokoski 1991 |
|---|---:|---:|---:|
| ``u\bar s`` | +4.3° | +34° | +5° |
| ``c\bar u`` | −25.8° | −41° | −26° |
| ``c\bar s`` | −39.6° | −44° | −38° |
| ``b\bar u`` | −27.9° | −43° | −31° |
| ``b\bar s`` | −41.1° | −45° | −40° |
| ``b\bar c`` | +69.1° | −53° | +68° |

The GIModel column is from the recorded report, computed with the
finite-difference solver; the oscillator solver gives values within 2° of it
(see [Heavy-light mesons and mixing](@ref)). Most higher states (``1D``,
``2P``, ``1F``, ``2D``, ``1G``) show the same pattern; the [mixing-angle report](reports/mixing_angles.md) lists all thirteen
caption angles.

## Why the implementation is not at fault

**1. A real bug was found and fixed first.** The first version of the code
applied the smeared contact interaction only to S waves. Appendix A and
Eqs. (23)–(26) of the paper keep it in every partial wave, because smearing
gives it a finite range. After the fix, the singlet–triplet gaps of all eight
comparable two-state systems matched the paper and the median angle
difference fell from 33° to 15°. What remained was concentrated in the
antisymmetric spin-orbit term.

**2. The same authors' later calculations agree with GIModel.** The
Godfrey–Isgur model with the 1985 parameters was used again in later papers.
Their angles agree with GIModel within a few degrees, including the sign of the
``b\bar c`` angle:

- S. Godfrey and R. Kokoski, *Phys. Rev. D* **43**, 1679 (1991), Table I:
  ``K``, ``D``, ``D_s``, ``B``, ``B_s``, ``B_c``. An early draft of this paper
  (July 1986) already contains these angles, together with the same
  contact, spin-orbit and tensor expectation values that GIModel gives, and
  presents them as the 1985 calculation without comment.
- R. Blundell, S. Godfrey and B. Phelps, *Phys. Rev. D* **53**, 3712 (1996):
  ``\theta_K`` consistent with the Godfrey–Kokoski value.
- S. Godfrey, *Phys. Rev. D* **70**, 054017 (2004): ``B_c``.
- S. Godfrey and K. Moats, *Phys. Rev. D* **93**, 034035 (2016): ``D`` and ``D_s``.
- S. Godfrey, K. Moats and E. S. Swanson, *Phys. Rev. D* **94**, 054025 (2016):
  ``B`` and ``B_s``.

**3. The published masses are consistent with the model.** A P-wave multiplet
of an unequal-mass meson is described by four spin-dependent strengths:
symmetric spin-orbit, tensor, contact, and antisymmetric spin-orbit. Extracting
them from the paper's own printed masses and angles, the first three agree
with GIModel in all eleven systems. Only the antisymmetric strength implied by
the caption angles differs, by factors of 4 to 18. The masses barely depend on
that strength, because the diagonal gap dominates; the angles depend on it
directly.

## The conclusion

The captions of Figs. 4, 7 and 9 appear to contain angles that the 1985
calculation itself did not produce. Where they came from is unknown. GIModel
is **not** tuned toward them: its parameters and operators are the paper's, and
the reproduction reports keep both the 1985 captions and the later values as
comparison targets, never replacing one with the other.

When you compare with the 1985 paper, compare mixing angles with the later
Godfrey–Isgur values as well as with the captions.

## Still open

The small ``{}^3S_1``–``{}^3D_1`` admixtures in the ``\psi(3770)`` region,
which the paper also discusses, are a separate question and have not been
resolved by this investigation.
