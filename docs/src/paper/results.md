# Results at a glance

This page summarizes each comparison with the 1985 paper. The numbers come
from the recorded reports, which list every row; follow the links for detail.

## Spectra

All seven flavor sectors of the paper's figures, 209 states in total, are
reproduced within a few MeV on average:

| sector | states | mean \|Δ\| (MeV) | max \|Δ\| (MeV) |
|---|---:|---:|---:|
| isoscalar | 48 | 2.6 | 5.5 |
| charmonium | 28 | 2.9 | 8.6 |
| ``b``-flavored | 21 | 3.2 | 6.9 |
| isovector | 30 | 3.3 | 8.0 |
| bottomonium | 30 | 3.5 | 7.7 |
| strange | 30 | 4.8 | 14.4 |
| charmed | 22 | 6.1 | 24.8 |

Δ is the model mass minus the paper's value, read from its figures, computed
with the paper's oscillator method. The
figures are printed as level diagrams, so the reference values carry a
digitization uncertainty of a few MeV, comparable to the mean residuals. The
largest residuals are in states whose mass depends on singlet–triplet mixing
(see below). Source: [spectrum scorecard](reports/scorecard.md).

## Same-J mixing angles

The mixing angles quoted in the 1985 figure captions are **not** reproduced,
while the masses are. The model angles agree instead with the same authors'
later calculations with identical parameters. The evidence points to an error
in the 1985 caption values, and the implementation is not tuned toward them.
[Mixing angles: a paper erratum](@ref) gives the full account; the numbers are
in the [mixing-angle report](reports/mixing_angles.md).

## Isoscalar mixing (Table III)

The paper gives two pseudoscalar prescriptions, Eq. (18a) (`P1`) and
Eq. (18b) (`P2`). Evaluated on the solved wavefunctions, `P2` reproduces the
Table III flavor–radial compositions of ``\eta``, ``\eta'`` and their radial
excitations with a mean component RMS of 0.015 and a mean mass difference of
11 MeV. `P1` agrees for the ``\eta`` but less well for the ``\eta'`` and the
radial excitations (mean RMS 0.08, mean mass difference 56 MeV). The vector
(``\omega``–``\phi``) and tensor (``f_2``–``f_2'``) blocks of Eq. (16) match
Table III to its printed precision.

The isoscalar spectrum in the table above uses a `P1` block calibrated to the
paper's pseudoscalar masses; all other
isoscalar channels use the literal Eq. (16) or ideal mixing. Source:
[Table III report](reports/table_iii_mixing_audit.md).

## Strong decays (Tables IV–V)

Table V was computed in the paper with a single Gaussian wavefunction
(``\beta = 0.40`` GeV) for all mesons and two fitted strengths, ``A`` and
``S_0``. Refitting ``A`` and ``S_0`` to the same two decays gives
``A = 1.644`` and ``S_0 = 3.291`` (paper: 1.67 and 3.27). With them, every
Table V amplitude that has a paper value is reproduced, most to within a few
percent.

One detail is essential to reproduce the table: the numerical columns use the
leading constant ``S_0`` for the structure-dependent amplitudes (``S``, ``D``,
``P`` classes), not the full ``q``-dependent polynomial printed in Table IV.
Kinematics use current PDG masses, which accounts for most of the small
differences. Source: [Table V report](reports/table_v_reproduction.md).

This reproduction deliberately uses the paper's simplification. Evaluating
the same operator on the solved wavefunctions goes beyond it; see
[Strong decays beyond the paper](@ref).

## Radiative transitions (Table VI)

All 79 rows (42 M1, 35 E1, 2 M2) are computed from the solved wavefunctions
with the paper's operators. For the rows with a definite paper value, the
median ratio of model to paper amplitude is 1.00 for M1 and 1.01 for E1, with
the sign reproduced in 58 of 62 cases. The large deviations are concentrated
in transitions from radially excited ``\eta`` states, whose amplitudes involve
strong cancellations and depend on the isoscalar mixing. Source:
[Table VI report](reports/table_vi_photon_decays.md).

## Annihilation and electromagnetic properties (Table VII)

| quantity | rows | median model/paper | signs |
|---|---:|---:|---|
| gluonic widths | 16 | 1.05 | 16/16 |
| leptonic decay constants | 24 | 1.03 | 24/24 |
| two-photon widths, pure flavor | 8 | 1.06 | 8/8 |
| two-photon widths, mixed isoscalars | 3 | 0.78 | 3/3 |
| charge radii | 2 | 1.02 | — |

The only fitted constant is the charge-radius exponent ``f = 0.2``, which the
paper fits to the ``\pi^+``. Source: [Table VII report](reports/table_vii_annihilation_em.md).

## Reading these numbers

The reproduction is not a fit. The spectrum uses the Table II parameters and
the Appendix A operators as written; the decay tables add only the constants
the paper itself fits. Remaining differences come from three sources:
digitization of printed figures, updated experimental kinematics, and, for the
mixing angles, the published values themselves.
