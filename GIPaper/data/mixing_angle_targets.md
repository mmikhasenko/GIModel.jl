# Same-J mixing-angle comparison targets

`mixing_angle_targets.csv` retains all 13 Godfrey–Isgur (1985) caption angles
and adds the six ground-P angles in Godfrey–Kokoski (1991), Table I. These are
comparison inputs only: they do not enter GIModel or alter the 1985 authority
policy. No original number is replaced.

## Sources and scope

- `GI1985`: *Phys. Rev. D* **32**, 189, Figs. 4, 7, 9. The caption values are
  the same ones previously embedded in `scripts/audit_mixing_angles.jl`.
- `GK1991`: Godfrey–Kokoski, *Phys. Rev. D* **43**, 1679, Table I, transcribed
  from the revised preprint (July 1986 draft) with the same model parameters;
  later publications corroborate it. Table I supplies six directly
  convention-mappable targets; later literature is not silently merged into
  those six numbers.

Both use the GI85 Table-II model parameters (b=0.18 GeV², m_u=0.220,
 m_s=0.419, m_c=1.628, m_b=4.977 GeV). Angles are in degrees. The
`printed_step_deg` field records the one-degree published resolution, not an
experimental uncertainty or an acceptance tolerance.

## Convention and provenance fields

`f1 f2bar` is the comparison ordering; `source_f1 source_f2bar` records the
printed ordering. `printed_angle_deg` is transcribed without conversion;
`angle_deg` is mapped to the comparison ordering. The lower state is
`cos(theta) |n 1L_L> + sin(theta) |n 3L_L>`.

Only GK91 K changes sign: its s ubar angle -5° becomes +5° for u sbar.
The conversion is recorded as `constituent_exchange`; all other rows use
`identity`. State ordering is kept fixed. The projected singlet probability
`cos(theta)^2` is reported alongside angles, avoiding basis-sign dependence.
The calculation also reports probability outside the selected radial pair.

## Resolution and regeneration

The implementation investigation is **closed**: missing higher-L contact was
fixed, and the remaining large caption mismatch is classified as a historical
reference discrepancy supported by the later GI-model calculations. Why the
1985 captions differ remains unknown. Small quantitative residuals remain
reported, without tuning toward either source.

Run `julia GIPaper/checks/audit_mixing_angles.jl` from the repository root.
It computes the same production spectra as before and writes both
[`mixing_angles.md`](../reports/mixing_angles.md) and
[`mixing_reference_comparison.csv`](../reports/mixing_reference_comparison.csv).
The seven GI85 rows without a GK91 counterpart are retained, not treated as
validated by missing data. Tensor ψ(3.82) admixtures and isoscalar annihilation
are separate investigations and are not closed by this result.
