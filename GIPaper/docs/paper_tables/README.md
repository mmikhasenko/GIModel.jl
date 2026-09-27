# Paper tables and figures

Every numbered table and spectrum figure of Godfrey & Isgur,
*Phys. Rev. D* **32**, 189 (1985), recomputed with the current packages and
laid out row-for-row as printed. Each object has one driver and one generated
file; the original GI value is kept beside each computed value where the paper
prints one.

Regenerate after the canonical reports (for example after
`bash scripts/verify_project.sh`), then check the results:

```sh
julia GIPaper/scripts/paper_tables/generate.jl
julia GIPaper/scripts/paper_tables/check.jl
```

| Paper object | Driver | Output |
|---|---|---|
| Table I — importance of confinement | [`table_i.jl`](../../scripts/paper_tables/table_i.jl) | [`table_i.tsv`](table_i.tsv) |
| Table II — model parameters | [`table_ii.jl`](../../scripts/paper_tables/table_ii.jl) | [`table_ii.tsv`](table_ii.tsv) |
| Table III — isoscalar mixings | [`audit_table_iii_mixings.jl`](../../scripts/audit_table_iii_mixings.jl) | [`table_iii_mixing_audit.md`](../residual_reports/table_iii_mixing_audit.md) |
| Table IV — reduced partial-wave amplitudes | [`table_iv.jl`](../../scripts/paper_tables/table_iv.jl) | [`table_iv.tsv`](table_iv.tsv) |
| Table V — strong decays | [`reproduce_table_v.jl`](../../scripts/reproduce_table_v.jl) | [`table_v_reproduction.md`](../residual_reports/table_v_reproduction.md) |
| Table VI — radiative decays | [`audit_table_vi_photon_decays.jl`](../../scripts/audit_table_vi_photon_decays.jl) | [`table_vi_photon_decays.md`](../residual_reports/table_vi_photon_decays.md) |
| Table VII — leptonic, annihilation, radii | [`audit_table_vii.jl`](../../scripts/audit_table_vii.jl) | [`table_vii_annihilation_em.md`](../residual_reports/table_vii_annihilation_em.md) |
| Fig. 1 — annihilation coupling scales | [`fig_i.jl`](../../scripts/paper_tables/fig_i.jl) | [`fig_i.tsv`](fig_i.tsv) |
| Fig. 2 — running coupling | [`fig_ii.jl`](../../scripts/paper_tables/fig_ii.jl) | [`fig_ii.tsv`](fig_ii.tsv) |
| Figs. 3–9 — sector spectra | [`spectra.jl`](../../scripts/paper_tables/spectra.jl) | [`fig_iii.tsv`](fig_iii.tsv) … [`fig_ix.tsv`](fig_ix.tsv) |
| Figs. 3, 4, 6 captions — tensor-mixed states | [`caption_mixing.jl`](../../scripts/paper_tables/caption_mixing.jl) | [`caption_mixing.tsv`](caption_mixing.tsv) |
| Figs. 4, 7, 9 captions — same-J mixing angles | [`audit_mixing_angles.jl`](../../scripts/audit_mixing_angles.jl) | [`mixing_angles.md`](../residual_reports/mixing_angles.md) |

Two registries support the tables: [`mass_inputs.toml`](mass_inputs.toml), the
PDG 2026 kinematic masses with their historical counterparts
([`mass_inputs.jl`](../../scripts/paper_tables/mass_inputs.jl)), and
[`quantity_columns.toml`](quantity_columns.toml), the role, unit and method of
each extra numerical column in Tables III and V–VII
([`quantity_columns.jl`](../../scripts/paper_tables/quantity_columns.jl)).

Notes on individual objects:

- Table I: GI do not define their rounded confinement percentage; the driver
  reports the Hellmann–Feynman share of the 2S–1S splitting.
- Table II: values are read from `GIModel.default_parameters_path()`. Λ has no
  active value: the model uses the fixed Gaussian running coupling it motivates
  (see [model inputs](../../../docs/model_inputs.md)).
- Table IV: `A` and `S0` are the Table V calibration amplitudes, read from its
  report.
- Figs. 3–9: same-J singlet–triplet pairs are mixed with the spin-orbit
  off-diagonal element; the `mixed` column marks them.
