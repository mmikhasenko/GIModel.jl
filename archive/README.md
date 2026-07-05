# Archive

Concluded material from the mass-spectrum reproduction phase. Nothing here is
needed for day-to-day work with **GIModel** (computation) or **GIPaper**
(paper comparison); it is kept for provenance and can be revived if a question
reopens. Everything was moved here with `git mv`, so full history is available
via `git log --follow`.

## poster/

Chaptered write-up of the reproduction workflow (introduction, inputs, radial
solver, Appendix-A smearing, spin-dependent operators, results), with images
and theme. Superseded as documentation by `docs/` and the residual reports;
kept as a self-contained narrative of the project as of mid-2026.

## report/

Quarto/LaTeX report (`gi_reproduction.qmd` → `.pdf`) with curated sector
figures. A snapshot of the reproduction status; regenerate figures with
`GIPaper/scripts/plot_all_sector_spectra.jl` if a new edition is ever needed.

## docs/

- `deep-research-report.md` — early literature/context survey used to plan the
  reproduction.
- `extraction_notes.md` — running notes from digitizing the paper's figures
  and tables. The promoted results live in `GIPaper/data/` with provenance
  columns; consult this only when auditing an extraction decision.
- `source_inventory.md` — inventory of paper sources/OCR artifacts taken
  during extraction. Current authority rules are in the top-level `README.md`.

## table_iii_forensics/

Forensic investigation of the Table III isoscalar mixing amplitudes, concluded
when the annihilation implementation was validated (the surviving checks are
`GIPaper/scripts/audit_table_iii_mixings.jl` and the GIPaper test suite):

- `audit_annihilation_method.jl` + `annihilation_method_audit.md` — diagnosis
  of the FD-vs-HO wavefunction-at-origin scale in the annihilation matrix
  elements (led to the HO wave cache used by the current code).
- `investigate_table_iii_suspects.jl` + `table_iii_suspect_investigation.md` —
  row-by-row investigation of suspect digitized Table III entries.
- `infer_table_iii_mass_matrix.jl` + `table_iii_inverse_matrix.md` — inverse
  reconstruction of the paper's mass matrix from published amplitudes.

The scripts still run: they activate the GIPaper environment at the repository
root and write their reports next to themselves in this folder.
