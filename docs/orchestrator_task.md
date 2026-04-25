# Orchestrator task: reproduce and verify the Godfrey--Isgur meson spectrum

## Goal

Build a local, reproducible implementation of the Godfrey--Isgur relativized quark model for meson masses, with a verified machine-readable database of reference masses and parameters.

Primary target:

- S. Godfrey and N. Isgur, *Mesons in a relativized quark model with chromodynamics*, Phys. Rev. D 32, 189--231 (1985), DOI: 10.1103/PhysRevD.32.189.

The attached `godfrey_isgur_seed_masses.csv` is only a bootstrap table assembled from later papers quoting GI values. It is not the final authority.

## Non-negotiable rules

1. Do not silently refit parameters.
2. Do not treat later quoted GI tables as primary truth.
3. Every numerical reference value must carry provenance:
   - source paper,
   - table number if available,
   - page number if available,
   - extraction method,
   - confidence flag.
4. Keep raw extraction and cleaned physics data separate.
5. When a discrepancy is found, classify it before changing code:
   - extraction error,
   - convention error,
   - solver/numerics error,
   - physics-model implementation error,
   - later-source mismatch with original paper.
6. No “looks close enough” without a tolerance policy.

## Repository layout

```text
godfrey-isgur-reproduction/
  paper/
    Godfrey-Isgur-1985.pdf
    pages/
      p001.png
      ...
  data/
    seed/
      godfrey_isgur_seed_masses.csv
      godfrey_isgur_sources.csv
    raw/
      tables_from_pdf_text.csv
      tables_from_page_images.csv
      parameters_from_pdf.csv
    clean/
      masses.csv
      parameters.toml
      mixings.csv
      extraction_audit.csv
  src/
    GIModel/
      Hamiltonian implementation
      basis / solver code
  test/
    test_extracted_tables.jl
    test_spin_independent.jl
    test_ccbar.jl
    test_bbbar.jl
    test_open_charm.jl
    test_open_bottom.jl
  docs/
    conventions.md
    formula_map.md
    extraction_notes.md
    residual_reports/
```

## Phase 0: inventory and source control

Deliverables:

- `docs/source_inventory.md`
- `data/seed/godfrey_isgur_seed_masses.csv`
- `data/seed/godfrey_isgur_sources.csv`

Tasks:

1. Add the original GI PDF to `paper/`.
2. Record DOI, page range, and bibliographic metadata.
3. Record every secondary source used for seed values.
4. Add the current seed table.
5. Run schema validation on the seed table:
   - no missing `sector`,
   - no missing `assignment`,
   - all masses numeric,
   - all masses in MeV,
   - all rows have `source_status`.

Verification gate:

- Seed database loads without schema errors.
- Duplicate check reports only intentional duplicates such as `unassigned` rows with different assignments.

## Phase 1: PDF extraction

Deliverables:

- `paper/text/pdftotext_layout.txt`
- `paper/text/pdftotext_bbox.html`
- `paper/pages/*.png`
- `data/raw/tables_from_pdf_text.csv`
- `data/raw/tables_from_page_images.csv`
- `docs/extraction_notes.md`

Recommended commands:

```bash
pdftotext -layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_layout.txt
pdftotext -bbox-layout paper/Godfrey-Isgur-1985.pdf paper/text/pdftotext_bbox.html
pdftoppm -r 300 -png paper/Godfrey-Isgur-1985.pdf paper/pages/gi
```

Tasks:

1. Identify all pages containing numerical tables.
2. Extract tables from `pdftotext -layout` first.
3. Extract/check the same tables from page images independently.
4. Compare text-based and image-based extraction row by row.
5. Any mismatch becomes an entry in `extraction_audit.csv`.

Required fields for raw rows:

```csv
source,page,table,row_index,raw_state_label,raw_assignment,raw_mass_model,raw_mass_exp,raw_notes,extraction_method,confidence
```

Verification gate:

- All table pages have been rendered as images.
- Every table row has a page and table label.
- All disagreements between text and image extraction are explicitly listed.

## Phase 2: clean physics schema

Deliverables:

- `data/clean/masses.csv`
- `data/clean/mixings.csv`
- `data/clean/parameters.toml`
- `docs/conventions.md`

Suggested mass schema:

```csv
source,page,table,sector,quark1,quark2,state_label,n,L,S,J,P,C,mass_model_MeV,mass_exp_MeV,assignment,source_status,confidence,notes
```

Tasks:

1. Normalize sectors:
   - `nnbar`, `ssbar`, `nsbar`,
   - `ccbar`, `bbbar`, `bcbar`,
   - `cqbar`, `csbar`, `bqbar`, `bsbar`.
2. Preserve original labels.
3. Add parsed quantum numbers where unambiguous.
4. For mixed heavy-light states, keep both:
   - original paper label,
   - cleaned label.
5. Do not force uncertain rows into a wrong schema. Use `confidence = low` and explain.

Verification gate:

- Every cleaned row points to a raw row.
- Every cleaned mass has units MeV.
- Ambiguous mixed states are not silently assigned to pure singlet/triplet states.

## Phase 3: parameter and formula map

Deliverables:

- `data/clean/parameters.toml`
- `docs/formula_map.md`
- `docs/conventions.md`

Tasks:

1. Extract all published parameters:
   - quark masses,
   - string tension,
   - constant offset,
   - strong-coupling prescription,
   - smearing parameters,
   - spin-dependent interaction parameters.
2. Map every formula in the implementation to an equation or paragraph in the paper.
3. Write explicit conventions:
   - spectroscopic notation,
   - spin operator normalization,
   - tensor operator,
   - spin-orbit decomposition,
   - units,
   - basis normalization,
   - radial wave function convention.

Verification gate:

- No physics formula appears in code without an entry in `formula_map.md`.
- Parameter file is read by code; parameters are not hardcoded.

## Phase 4: minimal solver

Deliverables:

- spin-independent solver,
- test script for heavy quarkonia,
- first residual plots/tables.

Start with sectors:

1. `bbbar`
2. `ccbar`
3. `bcbar` if available
4. open charm/bottom
5. light sectors last

Tasks:

1. Implement the semirelativistic kinetic term.
2. Implement the spin-independent confinement + Coulomb-like potential.
3. Use a numerically stable basis:
   - harmonic oscillator basis, Gaussian expansion, or momentum-space discretization.
4. Compare level spacings before absolute masses.
5. Fix numerical convergence before adding spin splittings.

Verification gate:

- Eigenvalues stable under basis-size increase.
- Heavy-quarkonium spin-averaged levels agree with the reference table within a documented tolerance.

## Phase 5: spin-dependent structure

Deliverables:

- spin-spin term,
- tensor term,
- spin-orbit terms,
- mixing treatment,
- sector tests.

Tasks:

1. Add one interaction at a time.
2. Use diagnostic splittings:
   - `J/psi - eta_c`,
   - `Upsilon - eta_b`,
   - `chi_c2 - chi_c1 - chi_c0`,
   - `chi_b2 - chi_b1 - chi_b0`.
3. Implement mixing only after pure states are correct.
4. For each bad splitting, classify whether the issue is:
   - spin algebra,
   - radial matrix element,
   - smearing,
   - relativization factor,
   - table extraction.

Verification gate:

- `ccbar` and `bbbar` low-lying S/P/D states pass sector-specific tolerances.
- Mixing states are covered by explicit tests.

## Phase 6: sector expansion

Suggested order:

```text
bbbar -> ccbar -> bcbar -> bsbar -> bqbar -> csbar -> cqbar -> ssbar/nsbar/nnbar
```

For each sector produce:

- `docs/residual_reports/<sector>.md`
- comparison table:
  - paper mass,
  - reproduced mass,
  - residual,
  - tolerance,
  - pass/fail,
  - likely issue if failed.

Tolerance policy:

- heavy-heavy low states: start with 5--10 MeV target after full implementation;
- heavy-light: 10--30 MeV depending on state and mixing;
- light mesons: use broader tolerance and report qualitative failures separately.

Verification gate:

- No sector marked complete unless all rows are either pass or have a classified failure.

## Phase 7: final audit

Deliverables:

- `docs/final_reproduction_report.md`
- `data/clean/masses.csv`
- `data/clean/parameters.toml`
- reproducible scripts/notebooks
- CI tests

Final checks:

1. Re-run all extraction tests.
2. Re-run all model tests.
3. Confirm no accidental parameter refits.
4. Confirm all secondary-source seed values have been replaced or confirmed against the original paper.
5. Produce a final table of:
   - exact matches,
   - small residuals,
   - large residuals,
   - unresolved extraction ambiguities,
   - known physics limitations.

## Agent roles

Recommended multi-agent split:

- **Orchestrator**: plans phases, enforces gates, writes progress reports.
- **Extractor**: parses PDF text/images and creates raw tables.
- **Physics reviewer**: checks Hamiltonian, conventions, and quantum numbers.
- **Numerics implementer**: builds solver and convergence tests.
- **Verifier**: compares reproduced values against tables and opens issues.

The orchestrator should not let implementers patch the model before the verifier classifies discrepancies.

## First milestone

A useful first milestone is not the full paper. It is:

```text
parameter table + ccbar + bbbar + spin-independent solver + first spin splittings
```

Acceptance criteria:

- Original GI PDF table values for `ccbar` and `bbbar` are digitized and checked.
- `parameters.toml` exists.
- The solver reproduces the gross `1S`, `2S`, `1P`, `1D` structure.
- Residuals are tabulated.
- Every discrepancy above tolerance has a proposed diagnosis.
