# Godfrey-Isgur Reproduction

Local, reproducible reproduction of the Godfrey-Isgur relativized quark model
for meson masses.

The working plan lives in `docs/orchestrator_task.md`. The seed database in
`data/seed/` is only a bootstrap from later sources and must be verified against
the original 1985 paper before being treated as reference data.

## Immediate Workflow

1. Add the primary PDF as `paper/Godfrey-Isgur-1985.pdf`.
2. Run seed validation:

   ```bash
   python3 scripts/validate_seed.py
   ```

3. Generate text/page extractions from the PDF.
4. Create raw extraction CSVs in `data/raw/`.
5. Promote verified rows into `data/clean/` with provenance preserved.

## Authority Rules

- The original 1985 paper is the primary authority.
- Later quoted tables are seed data only.
- Every numerical value must carry source, table/page when available,
  extraction method, and confidence.
- Extraction data and cleaned physics data stay separate.
- Discrepancies are classified before code or data changes are made.

