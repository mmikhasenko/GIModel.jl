# Paper extraction utilities

These Python programs create or promote source artifacts from the scanned
paper. They are deliberately outside `GIPaper/src`, `GIPaper/test`, and the
Julia-only `GIPaper/scripts` reproduction environment:

- `vision_ocr_paper.py` renders and transcribes scanned PDF pages.
- `plot_figure3_digitization.py` and `plot_spectrum_digitizations.py` rebuild
  visual checks for digitized spectrum labels.
- `data_checks.py promote-clean` writes audited rows into `data/clean/`.

They are provenance tooling, not dependencies of either Julia package and not
part of `Pkg.test()`.
