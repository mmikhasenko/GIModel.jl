# Digitized Figures

One folder per figure. These files transcribe plotted labels from the spectrum
figures: spectroscopic composition and predicted model mass.

Rules:

- Keep each figure in its own folder.
- Record PDF page, figure number, sector, and screenshot provenance.
- Transcribe model-state labels separately from experimental shaded bands.
- Mark uncertain labels explicitly.
- Promote only audited rows to `data/clean/`.

Current folders:

- `fig03_isovector_mesons/`: Fig. 3 labels and clean replot.
- `fig04_strange_mesons/`: Fig. 4 labels and clean replot.
- `fig05_isoscalar_mesons/`: Fig. 5 labels and clean replot.
- `fig06_charmonia/`: Fig. 6 labels and clean replot.
- `fig07_charmed_mesons/`: Fig. 7 panel labels and clean replot.
- `fig08_bottomonia/`: Fig. 8 labels and clean replot.
- `fig09_b_flavored_mesons/`: Fig. 9 panel labels and clean replot.

Use `python3 scripts/plot_spectrum_digitizations.py` to regenerate the Fig.
4-9 SVG/PNG comparison plots from their CSV files.
