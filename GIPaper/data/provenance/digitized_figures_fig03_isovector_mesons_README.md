> Archived transcription provenance. Paths and tool descriptions below refer to the original extraction; reference values are preserved.

# Figure 3: Isovector Mesons

Source: Godfrey-Isgur 1985, Fig. 3, PDF page 6.

Screenshot:

- `gi_spectrum_page-06.png` (render of the 1985 paper PDF)

Scope:

- Transcribes the model-state labels printed next to solid mass bars.
- Records spectroscopic composition and plotted mass in GeV.
- Excludes axis tick labels and experimental shaded bands.
- Includes a clean comparison replot generated from the CSV:
  `figure_03_replot.svg` and `figure_03_replot.png`.

Status: first-pass manual visual transcription. The labels are mostly clear, but
the high-mass D/G-wave region should receive one visual audit pass before
promotion to clean data.

Rebuild the comparison plot with:

```bash
python3 extraction/plot_figure3_digitization.py
rsvg-convert data/raw/digitized_figures/fig03_isovector_mesons/figure_03_replot.svg \
  -o data/raw/digitized_figures/fig03_isovector_mesons/figure_03_replot.png
```

Caption note:

- Sector: isovector mesons, `[-u dbar, 1/sqrt(2)(u ubar - d dbar), d ubar]`.
- Caption reports significant mixing in the `1--` sector:
  `1--(1.45) ~= 1.00(2^3S_1) + 0.04(1^3D_1)`.
