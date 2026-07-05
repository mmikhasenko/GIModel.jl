# Report: A Computational Reproduction of the Godfrey–Isgur Model

`gi_reproduction.qmd` is one of the three first-class deliverables of this
repository, next to the **GIModel** (computation) and **GIPaper** (paper
comparison) packages. It is the high-level, code-free account of the project:
it walks the physics of the original paper section by section — the
relativistic Hamiltonian, the Appendix-A smearing and momentum-dependent
operator sandwiches, the running coupling, the spin-dependent operators, and
the isoscalar annihilation machinery — teaching the quantum-mechanical
computations while demonstrating, sector by sector, that the published spectrum
is reproducible from the paper's own ingredients.

The rendered `gi_reproduction.pdf` is tracked so readers never need the
toolchain.

## Building

```bash
quarto render report/gi_reproduction.qmd
```

(Requires Quarto with a LaTeX toolchain; `report/*.tex` build artifacts are
git-ignored, the PDF is tracked.)

## Regenerating figures

The curated figures under `figures/` are copies of generated output:

- Sector ladder plots (`charmonium.png`, `isovector.png`, …):
  `julia GIPaper/scripts/plot_all_sector_spectra.jl`, output in
  `GIPaper/scripts/spectrum_plots/`.
- Wavefunction/basis figures (`wavefn_*.png`):
  `julia GIPaper/scripts/plot_basis_wavefunctions.jl` and
  `julia GIPaper/scripts/plot_offdiagonal_blowup.jl`, same output directory.

Copy refreshed plots into `figures/` deliberately — the report should only
change when the model or the comparison meaningfully changes.

## Numbers

Every quantitative claim in the report traces to a generated artifact:
residual tables come from `GIPaper/scripts/run_all_spectrum_checks.jl`
(reports under `GIPaper/docs/residual_reports/`), diagnostics from the
`GIPaper/scripts/audit_*.jl` scripts. Run `bash scripts/verify_project.sh`
to regenerate the full evidence chain before editing report numbers.

## Scope

The current edition covers Secs. II–III of the original paper (model and
spectroscopy). The coupling analysis of Sec. IV (strong and electromagnetic
transitions) is the natural second installment as the Table V–VII work
matures.
