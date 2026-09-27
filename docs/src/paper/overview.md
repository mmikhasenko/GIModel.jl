# Reproducing Godfrey–Isgur

GIModel was built by reproducing S. Godfrey and N. Isgur, *Phys. Rev. D*
**32**, 189 (1985), table by table. This section describes that comparison:
how it is organized, what agrees, where it differs, and why.

## Two packages, one direction

The repository holds two Julia packages with a one-way dependency:

- **GIModel** (with its submodule `QuarkModelTransitions`) is the physics.
  It knows nothing about the paper's tables. You give it quark flavors and
  quantum numbers; it returns states and amplitudes.
- **GIPaper** (in `GIPaper/`) is the comparison. It holds the digitized paper
  data, maps each table row to a GIModel calculation, and writes the reports
  reproduced in this section.

GIPaper depends on GIModel, never the other way around. A reference row can
therefore never steer a calculation: the model decides which state is which,
and the comparison layer only assigns names.

## Rules of the comparison

- **The paper is the authority.** Every reference number is transcribed from
  the published article with its page and table; the raw extraction is kept
  separately from the cleaned data (`GIPaper/data/provenance/`). The article
  itself is not distributed with the repository.
- **No tuning toward the paper.** The model uses the paper's Table II
  parameters and the literal Appendix A operators. Where the reproduction
  differs, the difference is reported, not fitted away.
- **Fitted quantities are labeled.** Only the constants the paper itself fits
  are refitted, in the same way: the strong-decay strengths ``A`` and ``S_0``
  from the same two decays, and the charge-radius exponent ``f = 0.2`` from
  the ``\pi^+``.
- **Experimental inputs are pinned.** Decay kinematics use a fixed table of
  PDG 2026 masses (`GIPaper/data/mass_inputs/`), with no fallback to model
  masses.

## Running the comparison

From the repository root:

```bash
julia --project=GIPaper -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia GIPaper/checks/run_all_spectrum_checks.jl     # spectra, all sectors
julia GIPaper/checks/audit_mixing_angles.jl         # same-J mixing angles
julia GIPaper/checks/audit_table_iii_mixings.jl     # isoscalar mixing
julia GIPaper/checks/reproduce_table_v.jl           # strong decays
julia GIPaper/checks/audit_table_vi_photon_decays.jl
julia GIPaper/checks/audit_table_vii.jl             # annihilation, leptonic, radii
```

Each script writes a Markdown report to `GIPaper/reports/`, stamped with the
solver settings that produced it. The pages under
**Recorded reports** are copies of those files, made when this documentation was
built. `bash scripts/verify_project.sh` runs the whole chain together with both
test suites.

To use the comparison layer in your own code:

```julia
using GIModel, GIPaper
params, mq = load_parameters_and_quark_masses(GIPaper.model_parameters_path())
reference = load_reference_spectrum(reference_spectrum_path("charmonium"))
rows = compare_reference(params, mq, reference; solver = OscillatorSolver())
```

Each row pairs a paper value with the model value and their residual.

## What the paper computes, and where it is here

| paper | content | reproduced by | report |
|---|---|---|---|
| Figs. 4–9, Table II | meson spectra in all flavor sectors | [`compute_spectrum`](@ref) | [scorecard](reports/scorecard.md) |
| Figs. 4, 7, 9 captions | same-``J`` mixing angles | spectroscopic mixing | [mixing angles](reports/mixing_angles.md) |
| Table III, Fig. 5 | isoscalar annihilation mixing | [`compute_isoscalar_spectrum`](@ref) | [Table III](reports/table_iii_mixing_audit.md) |
| Tables IV–V | strong decays ``A \to B\pi`` | the paper's single-Gaussian reduction | [Table V](reports/table_v_reproduction.md) |
| Table VI | radiative transitions | [`PhotonEmission`](@ref) on solved waves | [Table VI](reports/table_vi_photon_decays.md) |
| Table VII | ``gg``, ``\gamma\gamma``, ``e^+e^-``, weak decay constants, charge radii | annihilation operators on solved waves | [Table VII](reports/table_vii_annihilation_em.md) |

[Results at a glance](@ref) summarizes each comparison. [Equations to code](@ref)
maps each equation of the paper to the code that implements it.
