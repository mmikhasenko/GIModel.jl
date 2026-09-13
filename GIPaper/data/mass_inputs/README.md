# Experimental mass inputs for the modern GI decay tables

`pdg-2026.csv` is a pinned selection from the official production PDG 2026 SQLite
snapshot. `source.json` records its URL, SHA-256, retrieval date, citation and
CC BY 4.0 license. `selection.csv` explicitly chooses each PDG mass property;
this matters for rho and K* where several incompatible mass definitions exist.
The PDG display field preserves published uncertainties (including multiple errors).
Calculations use the central database values, converted from MeV to GeV.

`assignments.csv` maps table-context/label pairs to these particles. A blank key
means unresolved experimental correspondence, not zero mass. `averages.csv`
defines arithmetic charge averages, repeating pi+ twice to give the triplet
average. PDG spectral assignments and the remix's GI-state correspondence are
separate statements. Broad scalar/effective channels and ambiguous excited states
are left unavailable rather than silently assigned a nearby resonance.

`historical.csv` preserves old audit masses and their provenance. They are used
only for historical tooltips and reconstructed historical q. A legacy script
constant is not evidence of an original GI experimental input. In particular,
several old bottomonium inputs were GI model estimates. The eta_r mapping is now
consistently eta(1295); eta'_r is deferred because historical eta(1440) does not
specify a unique modern state. Charged weak pi and neutral two-photon pi are distinct.

GIPaper's `experimental_mass(context, label)` returns `nothing` for registered
unavailable states and throws for unknown labels. It never selects a model or
historical mass. The spectrum still uses model eigenvalues; mock-meson masses
remain wavefunction expectations. PDG mass uncertainties are recorded but not yet
propagated into observables. Legacy measured-width comparisons are not part of
this mass update and must not be described as PDG 2026 widths.

To refresh the reviewed selection reproducibly:

```
python3 GIPaper/scripts/import_pdg_masses.py /path/to/pdg-2026.0.sqlite
```

The importer does not modify assignments, averages or historical data. A new PDG
edition requires reviewing the selected mass definitions and updating the edition.
