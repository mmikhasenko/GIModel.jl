# Code architecture (GIModel)

This note is for anyone opening `src/` or `scripts/` after refactors that split
**solver switches** from **constituent masses** and threaded **`ConstituentMasses`**
through the radial Hamiltonian and spin-dependent helpers. It replaces an older
mental model where masses lived inside `GIParameters` or were passed only as raw
`(m1, m2)` floats at every call site.

## Module layout (`src/GIModel.jl`)

Includes are grouped intentionally:

1. **Constants and core types** — `constants.jl`, `model_objects.jl` (`ConstituentMasses`,
   `RadialWaveOnUniformMesh`, fine-structure helpers, etc.).
2. **Setup / bookkeeping** — TOML → `GIParameters` (`parameters.jl`), `[masses]` →
   `QuarkMassTable` (`quark_mass_table.jl`), reference CSV rows (`reference_spectrum.jl`),
   parsing `quark_content` and attaching masses (`masses_from_content.jl`).
3. **Numerics** — potentials, Hamiltonian, `channel_solution`, contact / fine structure,
   `sector_workflow.jl`.

## Loading inputs

- **`load_parameters(path)`** → **`GIParameters`** only (potential, smearing switches,
  relativistic factors, fine-structure flags). It does **not** carry quark masses.
- **`load_quark_masses(path)`** → **`QuarkMassTable`** (`Dict{String,Float64}` with keys
  `"u"`, `"d"`, `"q"`, `"s"`, `"c"`, `"b"` in GeV).
- **`load_parameters_and_quark_masses(path)`** → `(GIParameters, QuarkMassTable)`. This is
  what **scripts** and most tests use: one TOML file with both `[potential]` / `[masses]` /
  related tables.

Quark masses are defined under **`[masses]`** in `data/parameters.provisional.toml`
(Table II–style MeV fields); they are **not** fields on `GIParameters`.

## Reference rows vs rows with masses

1. **`load_reference_spectrum(csv)`** → `Vector{ReferenceState}` (labels, `n`, `L`, `J`,
   reference mass, `quark_content`, etc.—no masses yet).
2. **`attach_constituent_masses(mq, reference_states, fallback_mass_GeV)`** →
   **`Vector{ReferenceStateWithMasses}`**, using `parse_quark_masses` /
   `resolve_constituent_masses` from `quark_content` and the flavor table `mq`.
   Equal-mass quarkonia scripts typically pass `mq["c"]` or `mq["b"]` as the fallback
   when the CSV row already pins both flavors.

## Sector workflow (`src/sector_workflow.jl`)

- **`compute_sector(params, annotated::Vector{ReferenceStateWithMasses}; …)`** builds one
  finite-difference radial solve per distinct **`RadialChannelKey`**:
  **`ConstituentMasses` + orbital letter `L`** (rounded masses define cache identity).
- **`compare(computed::SectorComputation, annotated::Vector{ReferenceStateWithMasses}; …)`**
  maps each reference row to the cached channel, picks radial level `n`, and adds contact /
  fine-structure shifts.

So batch drivers never pass “sector name” or flavor enums into `compute_sector`; all mass
information is already on each **`ReferenceStateWithMasses`**.

## Radial and central-potential API

- **`channel_solution(params, masses::ConstituentMasses, L::Integer; …)`** — spin-independent
  FD solve for one channel.
- **`relativistic_hamiltonian` / `nonrelativistic_hamiltonian`** take **`ConstituentMasses`**.
- **`central_potential_values`** and **`potential_diagonal`** accept either **`ConstituentMasses`**
  or **`(m1, m2)`** reals; prefer **`ConstituentMasses`** in new code for consistency with the
  Hamiltonian and **`RadialChannelKey`**.

## Reduced mass

**`reduced_mass(m::ConstituentMasses)`** (and related helpers) live next to **`ConstituentMasses`**
in `model_objects.jl`, not as stray utilities.

## What *not* to assume anymore

- Masses are **not** stored on **`GIParameters`**.
- **`compute_sector`** is not driven by a separate “flavor” argument; use annotated rows.
- Prefer **`ReferenceStateWithMasses`** at the boundary between CSV bookkeeping and numerics;
  avoid re-parsing `quark_content` inside the solver loop.

Keeping **`docs/formula_map.md`** aligned with this file is part of maintaining a reproducible
audit trail.
