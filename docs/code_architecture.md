# Code architecture (GIModel)

This note is for anyone opening `src/` or `scripts/` after refactors that split
**solver switches** from **constituent masses** and threaded **`ConstituentMasses`**
through the radial Hamiltonian and spin-dependent helpers. It replaces an older
mental model where masses lived inside `GIParameters` or were passed only as raw
`(m1, m2)` floats at every call site.

## Module layout (`src/GIModel.jl`)

Includes are grouped intentionally:

1. **Constants and core types** — `constants.jl`, `model_objects.jl`.
2. **Parameter bookkeeping** — TOML → `GIParameters` (`parameters.jl`), `[masses]` →
   `QuarkMassTable` (`quark_mass_table.jl`).
3. **Numerics** — potentials through Appendix-A status (`running_coupling.jl` … `appendix_a_status.jl`),
   Hamiltonian, `channel_solution`, contact / fine structure.
4. **State mixing layer** — **`state_mixing.jl`** (`BasisState`, `MixingBlock`,
   `MixingResult`, `diagonalize_mixing_block`) owns generic mass-matrix
   bookkeeping for same-`J`, tensor, flavor, or radial mixing blocks.
5. **Sector cache types + diagnostic sweep** — **`sector_solver.jl`** (`RadialChannelKey`,
   **`ChannelRadialSolution`**, **`SectorComputation`**, **`solve_sector`**).
6. **Sector batch solves + comparison + reports** — **`sector_comparison.jl`**
   (**`compute_sector`**, **`compare`**, **`write_residual_report`**).
7. **IO (last includes)** — **`reference_state.jl`** (**`ReferenceState`**, **`ReferenceStateWithMasses`**,
   **`load_reference_spectrum`**),
   **`masses_from_content.jl`** (**`parse_quark_masses`**, **`resolve_constituent_masses`**, **`attach_constituent_masses`**).

**`compute_sector`** / **`compare`** take **`AbstractVector`** rows with `.constituent_masses` and `.state`
(typically **`ReferenceStateWithMasses`** from **`attach_constituent_masses`**), so **`sector_comparison.jl`**
can load before reference structs are defined; CSV reading and string→mass helpers stay grouped here at the end.

## Loading inputs

- **`load_parameters(path)`** → **`GIParameters{FiniteDifferenceBasis}`** by default
  (potential, smearing switches,
  relativistic factors, fine-structure flags). It does **not** carry quark masses.
  The basis is a type parameter; `HarmonicOscillatorBasis` is already
  implemented and tested as an alternate paper-style oscillator expansion path,
  while `FiniteDifferenceBasis` remains the default production basis.
- **`load_quark_masses(path)`** → **`QuarkMassTable`** (`Dict{String,Float64}` with keys
  `"u"`, `"d"`, `"q"`, `"s"`, `"c"`, `"b"` in GeV).
- **`load_parameters_and_quark_masses(path)`** → `(GIParameters, QuarkMassTable)`. This is
  what **scripts** and most tests use: one TOML file with both `[potential]` / `[masses]` /
  related tables.

Quark masses are defined under **`[masses]`** in `data/parameters.provisional.toml`
(Table II–style MeV fields); they are **not** fields on `GIParameters`.

## Reference rows vs rows with masses

1. **`load_reference_spectrum(csv)`** (**`reference_state.jl`**) → `Vector{ReferenceState}` (labels, `n`, `L`, `J`,
   reference mass, `quark_content`, etc.—no masses yet).
2. **`attach_constituent_masses`** (**`masses_from_content.jl`**) → **`Vector{ReferenceStateWithMasses}`**,
   calling **`resolve_constituent_masses`** (**`masses_from_content.jl`**) with `String(state.sector)` and
   `String(state.quark_content)` per row. Equal-mass quarkonia scripts typically pass `mq["c"]` or `mq["b"]` as the fallback
   when the CSV row already pins both flavors.

## Sector workflow (`src/sector_solver.jl`, `src/sector_comparison.jl`)

- **`compute_sector`** (**`sector_comparison.jl`**) builds one finite-difference radial solve per distinct
  **`RadialChannelKey`**: **`ConstituentMasses` + orbital letter `L`** (rounded masses define cache identity).
- **`compare`** (**`sector_comparison.jl`**) maps each reference row to the cached channel, picks radial
  level `n`, and adds contact / fine-structure shifts.
- **`write_residual_report`** (**`sector_comparison.jl`**) turns **`compare`** output into markdown; script
  callers pass booleans that mirror the active **`GIParameters`** path for the prose header.

So batch drivers never pass “sector name” or flavor enums into `compute_sector`; all mass
information is already on each **`ReferenceStateWithMasses`**.

## Mixing layer (`src/state_mixing.jl`)

Mixing is intentionally above the pure radial/basis-state calculations. A
**`MixingBlock`** stores a list of **`BasisState`** labels and an arbitrary-size
Hermitian mass matrix in GeV; **`diagonalize_mixing_block`** returns sorted
eigenmasses and eigenvectors with stable phases for reports. This keeps
mechanism-specific matrix construction separate from the linear algebra:

- unequal-mass antisymmetric spin-orbit can build `^1L_L`/`^3L_L` blocks;
- tensor interactions can later build same-`J`, different-`L` blocks such as
  `^3S_1`/`^3D_1`;
- isoscalar or radial mixings can use the same block/eigenstate reporting path.

The current `same_j_mixing` helper in `spin_fine_structure.jl` is a two-state
convenience wrapper around this generic layer. The gap is not the absence of a
mixing layer; it is wiring mechanism-specific block builders into `compare`
and physical sector assignment.

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
