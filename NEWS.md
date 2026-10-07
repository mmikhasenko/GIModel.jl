# News

## v0.4.2 (unreleased)

### Fixed

- **Singlet E1 rates were 18× too large.** `e1_angular_coefficient` gave
  ¹P₁ → ¹S₀ an angular factor of √2 instead of 1/3. E1 does not act on spin, so
  the singlet factor equals the triplet ³P_J → ³S₁ one. This affects every
  singlet E1 width from `PhotonEmission`: h_c → η_c γ drops from 5.7 MeV to
  about 0.32 MeV, and the same applies to h_b, b₁ and h₁. The √2 came from the
  1985 Table VI row B → πγ = (√2 q/3) E1^u, which is inconsistent with the
  paper's own A₂ → ργ = (q/9) E1^u. The GIPaper Table VI audit now gives
  +0.149 for B → πγ against the printed +0.63 and explains the difference.
  ([#16](https://github.com/mmikhasenko/GIModel.jl/pull/16))
- **E1 widths ignored the quark masses for unequal-mass mesons.** Each quark's
  charge entered the E1 dipole with weight ½. The dipole is measured from the
  centre of mass, so the weight is m_j/(m_i+m_j). This affects every heavy-light,
  B_c and kaon E1 width. For example, B⁺(1³P₂) → B*⁺γ goes from 0.4 to 133 keV
  and D_s(1³P₂) → D_s*γ from 45 to 4 keV. Quarkonium and the GI 1985 Table VI
  audit are unchanged. With the m/E exponent off, the widths now match
  Godfrey's later GI-model papers on D, D_s, B and B_c to within 1.2%
  (`GIPaper/checks/audit_heavy_light_e1.jl`).
  ([#24](https://github.com/mmikhasenko/GIModel.jl/issues/24))

### Added

- Inverse singlet E1 transitions (¹S₀ → ¹P₁) are supported. Before, they threw
  an error.

### Changed

- Photon angular factors are derived from L⊗S coupling instead of being
  transcribed from the paper's Table VI. The E1 coefficients and the J
  dependence of the spin-flip denominators (√60 m, 6 m) now come from the same
  Clebsch–Gordan and spin algebra as strong decays. No computed value changes.
  The overall spin-flip normalization keeps the GI constant; it is an open
  question in [#18](https://github.com/mmikhasenko/GIModel.jl/issues/18).
  ([#19](https://github.com/mmikhasenko/GIModel.jl/pull/19))

## v0.4.1 (2026-09-29)

- Agent skill for using GIModel.jl (`skills/gimodel`).
  ([#14](https://github.com/mmikhasenko/GIModel.jl/pull/14))
- Removed scratch notebooks and a stale tutorial copy.
  ([#15](https://github.com/mmikhasenko/GIModel.jl/pull/15))

## v0.4.0 (2026-09-29)

First public release. See the
[release notes](https://github.com/mmikhasenko/GIModel.jl/releases/tag/v0.4.0).
