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
- **Spin-flip photon rates were 2× too large.** The ³P₂ → ¹S₀ (M2) and
  ³P₁ → ¹S₀ (E1) kernels used GI 1985's denominators √60 m and 6m. The
  magnetization current gives √120 m and √72 m. The same operator, evaluated by
  the new `MultipolePhotonEmission`, reproduces the standard M1 and E1 widths and
  the Karl–Meshkov–Rosner M2/E1 ratios, and the CLEO and BESIII M2 fractions in
  charmonium prefer it to the √2-larger normalization (χ² 20 against 38 for 8
  measurements). Each emitter also gets its centre-of-mass weight 2m_j/(m_i+m_j),
  as for E1. In the GI paper audit, A₂ → πγ and A₁ → πγ drop by 1/√2 and
  K*(1420) → Kγ by 0.82. GI fitted their m/E exponent to A₂ → πγ, so that fit
  absorbed the factor. ([#27](https://github.com/mmikhasenko/GIModel.jl/issues/27))
- The internal spherical Bessel function lost relative precision for orders
  l ≥ 9 at small arguments and overflowed at l ≥ 17. No published number changes.

### Added

- `MultipolePhotonEmission`: photon emission from the quark convection and spin
  currents, exact in the photon momentum and resolved into every multipole
  (E1, M1, E2, M2, E3, …). It handles any orbital angular momentum and mixed
  states, including the ³S₁–³D₁ admixtures that `PhotonEmission` rejects.
  `multipole_fractions` returns the normalized amplitudes measured in
  χ_cJ → J/ψγ and ψ(2S) → γχ_cJ. New GIPaper check
  `audit_charmonium_multipoles.jl` compares them with CLEO and BESIII.
  ([#27](https://github.com/mmikhasenko/GIModel.jl/issues/27))
- Inverse singlet E1 transitions (¹S₀ → ¹P₁) are supported. Before, they threw
  an error.

### Changed

- Photon angular factors are derived from L⊗S coupling instead of being
  transcribed from the paper's Table VI. The E1 coefficients and the J
  dependence of the spin-flip denominators (√60 m, 6 m) now come from the same
  Clebsch–Gordan and spin algebra as strong decays. No computed value changes.
  The overall spin-flip normalization kept the GI constant, corrected in #27
  above. ([#19](https://github.com/mmikhasenko/GIModel.jl/pull/19))

## v0.4.1 (2026-09-29)

- Agent skill for using GIModel.jl (`skills/gimodel`).
  ([#14](https://github.com/mmikhasenko/GIModel.jl/pull/14))
- Removed scratch notebooks and a stale tutorial copy.
  ([#15](https://github.com/mmikhasenko/GIModel.jl/pull/15))

## v0.4.0 (2026-09-29)

First public release. See the
[release notes](https://github.com/mmikhasenko/GIModel.jl/releases/tag/v0.4.0).
