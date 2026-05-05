# Open-Flavor Residual Investigation

Generated from the local investigation on 2026-05-05. This note records the
paper clues and scratch adjustment scans used to diagnose why the open-flavor
sectors are worse than heavy quarkonium, even though the central spectrum is
already close.

## Scope

Open-flavor sectors considered:

- `strange`: `data/reference_spectrum_strange.csv`
- `charmed`: `data/reference_spectrum_charmed.csv`
- `b_flavored`: `data/reference_spectrum_b_flavored.csv`

Heavy-heavy control sectors:

- `charmonium`: `data/reference_spectrum_charmonium.csv`
- `bottomonium`: `data/reference_spectrum_bottomonium.csv`

The baseline scores are those from the active GI-style path:

- closed-form smeared central `G~(r)` and `S~(r)`
- central Coulomb momentum sandwich
- contact hyperfine momentum sandwich
- fine-structure momentum sandwiches
- smeared fine-structure radial kernels
- Table II masses and parameters from `data/parameters.provisional.toml`

No source-code changes were made during the scans. The numbers below were
computed by scratch Julia scripts under `/private/tmp`.

## Paper Clues

The relevant paper guidance is concentrated in the model setup and Appendix A.

1. The paper warns that the nonrelativistic spin-dependent terms are too
   singular and must be regulated by the full relativistic potential. In the
   prose extraction:

   > The smearing of the potentials has the consequence of taming all of their
   > singularities...

   Nearby text also says the coefficients of the potentials become
   momentum-dependent.

2. The same section gives the most direct clue for the open-flavor problem:
   light-quark hyperfine should not behave like the naive `1/(m1 m2)` form.
   The prose extraction states that the hyperfine interaction of a light quark
   should not blow up like the nonrelativistic mass factor; it should have a
   finite limit controlled by confinement-scale momenta.

3. Appendix A classifies the relativistic corrections into:

   - momentum-dependent strengths,
   - nonlocality,
   - new `r` dependences through `Q` dependence.

   It then introduces Gaussian smearing via Eqs. (A7)--(A9) and the
   momentum-dependent factors in Eq. (A10).

4. Appendix A explicitly says each interaction could in principle have a
   distinct smearing function and more complicated energy-dependent factors.
   This matters because the current implementation uses the same Table II
   `sigma0, s` smearing family for contact and for the smeared kernels.

5. Table II gives the contact ambiguity parameter as `epsilon_c = -0.168`.
   The local extraction shows the sign as negative:

   - `data/table_ii_parameters.csv`
   - `data/raw/digitized_tables/table_ii_parameters/raw_page_text.txt`

   Therefore, a sign flip should be treated as a diagnostic clue, not a
   faithful parameter update.

## Baseline Pattern

The central spin-averaged open-flavor levels are already good. The residuals
are dominated by spin-dependent splittings.

### Charmed

From `docs/residual_reports/charmed_residuals.md`:

| multiplet | spin-averaged residual |
|---|---:|
| `1S` | `+0.0 MeV` |
| `2S` | `+4.0 MeV` |
| `1P` | `+2.3 MeV` |
| `1D` | `-12.3 MeV` |
| `1F` | `-43.2 MeV` |

But individual S-wave residuals show too much hyperfine splitting:

| state | residual |
|---|---:|
| `1^1S_0` (`c qbar`) | `-79.5 MeV` |
| `1^3S_1` (`c qbar`) | `+28.2 MeV` |
| `1^1S_0` (`c sbar`) | `-47.0 MeV` |
| `1^3S_1` (`c sbar`) | `+14.0 MeV` |

### B-Flavored

From `docs/residual_reports/b_flavored_residuals.md`:

| multiplet | spin-averaged residual |
|---|---:|
| `1S` | `-1.8 MeV` |
| `2S` | `+2.7 MeV` |
| `1P` | `-13.3 MeV` |
| `1D` | `-23.5 MeV` |
| `1F` | `-33.0 MeV` |

The b-flavored sector is already substantially better than strange/charmed.
The largest misses are in bottom-light higher-L plotted states:

| state | residual |
|---|---:|
| `1^3D_3` (`b qbar`) | `-56.4 MeV` |
| `1^3F_4` (`b qbar`) | `-78.0 MeV` |

These are more suggestive of spin-orbit / tensor sensitivity than of a central
potential failure.

### Strange

From `docs/residual_reports/strange_residuals.md`:

| multiplet | spin-averaged residual |
|---|---:|
| `1P` | `+4.6 MeV` |
| `2P` | `-1.7 MeV` |
| `1D` | `+5.1 MeV` |
| `2D` | `+0.3 MeV` |
| `1F` | `+2.3 MeV` |
| `1G` | `+1.6 MeV` |

But the S-wave hyperfine residuals are very large:

| state | residual |
|---|---:|
| `1^1S_0` | `-229.8 MeV` |
| `2^1S_0` | `-258.9 MeV` |
| `3^1S_0` | `-230.3 MeV` |
| `1^3S_1` | `+105.4 MeV` |
| `2^3S_1` | `+83.2 MeV` |
| `3^3S_1` | `+67.1 MeV` |

This is the clearest evidence that the contact hyperfine term is too strong for
light constituents.

## Mixing Check

For Fig. 9 b-flavored mesons, the plotted states do not require same-`J`
`^1P_1`/`^3P_1` mixing in the residual score.

The rendered Fig. 9 page states:

- significant spectroscopic mixing exists in the sector,
- `theta_1P^bq ~= -43 deg`,
- `theta_1P^bs ~= -45 deg`,
- `theta_1P^bc ~= -53 deg`,
- the `1+` states are not shown in the figure and are within `40 MeV` of their
  respective `2+` partners.

Therefore the current plotted-label comparison for Fig. 9 is legitimate without
mixing corrections. A separate mixing validation should reproduce the quoted
angles and the `1+` proximity statement, but those masses are not directly
available as plotted labels in Fig. 9.

## Adjustment Scans

### 1. Linear Spin-Term Rescaling

Keeping central eigenvalues fixed, fit the already computed spin terms as

```text
prediction = central
           + a_contact * contact_shift
           + a_so      * spin_orbit_shift
           + a_tensor  * tensor_shift
```

Baseline aggregate scores:

| set | rows | mean abs | max abs | rms |
|---|---:|---:|---:|---:|
| `open all` | 73 | `39.85` | `258.93` | `65.41` |
| `open heavy only` | 43 | `25.21` | `79.53` | `33.14` |
| `heavy-heavy` | 58 | `5.33` | `18.01` | `7.09` |

Least-squares multipliers:

| fit set | contact | spin-orbit | tensor |
|---|---:|---:|---:|
| `open all` | `0.478` | `0.356` | `0.428` |
| `open heavy only` | `0.596` | `0.204` | `0.874` |

The contact-only fit is:

| fit set | contact |
|---|---:|
| `open all` | `0.478` |
| `open heavy only` | `0.596` |

This points to over-strong contact hyperfine as the most robust first problem.

### 2. Simple Contact Multiplier

Trial: multiply only the contact hyperfine contribution by `0.60`.

| set | baseline mean abs | contact `0.60` mean abs |
|---|---:|---:|
| `open all` | `39.85` | `26.29` |
| `open heavy only` | `25.21` | `17.30` |
| `heavy-heavy` | `5.33` | `5.07` |

Per-sector:

| sector | baseline mean | contact `0.60` mean | baseline max | contact `0.60` max |
|---|---:|---:|---:|---:|
| `strange` | `60.83` | `39.16` | `258.93` | `169.93` |
| `charmed` | `33.02` | `21.13` | `79.53` | `69.52` |
| `b_flavored` | `17.02` | `13.29` | `78.03` | `78.03` |
| `charmonium` | `6.75` | `6.49` | `18.01` | `30.64` |
| `bottomonium` | `4.01` | `3.74` | `12.56` | `14.08` |

This is a good diagnostic improvement, but it should not be promoted as a
physics parameter by itself. It says the active contact implementation is too
large by roughly a factor of `1.7--2.1` in open-flavor channels.

### 3. Contact Momentum-Exponent Scan

The active contact momentum sandwich uses side exponent

```text
0.25 + 0.5 epsilon_c
```

with Table II `epsilon_c = -0.168`, so the current side exponent is `0.166`.

Keeping all other settings unchanged:

| `epsilon_c` | side exponent | open mean | strange | charmed | b_flavored | heavy-heavy mean |
|---:|---:|---:|---:|---:|---:|---:|
| `-0.300` | `0.100` | `53.98` | `86.17` | `41.87` | `20.69` | `6.14` |
| `-0.168` | `0.166` | `39.85` | `60.83` | `33.02` | `17.02` | `5.33` |
| `-0.100` | `0.200` | `34.32` | `51.28` | `29.24` | `15.41` | `5.00` |
| `+0.000` | `0.250` | `27.80` | `40.36` | `24.47` | `13.34` | `4.58` |
| `+0.100` | `0.300` | `23.44` | `33.73` | `20.49` | `11.82` | `4.39` |
| `+0.200` | `0.350` | `22.20` | `32.04` | `18.88` | `11.62` | `4.37` |
| `+0.300` | `0.400` | `24.28` | `34.91` | `20.79` | `12.75` | `4.63` |
| `+0.400` | `0.450` | `26.54` | `38.30` | `22.72` | `13.74` | `4.91` |
| `+0.600` | `0.550` | `29.96` | `42.79` | `26.32` | `15.44` | `5.38` |

This is the most informative scan. The score wants a much stronger contact
momentum suppression, roughly equivalent to `epsilon_c ~= +0.1--0.2`.

However, Table II gives `epsilon_c = -0.168`. The result should therefore be
interpreted as evidence that the active contact implementation is not yet
GI-equivalent, not as permission to change the published parameter.

Likely audit targets:

- sign or placement of the contact `m/E` ambiguity factor,
- whether the side exponent should be applied in the current FD `p^2` basis in
  the same way as in the paper's HO basis,
- missing nonlocal/operator-ordering terms that further suppress light-quark
  contact matrix elements,
- whether the current contact kernel should use a distinct effective smearing
  from the central/fine-structure smeared kernels.

### 4. Smearing Width Scan

The paper says the mass dependence of the `sigma0` term in Appendix A is
significant mainly for pseudoscalar mesons. This motivated a direct scan of
`sigma0` and `s`.

| `sigma0` | `s` | open mean | strange | charmed | b_flavored | heavy-heavy |
|---:|---:|---:|---:|---:|---:|---:|
| `1.20` | `1.55` | `32.39` | `44.16` | `29.19` | `18.92` | `5.16` |
| `1.40` | `1.55` | `33.29` | `50.22` | `27.39` | `15.27` | `5.18` |
| `1.60` | `1.55` | `36.34` | `55.80` | `29.82` | `15.35` | `5.25` |
| `1.80` | `1.55` | `39.85` | `60.83` | `33.02` | `17.02` | `5.33` |
| `2.00` | `1.55` | `43.11` | `65.43` | `36.01` | `18.65` | `5.47` |
| `1.80` | `1.20` | `38.78` | `60.17` | `31.02` | `16.35` | `5.46` |
| `1.80` | `1.80` | `41.18` | `61.39` | `34.56` | `19.24` | `5.88` |
| `1.60` | `1.20` | `35.68` | `54.99` | `27.50` | `16.65` | `5.72` |

Lowering `sigma0` helps, but less cleanly than increasing contact momentum
suppression. It also changes the central smeared potential because the active
central path uses the same smearing family.

## Interpretation

The central potential is not the primary open-flavor failure.

Evidence:

- spin-averaged centers are already close in strange and charmed sectors;
- b-flavored centers are also acceptable;
- heavy-heavy remains excellent.

The dominant mismatch is the spin-dependent sector:

1. Contact hyperfine is too strong for light constituents.
2. Higher-L bottom-light misses probably involve spin-orbit / tensor details,
   especially the delicate vector--Thomas cancellation.
3. Same-`J` mixing is implemented as a diagnostic helper, but it is not the
   explanation for the plotted Fig. 9 residuals because those mixed `1+` states
   are explicitly not plotted.

The most promising next step is a contact hyperfine audit, not a central
potential retune.

## Recommended Next Steps

1. Audit the contact momentum sandwich against Appendix A Eq. (A10), including
   exponent placement, side-factor convention, and whether the current FD
   `p^2` eigenbasis implements the same operator ordering intended in the
   paper's HO calculation.

2. Add an explicit diagnostic option for a contact-strength multiplier, but keep
   it marked non-physical. Use it only to demonstrate that a `~0.5--0.6`
   reduction repairs the open-flavor pattern.

3. Add a targeted score report that splits residuals into:

   - spin-averaged central centers,
   - contact-only S-wave splittings,
   - spin-orbit/tensor splittings,
   - mixed-state diagnostics where the paper actually quotes angles.

4. Use the Fig. 9 caption as a separate validation target:

   - reproduce `theta_1P^bq ~= -43 deg`,
   - reproduce `theta_1P^bs ~= -45 deg`,
   - reproduce `theta_1P^bc ~= -53 deg`,
   - check that the resulting `1+` eigenmasses lie within `40 MeV` of their
     `2+` partners.

5. Only after contact is audited, revisit global `k_spin_orbit` / `k_tensor`.
   Least-squares spin-term fits can improve open flavor, but they degrade the
   heavy-heavy control more than the contact-focused adjustment does.
