# Discrepancies beyond orbital P-wave fine structure

## Statement

There is no evidence for one common numerical cause. Two controlled substitutions
localize major decay residuals to the isoscalar composition used by the comparison
pipeline. Other discrepancies remain unresolved or compare different definitions
and inputs. “Cancellation” identifies sensitivity; it does not by itself identify
an implementation error or establish what GI computed.

The largest independent spectrum problem is the **pseudoscalar P1 annihilation
model**, involving radial S waves, not orbital P waves. Its printed compositions
also fail a consistency check for eigenvectors of a single symmetric matrix.
Consequently we cannot attribute its mismatch solely to GI's numerical precision.

## Excited eta moments: composition tested independently of kinematics

The executable [composition substitution](excited_eta_moment_results.md) replays
the exact native HO waves recorded by `trace_rate_inputs.jl`. All six original
moments are recovered to 1e-9 absolute / 1e-8 relative tolerance. It then replaces
only the parent's coefficients with the canonical Table III P1 entries. Both
literal and normalized versions are reported; no transition strength is fitted.

With literal printed coefficients, eta_r→rho changes from +0.297925 to +0.565579
μN, against +0.57; eta_r→phi changes from +0.502568 to +0.371551, against +0.37.
The eta'_r→phi result changes from −0.026182 to −0.293162, against −0.29.
The tiny eta'_r→rho result improves from −0.218603 to +0.011777, against +0.008,
but retains a material relative residual. Its native sum cancels by a factor
7.94; the phi sum cancels by 21.77. Large relative errors here do not imply
large errors in each individual radial overlap.

**Conclusion:** the composition difference accounts for most of these magnitude
residuals using unchanged local waves. Changing experimental photon momentum
cannot fix them: none of these six M1 moments includes a recoil correction or
depends on q. A width inferred from the moment does depend on q³. Missing parent
masses therefore obstruct widths, not these particular moment comparisons.
This corrects the earlier kinematic explanation in the paper-remix overview.

One sign remains: eta_r→omega gives +0.179867 against the printed −0.18, while
eta_r→rho has the correct positive sign. A global parent rephasing cannot repair
both. The signs were checked against the Table VI page image (PDF page 24,
journal page 212). No source correction is made: a channel convention or a
printed sign error remains possible; the current test does not decide between
them. Individual transition signs are convention dependent, so this is not by
itself an observable width failure.

The existing [HO/FD control](../residual_reports/table_vi_solver_comparison.md)
also retains the large excited-eta residuals. That is independent evidence
against a basis-choice explanation, though it shares the physical operator and
mixing prescription. The archived control has slightly older heavy-state
kinematics and is not a fresh run of every current row.

## What can Table III actually constrain?

The source image on PDF page 11 (journal page 199) confirms the canonical numbers.
For the two excited P1 states, the seven printed coefficients have squared norms
0.914210 and 0.846110, and mutual inner product −0.111790. These deficits are much
larger than two-decimal rounding. Normalization alone does not restore
orthogonality. Tiny listed heavy components cannot supply the missing norm.

The caption calls the compositions approximate. It explicitly exempts the
mass-dependent **P2** construction from orthogonality; that does not explain the
same defect in P1. Thus these printed P1 vectors cannot all be exact eigenvectors
of one symmetric matrix in the stated orthonormal unmixed basis. Additional
unlisted components, a different approximation, or a source inconsistency must
be considered before treating every component as an exact solver target.

The [canonical Table III audit](../residual_reports/table_iii_mixing_audit.md)
still shows a 169 MeV discrepancy for the fourth P1 mass. The coefficient
substitution neither repairs nor explains that mass. The next discriminating
investigation is the original annihilation-block construction and the meaning
of its approximate compositions, not a parameter scan. No new Hamiltonian bug
has been demonstrated by this study.

## Tensor two-photon amplitude: ideal flavor is the dominant omission

The [two-flavor substitution](tensor_photon_results.md) uses the existing
Table VII central-wave amplitudes. The operator's exact M^(3/2) dependence
transports both flavor amplitudes to the same external mass, then combines
them with normalized printed Table III coefficients.

For f'_2 the nonstrange contribution is +0.172524 and the strange contribution
is −0.432574 keV^(1/2). Their sum is −0.260050, compared with the ideal-flavor
baseline −0.433639 and GI −0.250000. This removes 94.5% of the absolute residual.
The small nonstrange coefficient matters because its electromagnetic amplitude
is large. The f_2 control changes only from −1.898367 to −1.915014, against −1.90.

**Conclusion:** the large f'_2 discrepancy is predominantly explained by the
baseline's ideal-flavor approximation in this controlled test. It is not evidence
for a 73% error in the strange radial integral. This is not yet a full prediction:
the test imports paper coefficients, omits heavy admixtures, and retains central
P waves. A production change should compose native physical tensor states
through the shared transition API and verify the remaining wave-treatment
dependence. The canonical baseline has not been silently replaced.

## Remaining families and limits of the conclusions

| Family | Evidence and clear status |
|---|---|
| Strange D/F/G masses and mixing | The strange-sector ledger gives centroids within 2 MeV but singlet-labelled residuals of 10–15 MeV and discrepant angles. This extends the splitting/mixing issue beyond P waves. It does not establish the exact faulty term; a new high-L operator/convergence study has not been performed here. |
| Bottomonium hindered E1 | Both existing HO and FD results retain the small positive Upsilon''→chi_b0 amplitude against −0.002. A node-sensitive overlap is implicated; historical GI precision or a particular wave distortion has not been demonstrated as the cause. Mass convergence alone cannot certify such a near-zero integral. |
| Strong-decay misses | The [13-row ledger](../residual_reports/table_v_reproduction.md#investigated-misses) tests threshold, momentum and mixing sensitivity. An inverted mass that matches a target is a sensitivity diagnostic, not proof of the historical mass used. This baseline uses fitted analytic amplitudes and externally prescribed angles; these misses do not independently test the GIModel eigensolver. |
| Neutral-kaon radius | −0.1029 versus −0.0900 fm² remains unresolved. Opposite-charge contributions explain sensitivity but have not quantitatively explained the residual. The smearing exponent was calibrated on the pion. |
| Table VII D-vector constants | psi'' is 0.0070 versus 0.0110; rhoD is 0.0154 versus 0.0190. The baseline uses central D waves. S–D composition and operator sensitivity require a dedicated controlled comparison; a good overall median does not resolve these rows. |
| Leptonic widths versus experiment | These are a separate comparison from reproducing GI's printed constants. Squaring constants amplifies residuals. Omitted radiative corrections are not demonstrated here to account quantitatively for the discrepancy and must not be used as a blanket explanation. |
| Table I confinement share | The local Hellmann–Feynman definition differs from an unspecified historical estimate. Failure to match is definition-dependent; inverse string tensions do not uniquely reconstruct the author's method. |
| Unavailable rows | Missing assignments or an unimplemented historical sector are coverage limitations, not numerical discrepancies. |

No new spectrum solve or parameter scan was needed for the two substitutions.
The generated files record input hashes. Refresh the native input traces first
if the solver, parameters, or comparison pipeline changes; then run the two
scripts linked from this directory's README. Full row-by-row comparison remains
owned by `residual_reports/`.
