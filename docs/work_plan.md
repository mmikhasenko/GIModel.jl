# Table VI completion record

The radiative-decay implementation now covers all **79 canonical rows**:
42 M1, 35 E1, and 2 M2. Numerical results are in
`GIPaper/reports/table_vi_photon_decays.md` and its CSV companion.

Completed work:

1. Canonical identity mapping, unique-row validation, and explicit state/mass inputs.
2. Excited bottomonium E1 and hindered M1 transitions, including 3S -> 2S.
3. Excited eta/eta-prime and mixing-induced channels, through the final
   four-flavor MixedSpectrum rather than report-local vectors.
4. Shared native fixed-channel HO treatment, including spin-distorted P waves;
   an independent FD run uses the same observable pipeline.
5. Pure-flavor charge normalization, shared precursor-based phase conventions,
   recoil/form-factor prescriptions, and M1 width conversion.
6. Encoding corrections: the neighboring bottomonium target, the A2 mass
   outside the square root, and the spin-flip E1 classification of A1 -> pi.

The paper's +0.01 μN pi0-eta contribution is retained as an explicitly labelled
external input. No parameter is fitted to the new rows.

Completion means all rows are computed and their conventions and residuals
reported. It does not mean exact numerical reproduction of the excited eta
channels or every cancellation-sensitive transition. Further study of those
residuals is physics follow-up, not missing photon-decay coverage.

Run `julia GIPaper/checks/audit_table_vi_photon_decays.jl` for the HO audit,
or set `GI_TABLE_VI_SOLVER=fd` for its independent comparator.

## Verification

`scripts/verify_project.sh` passed on 2026-09-09, including 1,145 GIModel tests,
2,505 GIPaper tests, the native-HO and independent-FD convergence gates,
all reproduction audits, and the 103-unit manifest/link check.

The complete Table VI FD cross-check has median differences from HO of
0.099% (M1), 0.049% (E1), and 0.092% (M2). The detailed comparison is in
`GIPaper/reports/table_vi_solver_comparison.md`.
The shared phase correction also resolves the former Table VII excited-eta
two-photon sign discrepancies (now 4/4 signs), without changing magnitudes.
