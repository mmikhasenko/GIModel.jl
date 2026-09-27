# Mixing: what must be matched?

> **Closed / historical investigation brief (2026-09-27).** The pre-fix
> numbers and proposed checks below describe the original investigation, not
> current failures. The all-L contact defect is fixed; the remaining large
> GI85 caption mismatch is a documented historical reference discrepancy.
> Current targets retain GI85 and add GK91 side by side:
> [target ledger](../../data/mixing_angle_targets.md),
> [comparison](../residual_reports/mixing_angles.md),
> [resolution](mixing_composition_investigation.md).
> The separate ψ(3.82) tensor-composition question remains open.

## Decision

Use published eigenstate compositions as the primary mixing targets, alongside
final masses. Treat reconstructed two-level repulsion as a conditional diagnostic,
not a published GI before/after energy shift. A failure to match composition is
already evidence of a reproduction discrepancy without using any digitized mass.
It is not, by itself, proof of a programming error.

## Source and mechanism ledger

The local source is Godfrey–Isgur (1985), Sec. II, p. 192 and Figs. 3–9.
The original page images for Figs. 3, 4, and 6 were visually checked for this note,
not just their OCR transcriptions.

| Mechanism | Direct publication target | Comparison quantity |
|---|---|---|
| Antisymmetric spin–orbit, singlet/triplet at L=J | Fig. 4: six angles; Fig. 7: four; Fig. 9: three | Lower-state singlet probability, projected angle with fixed conventions, and norm outside the quoted pair |
| Tensor, triplet L=J−1 ↔ J+1 | Fig. 3: vector at 1.45 GeV ≈ 1.00(2S)+0.04(1D); Fig. 4: vector at 1.58 GeV has the same quoted amplitudes | Individual component magnitudes/ratios; dominant component identifies the state |
| Tensor, charmonium | Fig. 6: vector at 3.82 GeV ≈ 1.00(1D)+0.01(1S)−0.03(2S)−0.01(3S) | Three separate S-wave amplitudes relative to the dominant D component; not a single angle |
| Isoscalar annihilation | Table III, explicitly in eigenstates in the absence of annihilation | Flavor/radial amplitude vectors, separately for P1 and P2; see the existing Table III audit |

The 13-angle ledger remains authoritative at
[mixing_angles.md](../residual_reports/mixing_angles.md). The repulsion scatter
uses only the first mechanism: eight of those angle rows also have both final
mass labels. It does not use the tensor coefficients or Table III.

The tensor numbers are amplitudes: 0.04 corresponds to approximately 0.16%
probability, not 4%. Printed 1.00 does not imply exactly zero other components.
Rounded, incomplete vectors must not be treated as exactly normalized. Relative
signs require a consistent basis-phase convention; component magnitudes are the
first robust target. The existing census does not establish whether these three
specific tensor vectors agree: that needs a component-by-component comparison.
Table III disagreement is a separate annihilation question, not corroborating
evidence for an antisymmetric spin–orbit bug.

## Does our procedure match?

Sec. II first diagonalizes fixed (J,L,S) sectors, then diagonalizes tensor and
antisymmetric spin–orbit couplings in those eigenstates. “Perturbatively” here
does not mean replacing the second diagonalization by a second-order energy
formula. GI says both bases are expanded to convergence.

The current production path has this structure:

- `fixed_spectrum` calls `fixed_channel_solution` for each (L,S,J), retaining
  its eigenvalues and native radial eigenvectors, not central-only waves.
- `_apply_same_j_spin_orbit_mixing!` constructs diagonal fixed-sector masses
  and all requested singlet/triplet cross-radial couplings at fixed L=J.
- `_apply_tensor_mixing!` constructs the corresponding triplet L=J±1 block.
- Both call `diagonalize_mixing_block` and retain the complete eigenvectors.

These mechanisms act in disjoint spectroscopic subspaces here. Their sequential
application is not an ordering approximation between competing couplings.
This verifies architectural agreement, not independent correctness of kernels.

Important limitation: second-stage blocks contain the requested `levels`, not
automatically every radial eigenstate available to the solver. Increasing
`nlevels_per_channel` alone does not enlarge these blocks. The ten-panel census
is a plotting envelope, not a demonstrated converged second-stage basis.

## Clear evidence without a mass inversion

For the paper's lower-state convention, the projected singlet probability is
cos²(theta). Using the existing full-envelope computation:

| Lower 1P state | GI projected singlet probability | Computed projected singlet probability |
|---|---:|---:|
| Kaon | 68.7% | 0.20% |
| D | 57.0% | 4.60% |
| Ds | 51.7% | 3.80% |

These are direct composition discrepancies, independent of mass digitization.
The GI two-state expressions are approximate; our probabilities are conditional
on the same-n pair. Report outside-pair norm separately, not silently normalize
it away. Basis sign changes cannot repair these probability discrepancies.
Compare lower with lower: swapping eigenstates is not a permissible fit choice.

## What is and is not determined by composition?

For an isolated real symmetric two-state matrix with gap Δ=a−b and coupling V,
an angle fixes a mixing ratio, not an energy scale. Multiplying both Δ and V
by the same factor leaves the eigenvectors unchanged. With physical splitting
S supplied as well, the reconstruction is unique up to basis conventions:

    |V| = S |sin(2θ)| / 2
    |Δ| = S |cos(2θ)|
    δ = S [1 − |cos(2θ)|] / 2

Here δ is the outward shift relative to the ordered diagonal energies. A useful
mass-independent strength measure is δ/S = min(cos²θ, sin²θ), conditional on
two-state isolation. It must not be labeled a measured shift in MeV.

In a larger block, one composition vector—even known exactly—and its energy
do not determine the Hamiltonian. A complete orthonormal eigenvector matrix U
and all eigenvalues would give H=U diag(E) U†, but GI's caption summaries are
not such a dataset. Missing components and eigenvalues preclude a unique
before/after reconstruction. Small omitted probability alone is not an energy
error bound without information about the omitted energy scales.

## Investigation acceptance targets

1. Match all 13 published spin–orbit compositions, initially using probabilities
   and then signed angles after one global convention check. Keep all rows,
   including the five with no printed mass pair.
2. Compare the three explicitly quoted tensor vectors component by component.
   Do not infer failure from aggregate tensor-shift statistics.
3. Separately match published final masses and splittings; do not fit inferred
   shifts as if they were extra independent observations.
4. Establish convergence of both fixed-sector waves and second-stage radial
   blocks for these states, tracking identity by composition/overlap rather
   than relying only on precursor labels.
5. Only then use reconstructed |V|, |Δ|, and δ to localize disagreement. Printed
   rounding is a sensitivity scale, not a statistical error bar or a bound on
   the omitted-state approximation.

No new solver run or kernel correction was made for this target-definition
review. The probability examples are transformations of the existing census.
