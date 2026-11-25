# A penetrating probe distinguishes the source currents

The compact vacancy prepared in either stationary current mode gives a source-reversal-sensitive **exterior load reading when the probe can enter through two gateways**. The interaction disturbs that source. At supplied exposure 32, an independent excess is detected at port 48 with probabilities **0.406115753 and 0.428729753** for the opposite initial currents. No external flux is installed. A source/probe number aperture retains owners 0, 4 and 6, with the original native diagonal terms and all contacts inside that region. The coherent probe preparation is supplied. [Design and disclosures](PROTOCOL.md) · [Calculation](measure.py) · [Complete records](results.json)

| Probe geometry and law | Positive source: excess at fixed receiver | Negative source | Difference |
| --- | ---: | ---: | ---: |
| One opening, contrast 8 | Same reading | Same reading | Less than 2×10⁻¹¹ at every tested exposure |
| Two openings, contrast 8 | 0.406115753 | 0.428729753 | −0.022614000 |
| Two openings, contrast zero | 0.059842419 | 0.052143699 | +0.007698719 |

The receiving cell's native signed seam pointers resolve the local word and the complete load-0/host/load-2 outcomes at its fixed port. The contrast-zero law keeps the graph, alphabet, preparations and reader and reverses this particular detector contrast. All reported probabilities are ideal finite-instrument predictions with no sampled trials. The two-opening setting was exposed during exploration before the supplementary design. The maximum residual-based double-precision probability allowance is 1.76×10⁻¹⁰; independent raw-board propagation reproduces the headline readings.

## One opening gives an exact obstruction even during penetration

With only owner 6 beside source owner zero, the open bridge is (1,82). The source-cell permutation

```text
[9, 1, 2, 10, 11, 5, 6, 8, 7, 0, 3, 4]
```

fixes gateway 1, swaps the two closed gateways and preserves the internal graph and reference pattern. In a rational basis of the compact source space it acts as [[0,1],[1,0]], with determinant −1. Every imaginary antisymmetric circulation operator on this two-dimensional real space reverses under that reflection, so the two opposite-current sources are exchanged up to a global phase. The implemented source states satisfy that relation within 8.7×10⁻¹⁶.

Extend the reflection by the identity on the exterior owner and apply it to both defect positions. The complete one-opening generator commutes with it entrywise, at contrast 8 and zero. Every exterior-only instrument does too. Unitary evolution followed by an arbitrary exterior outcome therefore gives equal probabilities for the two sources; induction extends the equality to adaptive exterior protocols and untouched references. This null permits the probe to enter the source. Its scope is the declared two-owner aperture. Additional probes beside the other gateways can break the symmetry through retained contact energies.

## The positive reading includes source back-reaction

At contrast 8 and exposure 32, the complete source diagnostic on the exterior-probe branch gives:

| Outcome, all original launches retained | Positive initial current | Negative initial current |
| --- | ---: | ---: |
| Exterior probe and original source mode | 0.222108 | 0.235773 |
| Exterior probe and opposite compact circulation | 0.477882 | 0.504020 |
| Exterior probe and source outside the doublet | 0.114215 | 0.106596 |
| Probe inside the source owner | 0.185794 | 0.153611 |

These source-return tests use a supplied ideal coherent number flag and local mode analyzers. No bare source-mode projector is assigned when the excess occupies the source owner. Their complements are explicit, and no source-return postselection changes the exterior detector's denominator. The large opposite-mode return and probe penetration distinguish this result from an unchanged stationary source acting through a static exterior field. The calculation establishes intrusive material scattering. It supplies neither a native static dipole field nor an electron gyromagnetic factor. The full local OW-05 criterion is partial.

## Replay

From the RER root, use single-thread BLAS. Replay preserves the receipt. `--record PATH` explicitly creates a first record and refuses an existing path; serialization completes before the destination opens.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_05_dipole.penetrating_probe.measure
OPENBLAS_NUM_THREADS=1 python -B -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_05_dipole/penetrating_probe
```
