# Does one kinetic coefficient survive a change of probe?

The same neutral moving band has a positive, reproducible low-shell response in operational graph and clock units. After correction for independently measured owner populations, one calibration transfers to the other two launch sources and the second low-shell profile. It does **not** transfer to the third spatial shell. This rejects a particular scalar kinetic interpretation even when arbitrary on-site inhomogeneity and every transition line in the moving band are allowed. It strengthens the earlier dominant-line test; it supplies no universal no-go for composite inertia. [Complete counts and certificates](results.json) · [Derivation and limits](PROOF.md)

The scalar hypothesis is H = A L + V on the independently specified 20-owner graph, A > 0, with arbitrary real diagonal V. Its resting envelope is nodeless and its point-owner density equals the measured coarse density. The sum of all energy-weighted density transitions then obeys M₁(f) = A Σ_edges √(pᵤpᵥ)(fᵤ−fᵥ)². Independent occupation counts determine p. No single line, homogeneous population or zero potential is assumed. Native measurements use the actual arithmetic constituent-centre grating and a random constituent-owner population read. Identifying these with one point-owner algebra is part of the tested dictionary, not a theorem about an extended composite.

| Source and spatial profile | Inferred A, native exchange energy × owner-step² |
| --- | ---: |
| Calibration: launch 0, shell 1, profile column 0 | [2.0871, 2.1519] × 10⁻⁶ |
| Held: launch 5, same profile | [2.0924, 2.1468] × 10⁻⁶ |
| Held: launch 10, same profile | [2.1003, 2.1385] × 10⁻⁶ |
| Held: launch 0, shell 3, column 0 | [1.6845, 1.7254] × 10⁻⁶ |
| Held: launch 5, shell 3, column 0 | [1.6886, 1.7212] × 10⁻⁶ |
| Held: launch 10, shell 3, column 0 | [1.6948, 1.7149] × 10⁻⁶ |

Intervals include simultaneous 99% synthetic sampling allowances and conservative numerical source/spectral bounds. The held third-shell intervals remain disjoint after adding a further 5% tolerance to the calibration. The full receipt retains all six profiles and both clock settings: it records other compatible and rejected settings as well. The calibration would give m = 1/(2A) in [232348, 239572] native exchange-time / owner-step², conditional on that scalar interpretation. The failed extension prevents assigning this number as a single mass across the tested profiles. No metre, second or kilogram calibration is inferred.

## What the observer reads

Use the same original sharp excess–vacancy sources and select the lowest unresolved doublet of the existing 24-line moving band. Admissions are 0.012244%, 0.026168% and 0.306870%. This is the same rest selection as the earlier recoil experiment. On separate fresh starts read both terminal addresses and choose one constituent with a fair classical bit, giving the owner population. The six grating profiles are normalized projector columns 0 and 10 of the first three independent owner-Laplacian shells. Their definitions use only the graph.

For each grating, a coherent timing reference compares applying it before or after a common free flight. Recombine in Y, retaining its two exits, grating refusal, final moving-band refusal and initial rest-filter refusal. Flights of 4096 and 8192 exchange units sample the summed sine response. A bandwidth inequality bounds the full energy first moment without resolving or discarding any line in the 24-line band. Both arms have the same elapsed time, so a uniform energy shift cancels. Each count bank enrolls 2⁵⁸ original starts; no failed final filter is removed from the score denominator. [Acquisition and all observer formulas](measure.py)

An independent native 10/01 dimer clock has gap one: its quarter-unit load-read calibration and held half/three-quarter appointments agree. A parked excess compensates this clock input’s vacancy. Its reference response fixes the operational exchange unit; protected scheduling and long-time tick repetition remain supplied. The native generator, source filter and full-space residual bounds are recomputed. The coherent signed amplitude grating, ideal spectral filters, timing reference, controlled application and terminal sharp reads remain specified apparatus. Their native compiler is not supplied here. The labels “held” and “calibration” describe the split inference: existing spectra, earlier failures and a deterministic pilot were known when this exploratory test was designed. The 5% tolerance and chosen settings are not a prospective registration or an independent experimental discovery.

## Controls and scientific scope

Uniform and inhomogeneous scalar reference particles both pass all held native-rate profiles using one reference calibration. The latter has a deliberately supplied on-site potential and a nonuniform ground envelope, so population correction and multiple-line summation receive a positive check. These references validate the reader; their assigned kinetic coefficient earns no emergence claim.

The changed native law halves every free-flight generator rate with carrier, alphabet, routing, source and reader held fixed. All its measured coefficients reject the original native-rate calibration. Reading the timing arm destroys interference and returns a null response, with its complete failures retained. The finite energy shifts and timing fringe, rather than simulator state labels, supply the observer's kinetic response.

This supports a source-independent positive response on the tested first-shell profiles of a stable load-neutral band. The existing native material-source experiment separately tests neutrality on that band. A scalar point-owner kinetic law across the larger spatial range fails. Matrix-valued internal dynamics, nonlocal kernels, other centre coordinates and a momentum-resolved long-wavelength branch remain possible. A source-independent physical mass and electromagnetic identification remain open, so the full literal OW-19 assessment stays locally partial and officially untested. [Native charge measurement](../native_source_probe/README.md) · [Full-address recoil](../unfiltered_recoil/README.md)

## Replay

From the RER root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.sum_rule_inertia.measure
.venv/bin/python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.sum_rule_inertia.readback
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -W error -m pytest -q model/experiments/material_inertia/pair_spectrum/sum_rule_inertia
~~~

The headline recomputes the 57,360-configuration source sector and every count, leaving the receipt read-only. `--output /tmp/new-sum-rule.json` writes a new exclusive acquisition. The solver-free [readback](readback.py) reconstructs all intervals and verdicts from count arrays, geometry and the explicit numerical certificates. Four tests independently check the graph sum rule with an arbitrary potential, the complete reference interferometer, count inference and held controls, and a resolvable held-fringe mutation plus optimized-mode refusal.

Verified with Python 3.12.6, NumPy 2.3.5, SciPy 1.18.1 and NetworkX 3.7 on the existing macOS environment with single-thread BLAS. The final normal and optimized full replays exited zero in 9.82 s and 11.37 s; four tests passed in 2.37 s. Count-only replay, Ruff, source syntax, whitespace and local links also pass. Source receipt SHA256: `5664837c5d206c70193508ef927b09a4469161bfa9887775ebc2c587a34bc47e`.
