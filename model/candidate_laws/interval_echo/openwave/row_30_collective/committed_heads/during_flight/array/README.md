# What seven wall receivers hear after two connector readings

The observer in cell 2 reads the same connector coin at one and two exchange radians. On independently prepared copies, receivers on seven successive bonds to the right count unequal loads. Two localized launches give different records. Their signed difference oscillates at several rates. On independent counts the two strongest separated bins have angular frequencies **1.41126** and **1.10447** per exchange radian, with bin spacing 0.01227. Both seven-bond patterns exclude a single spatial wave number, even allowing an arbitrary reflected wave. These oscillations are density beats of the measured wall. Reading either as the energy of one freely traveling mode would place it far below the independently metered gap near 15. [Enrollment and wavelength test](PROTOCOL.md) · [Counts](counts.npz) · [Decoded results](results.json)

A third frequency distinguishes passage through the observed connector. At **0.353781** the held wall record on bond 11 has amplitude **0.010630**. For the nominal cut generator the complete spectrum of the nine right-hand dimers gives an upper bound **0.010103** there for any two incoming preparations. The joint count radius is **0.000164**. A supplementary bound using the independent contact-meter sensitivity box and the clock-spacing interval gives the stronger ceiling **0.006861**. Thus the response retains a contribution from the attached strand after both local coin readings. The reverse test, using the conservative nominal bound, identifies zero frequencies requiring the cut branch. The calculation pools every home record. It does not identify the passage of an individual launch from its terminal answer. [Derivation and bounds](SPECTRUM.md)

## What the stations record

On the commissioned 14-dimer strand the sources write sharp reversals of orientation at bond 4 or bond 1. The same head at ports 31–32 is read twice. Each receiver at bonds 7–13 records the two loads of its own connector on a fresh source copy. Across 4,096 clock appointments, the recorded alternatives are 00, 01, 10 or 11 for the two home coin readings, each with wall or no wall at the remote detector. Each of two independent count banks retains 524,288 simulated starts per site, source, appointment and law. The clock is a separate blank exchange dimer with 16,777,216 starts for each sine and cosine read. Across both contrasts, all seven sites and both replications there are **120,259,084,288 simulated wall starts**. They are seeded draws from computed joint probabilities. A native propagation of the initial two exchange radians supplies four amplitudes, and a quadratic-observable identity propagates them to the final read. [Material evolution](material.py) · [Sampler and complete exposure](measure.py)

| Quantity read or bounded | Result |
| --- | ---: |
| Clock spacing between neighboring terminal cohorts | 0.12500003 ± 0.00000096 exchange radians |
| Held amplitudes at the two discovery frequencies | 0.012181 and 0.011644 |
| Second singular value of the two complex spatial profiles | 0.011476 against a simultaneous error bound of 0.000614 |
| Largest signed remote contrast with connector strength set to zero | 0.003733 against a simultaneous bound of 0.005838 |

The amplitude is the modulus of a unit-sum Hann-windowed Fourier coefficient of the *difference* in hit rates, with all starts in each denominator. The earliest and latest clock readings calibrate the temporal scale. Uniform appointments and their ordering are supplied. The source injects no wavelength. An independent acquisition measures the spatial profile after temporal selection. The enclosure for a superposition of two opposite waves uses its spatial recurrence at every one of five overlapping triples. Both sets of allowed real wave numbers are empty. Each finite time window can include unresolved spectral lines. Coordinates count physical dimer links along this commissioned material. [Count-only decoder](readback.py) · [Checked recurrence](SpatialPair.lean)

OpenWave lists M10 as untested for the Klein–Gordon row. Local dispersion on a prepared strand with *unread* connections belongs to the [collective dispersion experiment](../../../README.md). This detector resolves propagation, multiple density-beat profiles and incomplete cutting under two local reads. Single-mode dispersion under that schedule remains open. The prospective public confirmation has its own frozen admission and acquisition conditions. [Official row](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L467) · [Filing decision](../../../../submission/ROW_30_DECISION.md)

Run the following from the repository root. The replay compares all raw counts, the clocks and the decoded report.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_30_collective.committed_heads.during_flight.array
python -B -W error -m "$p.measure"
python -B -m pytest -q -W error --pyargs "$p.test_array"
~~~
