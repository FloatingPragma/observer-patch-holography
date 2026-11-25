# Spatial charge records of the moving neutral pair

After the same filtered pair flies for one million exchange units, two charges occupy different cells in about **0.107%** of admitted runs. Complete native signed records determine each cell's charge. There is no arrival or common-cell herald. A legal hard-core model comparison with matching signed port marginals gives **29.33% to 50.54%** separated-cell probability. Four further exchange-only units, from the same arrived state, raise it to **61.03% to 61.93%**. Native continuation preserves the small value. [Protocol](PROTOCOL.md) · [Code and native decoder](measure.py) · [Complete probabilities](results.json)

| Original launch owner | Native arrival: different cells | Matched-marginal model | Exchange-only continuation |
| --- | ---: | ---: | ---: |
| 0 | 0.001069879 | 0.505369309 | 0.610348177 |
| 5 | 0.001070219 | 0.293308482 | 0.616889457 |
| 10 | 0.001068564 | 0.344084054 | 0.619293066 |

These are ideal probabilities conditional on the original filter, whose admission is about 13.27%. Every cell record also includes the filter refusal in its original-start distribution. The source's spectral and propagation error bounds are retained per launch. The control mixture forbids coincident excess and vacancy, preserves the computed one-point marginals within 2×10⁻¹⁴, and is a prepared model comparison with no native reachability claim. Its effect separation exceeds the source uncertainty by more than four orders of magnitude.

## What the field probe measures

The native meter reads all 30 internal signed seam pointers by echo, observation and undo. On the declared one-excess/one-vacancy domain their complete records resolve 157 allowed local load words and determine charge Q_C ∈ {−1,0,+1}. All three charge outcomes are retained. The sum of their uncentred second moments is Σ_C⟨Q_C²⟩ = 2 Pr(different cells). This is a joint measurement, while the summed signed mean vanishes for every globally balanced pair. The microscopic port second moment is two by the input sector, so the suppression concerns the specified cell scale.

A separate Ramsey qubit with a supplied phase coupling θQ_C gives X and Y readings from the characteristic function of these native charge records. Host and single-charge calibration, θ = π/2 and held π/3, complete plus/minus laws and filter misses are recorded. A finite coherent sample–qubit unitary checks that map. The qubit, its phase drive, coupling and clock are supplied apparatus. The calculation derives no electromagnetic field or physical electromagnetic coupling. It tests spatial neutrality under that declared probe and distinguishes the bound source from an uncorrelated charge mixture, with every final address retained.

The same 24-line source and original launches underlie the motion and recoil studies. The different-cell comparison supports a locally charge-balanced moving composite at this resolution. Physical electrical neutrality, a field interaction derived from the native law, and a calibrated inertial mass are separate obligations. The local full OW-19 assessment is partial; this experiment supplies no platform award.

## Reproduction

Run from the RER root with single-thread BLAS. Replay recomputes the band, all native charge laws, both probe phases and the changed-law continuation, then compares with the receipt. `--record PATH` explicitly creates a first record and refuses an existing path. Serialization completes before the exclusive destination opens.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.cell_field_probe.measure
OPENBLAS_NUM_THREADS=1 python -B -m pytest -q -W error model/experiments/material_inertia/pair_spectrum/cell_field_probe
```
