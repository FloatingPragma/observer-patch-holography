# Free dispersion of the doubled load

On the 20-cell sphere, five separately prepared waves of the same doubled excitation used on the guide gave a beat frequency of one half the spatial Laplacian eigenvalue. The spatial measurement uses signed mode profiles reconstructed from local repairs. A second set of preparations supplies the time series, read on three seams of one cell. All 630 seams evolve freely during flight. [Protocol](PROTOCOL.md) · [Measurements and raw clicks](results.json)

| Spatial eigenvalue from counts | Frequency from counts, radians per exposure unit | Held-time RMS error in port occupation |
|---:|---:|---:|
| 0.119026 | 0.059498 | 0.0000203 |
| 0.322777 | 0.161365 | 0.0000145 |
| 2.897034 | 1.448504 | 0.0000116 |
| 6.000002 | 2.999999 | 0.0000500 |
| 6.846201 | 3.423103 | 0.0000160 |

For the unweighted incidence definition of squared wave number, the exact free generator gives a Schrödinger dispersion in this sector. Adding any finite rest frequency produces a quadratic term in the squared frequency. That curvature prevents an exact Klein–Gordon fit with one common speed. The finite-count data alone can approximate the large-rest-frequency limit and require an independent mass or speed calibration to bound it. The result concerns a single doubled load, with service exposure as the time coordinate. The alternating guide has its own [dispersion experiment](../row_30/PROTOCOL.md).

From the repository root, reproduce the receipt with

```sh
python -m model.candidate_laws.interval_echo.openwave.row_30_free.experiment
```

The default run checks the archived JSON. Passing `--output PATH` writes a separate receipt to a fresh path. [Experiment](experiment.py)
