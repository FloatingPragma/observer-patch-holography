# Do the joint load–coin phases give a matter-wave dispersion?

The terminal phase reader detects load–coin coherence after free flight. Across five spatial modulations and two sharp releases, **none of the 20 responses** had a single frequency that explained 90% of both training and held variance. The imaginary component was read on all 630 seams of the sphere. A native reference phase also exposed the real component on 592 seams. Both quadratures passed separate sensitivity controls. [Protocol](PHASE_SPECTRUM_PROTOCOL.md) · [Counts and frequency fits](phase_spectrum.json)

| Squared graph wave number | 32-state current | 32-state real bond | 2,730-state current | 2,730-state real bond |
| ---: | ---: | ---: | ---: | ---: |
| 0.118992 | 0.380 | 0.375 | 0.365 | 0.105 |
| 0.322744 | 0.662 | 0.468 | 0.381 | 0.414 |
| 2.897020 | 0.141 | 0.199 | 0.209 | 0.193 |
| 6.000000 | 0.362 | 0.375 | 0.383 | −0.140 |
| 6.846204 | 0.179 | 0.237 | 0.184 | 0.375 |

Entries are the fractions of held variance explained by the fitted frequency, offset and two quadratures. The fits use exposures 0 through 16 and predict the following 384 settings through exposure 32. The best noiseless held fractions in the four columns are 0.659, 0.464, 0.382 and 0.415. Each address and terminal contact used 16,384 simulated preparations, including quiet. The two-position control recovered the native current frequency of 8 at 7.999812 with held explained variance 0.999963. Starting from a quarter oscillation and preparing opposite reference phases gave real contrasts of −1 and +1. Reading the coin before mixing erased both contrasts at the probability level.

For each physical edge the terminal contacts interfere the two accessible configurations that differ by moving the load and reversing that edge coin, with all spectator records matched. They measure the real or imaginary part of their joint coherence, summed over the spectators. Each scalar graph harmonic defines a spatial gradient for the imaginary read. The resulting current obeys the independently checked continuity equation for that harmonic's density. For the real read the same harmonic modulates the bond midpoint, counting the 592 edges for which a degree-five endpoint has a physical phase reference. Its scalar wave number describes the modulation of the reader. The resulting edge profile need not itself be an eigenmode of free matter.

The sharp source can excite several energy levels simultaneously. Both quadratures are bilinear in the state, so their frequencies can include differences between levels. Identifying an isolated matter-amplitude branch would require a mode-selective preparation and an independently calibrated phase reference on that branch. This experiment uses native service exposure for the frequency scan. Its controls establish sensitivity to coherent beats, while the tested spatial responses supply no rest gap or speed for a Klein–Gordon fit.

The receipt stores each click tensor in exposure–link–setting order. Decode the `clicks` field with Base64 followed by zlib decompression, interpret the bytes as little-endian unsigned 16-bit integers and reshape to the stored dimensions. Quiet is the difference between 16,384 and each click count. To reproduce all tensors and fits from the repository root, run

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_arrows.phase_spectrum
```
