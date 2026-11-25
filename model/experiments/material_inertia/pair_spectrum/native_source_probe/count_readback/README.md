# A charge bound read from the material meter counts

The two source-clock settings measure how often a selected region contains an uncompensated excess or vacancy. After one million exchange units, the strongest simultaneous 99% upper limit across all 60 regions and all three native launches is **0.1031%**. Each regional bound applies among launches passing both the original spectral filter and the unit-filled meter-site test for that region. The original band filter admits about 13.27% of launches. The result measures signed load on this board. Its relation to electromagnetic charge and physical mass remains open. [Original counts, including refusals](../results.json) · [All inferred intervals](results.json)

| Pair launched in cell | Largest native upper bound, any region | Native charge fraction at launch cell | Matched-marginal control at same cell | Four more units of exchange alone at same cell |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.1031% | 0.0831–0.0938% | 46.02–46.24% | 47.06–47.28% |
| 5 | 0.1005% | 0.0854–0.0963% | 27.23–27.42% | 49.24–49.46% |
| 10 | 0.1003% | 0.0822–0.0929% | 31.63–31.83% | 46.45–46.67% |

The probabilities in the source simulator are never read by the inversion. For each measured X setting, divide the count of dark pointer word 2 by the sum of all nine eligible meter-word counts. Refused and filtered launches stay identifiable as the other five bins. The registered dark fraction for each of the two compensated sharp charge standards and the neutral host is measured separately at each dwell. Conditional on eligibility, the dark fraction is a convex combination of these three calibration probabilities since the source coupling is diagonal in regional load. The signs may have different meter efficiencies.

For each count bank invert the one-sided binomial Chernoff bounds in both directions using a tail allocation of 0.01/(2 × 1,086). Intersect the two dwell estimates for the same region and preparation. The 1,080 sample banks and six standard banks have simultaneous coverage of at least 99% under the simulated experiment, including the randomness of the number of eligible trials. Dependence between the reported intervals from shared calibrations is allowed by the union bound. The numerical spectral and finite-flight error allowance for modeling the native source is reported separately in the original receipt. The intervals concern individual regional measurements on fresh starts. They do not imply a simultaneous observation of all 60 regions on one launch. [Count-only implementation](measure.py)

## Including every meter refusal at the launch cell

At each original source cell, a second count-only bound covers the **entire spectrally admitted band**. The signed-charge event commutes with the meter-site eligibility check, so every refused member can conservatively be assigned any charge. If the refused fraction has upper bound f and the eligible fraction's charge lies between l and u, the full-band charge lies between (1−f)l and f+(1−f)u. A separate 99% simultaneous family allocates the error across six calibration probabilities, six eligible dark fractions and six in-band refusal fractions. Intersecting the two dwell bounds gives **0.0840–0.0929%** for cell 0, **0.0864–0.0955%** for cell 5 and **0.0832–0.0920%** for cell 10. At some other placements up to 39.47% of the band fails the meter-site check; the all-region limits above retain their stated eligibility condition.

From the RER root, run `.venv/bin/python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.native_source_probe.count_readback.measure`.
