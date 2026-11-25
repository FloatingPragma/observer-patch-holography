# Does the row 4 pair persist in the row 30 spectroscopy materials?

Neither of the two archived 20-cell materials reproduces the M10.4 home-detection benchmark across the nine filed settings. One material prepares the same coherent pair as M10.4. At strength ½ and exposure 37 the original observer detects that pair on a home seam with probability **0.086409**, compared with **0.995945** on the eleven-write M10.4 material. On the second material a written contact changes the prepared state. Even supplying the M10.4 state directly gives only **0.311740** at that appointment. [Protocol fixed before the flights](PROTOCOL.md) · [All terminal records](results.json) · [M10.4 review](https://github.com/openwave-labs/openwave/pull/749#pullrequestreview-6060088147)

## Selection and observable

Row 30 recorded two fair local readings on each empty seam and wrote a neutral memory when both were positive. Its one-load experiment used the resulting material to measure spectral bands and transport. We used both archived histories on the board specified by M10.4 and reproduced their packed coin words with the original formation routine. The materials were archived on November 24, 2026, one day after the M10.4 source pin. A previous geometry check had established the absence of the eleven-write motif in these materials. This chronology establishes that the donor records preceded this pair flight. [Archived formation rule and measurements](../../row_30_formed/README.md) · [Donor revision](https://github.com/FloatingPragma/observer-patch-holography/tree/44104ae607c93382f6115181337e26846ac868d9/model/candidate_laws/interval_echo/openwave/row_30_formed)

The two loads start on ports 156 and 157. The pulse schedule, retry budget, home cell and thirty detector addresses come from M10.4. Both materials leave the center's meter seam blank. Copying its nonzero contrast after echo and undoing the echo measures the same center-occupation projector used by the filed neutral meter. The receipt includes all thirty home addresses and quiet on every trial. Quiet includes exhausted preparations and any accepted pair outside those thirty seam configurations. It leaves open whether the loads are still adjacent elsewhere. [Native meter check, acquisition and full-board flight](experiment.py)

## Results at exposure 37

Each probability in the two middle columns uses the full initial ensemble. The final column is the maximum over any normalized coherent input on the two possible accepted configurations (156,157) and (156,166). It allows arbitrary amplitudes, phases and perfect preparation. It bounds mixed inputs as well.

| Row 30 bank | Strength | Pair prepared on this material | Filed signed pair placed on this material | Maximum over the two input configurations |
| --- | ---: | ---: | ---: | ---: |
| 178121 | 0 | 0.116794 | 0.116794 | 0.131257 |
| 178121 | ½ | 0.086409 | 0.086409 | 0.120220 |
| 178121 | 8 | 0.611809 | 0.611809 | 0.636323 |
| 178122 | 0 | 0.360607 | 0.307775 | 0.369910 |
| 178122 | ½ | 0.448467 | 0.311740 | 0.580730 |
| 178122 | 8 | 0.998243 | 0.884497 | 0.998267 |

**Bank 178121** has 167 writes. Both pulse contacts transmit. The three occupied readings prepare the filed signed state, including the original 0.003864 exhaustion probability. Only one of nine appointments clears 0.94. Even at strength 8, which row 30 used for spectroscopy, the exposure-37 return for any input in this two-configuration source is bounded by 0.636323.

**Bank 178122** has 147 writes. Its (156,166) contact is inactive, so the pulses leave the initial pair on (156,157) and all preparations are accepted. This sharp input clears 0.94 at all three strength-8 exposures. At the same strength the historical signed input gives 0.869581, 0.924639 and 0.884497 at exposures 1, 7 and 37. At strength ½ and exposure 37, the maximum is 0.580730 even with perfect state preparation.

The upper bound is the larger eigenvalue of the two-by-two matrix formed by propagating both input configurations, restricting their amplitudes to the thirty home seams and taking their Gram matrix. The same-material stochastic control on the filed input returns 0.044959 and 0.129537 for the two banks at exposure 37. The finite results rule out the historical 0.94 all-strength, all-exposure claim on these two materials at the original address. Persistence on the filed eleven writes and the upstream partial award stand. [Filed connected-pair evidence](../README.md) · [Complete comparison records](results.json)

## Reproduction

From the repository root, replay both material histories, the local meters, the full flight and the archived bins.

```sh
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04.commissioned_medium.experiment --verify
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_04/commissioned_medium/test_transfer.py
```

The calculation reproduces all nine filed eleven-write distributions within $8\times10^{-12}$ per outcome and the matched blank-material return of 0.440929105 at strength ½ and exposure 37. On each donor material and at each strength, a separate sparse matrix exponential checks every home-seam effect on the two-dimensional input and the full-state norm at exposure 1. The second command builds the generator directly from the board edges and packed coin words. It checks both exposure-37 counterexamples at strength ½, including the thirty-one bins and the maximum over the two input configurations. The record also contains the SHA-256 digests of the fixed protocol, the board, the donor archive and both packed histories.
