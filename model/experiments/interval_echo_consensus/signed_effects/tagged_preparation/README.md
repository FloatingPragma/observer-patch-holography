# Preparing the common state on the material board

A single native source contact and a neutral register prepare the same reduced state for every context of the eight-slot phase table. The contact duration is calibrated to the older diagonal fibre census, which counted 111 and 68 records in its 179-member block, within a fibre of 512 and a total mass of 16,384. This experiment uses a fresh material source on the 20-cell Interval Echo board. The connection to the historical run is the recorded ratio. [Original census](../../../../../code/born_context_phase_lift/BORN_CONTEXT_WEB_PAYLOAD.v1.json) · [Archived context operands](../../../../../code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json)

Start with exactly one vacancy or excess at port 0 in a unit-filled data cell. Port 239 supplies a complementary charge reference so that superpositions of the two signs have a common total load. Isolate the contact between ports 0 and 1 for an exchange time of 2 arccos √(111/179), about 1.32820. For either sign the two address amplitudes, up to a common phase, are √(111/179) and i√(68/179). A native sequence copies occupation away from one at port 1 into a neutral (2,0) register at ports 76 and 77, using the unit host 79. Three π-contact swaps route port 1 along 1–82–72–73 and three reverse swaps restore the transit ports. Both orientations of the tag and the two source signs carry the same copy phase. The preparation code exponentiates the native contact on its full local alphabet, checks for seam-memory leakage and evolves all ternary occupations of the active ports. [Source implementation](source.py) · [Calibration and coherence tests](test_source.py)

Leaving the tag unread gives the data density diag(111/179,68/179) for either sign. Skipping the copy gives the same two diagonal weights and an off-diagonal amplitude of −i√(111·68)/179. Each real projector in the table has the same probability in both preparations. The +Y projector distinguishes them. The source, pulse isolation, clock, clean pointers and terminal probability rule are supplied. This construction does not identify the old fibre trajectory with the material source or derive the probability rule.

## The joint acquisition

For every named context, first apply that context's calibrated material operations on the source rails. For the weak version, address the extra rails, erase the unused record bit coherently, then read the four pointer words and conditionally restore the source. The experiment propagates both tag columns through each of the 17 full native output branches. Squaring the resulting amplitudes gives 34 joint terminal bins for each setting, including unused address words. The analyzer word avoids both tag ports, the tag host and the charge reference. [Instrument implementation](../phase_table/measure.py) · [Joint acquisition](measure.py) · [All 64 simulated cohorts](results.json)

| Named slots | Tagged sharp yes | Copy omitted, sharp yes | Tagged weak yes | Copy omitted, weak yes |
| --- | ---: | ---: | ---: | ---: |
| Diagonal and conjugated 0–1 | 0.620112 | 0.620112 | 0.560056 | 0.560056 |
| Conjugated 2–5 | 0.439944 | 0.439944 | 0.469972 | 0.469972 |
| Phase (+Y) | 0.500000 | 0.985359 | 0.500000 | 0.742679 |

The probabilities apply to either sign. For example, the vacancy +Y instrument recorded 130,728 sharp yes outcomes with the tag and 258,097 without it, each from 262,144 independent simulated launches. The weak instrument recorded 130,843 and 194,769. Across all 64 cohorts, the largest absolute difference between a sampled yes frequency and its circuit probability is 0.002856. The simultaneous 99% Hoeffding half-width is 0.004248. The counts are fresh model draws. The historical experiment supplies the calibration and comparison targets.

From the RER root, replay with `OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.signed_effects.tagged_preparation.measure`.
