# A cage in the flat branch of the scalar channel

Mark four seams that form a ring on a cell with plus, minus, plus, minus. At every port of the ring one plus
meets one minus, so the port passes nothing on to its other seams. An object on the ring can move between the
ring and the four ports and nowhere else: the signs of the law make a cage. This package builds that object on
the [scalar channel](../../candidate_laws/scalar_channel/README.md), law 2 of the
[construction manual](../../construction_manual/CONSTRUCTION_MANUAL.md) (sections 5.9, 6.6, 7.4 and 7.5), and
executes three experiments of an observer on it.

The candidate is `A-native-caged-exciton-v1`: the free law of the scalar channel is unchanged, and the
observer's menu gains a potential and a converter on the profile of the ring. The scalar channel is a candidate
under test, like the [contrast-memory law](../../candidate_laws/contrast_memory_spherical_v1/README.md), and the two laws act on different
state spaces.

## What is shown

| Statement | How | Result |
| --- | --- | --- |
| The law keeps the pair of a ring and couples it to nothing else. | theorem ([PROOF.md](PROOF.md), Proposition 4), in integers on one cell and on every ring with five seams at its ports on 20, 80 and 320 cells | 92, 312 and 1,272 such rings |
| The excitation of a cage is kept for every interval and every start. | theorem (Theorem 9) | the probability in the state is the same at five waits to rounding |
| A local detector reads the excitation. | three preparations, five waits, four analyzers, counts | sharp pulse: 0.9999 to 1.0002; coherent pulse: 0.7059 to 0.7060; no pulse: 0.000004 |
| After a finite pulse the seam record swings at one line with no drive. | the observer reads its own reference and the seam record in one trial, counts | 4.898982 and 4.898976 of the reference (standard error 0.000002) for steps of equal and of unequal length; the law gives $2\sqrt6 = 4.898979$ |
| The level difference of the pair, read separately, is the same number. | a fuel whose level is set by a dial, counts | 4.898980 (standard error 0.000008) |
| A changed law with the same alphabet, routing and symmetry loses both. | seam hopping round every port, counts | no single line (the best one carries 0.51 of the variance); its counts lie up to 0.56 from the record law of the law |
| The two operators are words on the 432 registers. | identities of integer polynomials | $n(e)$: 17 words; the converter: 8 words; both commute with all 260 loops |

Controls that test the readers: no pulse (amplitude 0.0019), the two modes of the cage stopped during the wait
while the reference runs (below 0.00001), a sharp excitation (no swing), the converter switched off (no record of a
full fuel), and counts that count against each claim ([test_cage.py](test_cage.py)).

Every reading is computed from counts. [counts.npz](counts.npz) holds them: 100 million trials per setting of
the clock and of the detector, 4 million per setting of the dial, drawn once from the record laws with a fixed
seed. A replay recomputes the record laws from the law and the prepared reference, checks that the stored
counts are draws from them, and recomputes every reading.

## Run

From the repository root, with `model/requirements.txt` installed:

```sh
python -B -O -m model.experiments.cage          # recompute and compare with results.json, a few seconds
python -B -m pytest -q -W error model/experiments/cage
```

| File | Content |
| --- | --- |
| [construction.py](construction.py) | the incidences, the rings, the pair, the changed law, the federations |
| [instruments.py](instruments.py) | the record laws of the three experiments, and the readers that take counts alone |
| [registers.py](registers.py) | the two operators as words on registers, with their identities |
| [replay.py](replay.py) | the recomputation against [results.json](results.json) and [counts.npz](counts.npz) |
| [PROOF.md](PROOF.md) | the statements with their proofs |

## Supplied, and where the results end

Supplied: the two coefficients of the law, the prepared reference, the two operators on a profile with their
guard, the pulses, the fuel dial, the observer's reference with the bounds on its step.

The pulse supplies the phase and the work of the swing; the reference by itself does not swing. A cage has a
flat level and does not move. The law has no term between two objects. A swing that starts with no supplied
phase and the motion of a caged excitation are work in progress.
