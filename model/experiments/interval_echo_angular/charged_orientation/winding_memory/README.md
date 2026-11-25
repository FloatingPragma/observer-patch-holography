# A vacancy writes a coherent winding record

A vacancy propagates across five cell boundaries in a guide whose contacts all carry unit-contrast marks. Returning the vacancy to its starting port can reverse the direction on every head. Two preparations with the same initial material positions and head populations give different arrival counts at a foreign observer. The observer needs only a local vacancy detector to read the difference.

| Preparation and read | Native probability | Simulated counts per 32,768 launches |
|---|---:|---:|
| Sharp vacancy at port 3, free flight to time 2, source vacancy present and winding detected | 0.14157749 | 4,722 |
| First winding phase, free flight to time 1, vacancy at foreign port 25 | 0.00943503 | 309 |
| Opposite winding phase, same appointment | 0.20916356 | 6,811 |
| Incoherent mixture, same appointment | 0.10929929 | 3,619 |

The two phase preparations occupy the same launch port. Their amplitudes span two head words that differ around the entire guide. A separate experiment prepares a phase on a **single launch contact**. At time 1.5 the remote observer accepts a vacancy on ports 25 or 32 with probabilities 0.113144 and 0.204743 for the two inputs. Respectively 90.90% and 19.85% of those accepted vacancies appear at port 25. Applying a local half contact before the read changes the first-port probabilities per original launch from 0.102853 and 0.040641 to 0.024170 and 0.021340. Every appointment includes a no-arrival outcome. An exact rational enclosure certifies that the two-position transfer matrix has determinant magnitude greater than 0.1. Selecting either single remote port gives a one-dimensional received state.

This experiment identifies a jointly propagating vacancy and direction-head signal. Its two receiver states occupy different material positions. On returning to the launch port, the two winding histories carry distinguishable records around the guide. A clockwise circuit with every forward hop confirmed takes the unwound word to the wound word. Another such circuit from the wound word has zero success probability. The physical spatial rotation and exchange sign of a charged lepton require separate measurements.

## Preparation and read

The 15 physical ports lie on the shortest cycle of five adjacent owners of the 20-cell federation. The signed material is one missing unit in a host of ones, balanced by an excess at port 141. Each ring contact receives its native unit-contrast write in even reversal parity. A blank half swap and a phase dwell on a borrowed private contact prepare the two ordered load inputs in phase. After the echo the measured vacancy endpoint heralds the direction bit. A full native transfer corrects the opposite herald. A blank bypass returns the vacancy to the stock after each write. Five zero private writes and 55 zero exterior writes close the guide. The code checks exterior paths, the compensation, the final shutter and both herald corrections. These steps are addressed apparatus operations.

Every contact then runs under registered Interval Echo at its registered coupling. The outside marks decouple the guide, and common idle energies cancel. Each observer reads the local ternary load with the reversible two-reference echo meter. To read the winding, first flag a source vacancy, pulse the first guide contact, and read the source vacancy again. The unflipped head transfers the charge to the next port. The flipped head leaves it at the source. Sharp addressed clockwise and reverse circuits and the inverse splitter provide the preparation and readback controls. The free runs use none of those pulses between launch and the appointed read.

The receipt includes all seven appointed times, independent seeded binomial samples, and marginal exact 99% count intervals. The probabilities come from the native matrix. The sampled counts are simulation records.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.winding_memory.measure --write
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/winding_memory
(cd Lean && lake env lean -j 1 Dynamics/WrittenHeadWinding.lean)
```

[Mixed-guide formation lines](MIXED_MASSES.md) · [Finite derivation](PROOF.md) · [Observer records](records.json) · [Experiment](measure.py) · [Lean source](../../../../../Lean/Dynamics/WrittenHeadWinding.lean)
