# A second reading during a wall's flight reverses the remote source contrast

An observer reads the same connector coin twice while a wall is moving. A receiver in a different cell reads only the two loads on its wall detector. At its six-radian appointment the receiver recorded 948,779 versus 877,062 wall hits for the two launches. A single coin read during the same flight gave 568,489 versus 874,065. A second reading changed which launch was more likely to produce a distant hit. Every figure is a count from 16,777,216 seeded simulated starts per setting. [Frozen choices](PROTOCOL.md) · [Complete joint counts and computed probabilities](records.json)

## What the two observers can compare

The experiment uses the 14 dimers and both sharp source preparations from the [original connector experiment](../README.md). The first setting launches from bond 4 and the second from bond 1. Owner 2 reads the coin of the magnitude-two head at ports 31 and 32 after one exchange radian. On the two-read menu it repeats that read after another radian of free flight. Owner 3 reads for unequal loads on ports 42 and 43. On matched copies, the cosine and sine counts of an independent blank exchange dimer resolve the ordered appointments. The widest displacement of a simultaneous clock interval from its assigned appointment is 0.00140 radians. The clock reads certify the exposures of this apparatus; scheduling and the projective coin instrument are supplied.

The source contrast is the first launch's wall-hit fraction minus the second's, with every head history retained. At four radians the twice-read contrast is **0.08428–0.08499**, compared with **0.09344–0.09418** for one reading. It remains positive when each of the four coin histories is examined separately. At six radians the twice-read contrast is **+0.00378 to +0.00477**, compared with **−0.01865 to −0.01777** for one reading. The second local record therefore changes the shape of the remote signal. At the instant of the second reading the difference between the two menus' remote contrasts is consistent with zero, as a complete local operation requires. The intervals have joint 99% coverage under the simulated finite sampling law. [Count-only decoder](readback.py)

The connector's own readings distinguish the launches too. In the four-radian cohorts its two coin labels disagreed on 441,950 starts for bond 4 and 1,864 for bond 1. Repeating the read immediately gave the same label every time. A magnitude-only reading, which leaves the commissioned coin coherent, gave a remote contrast of 0.18174–0.18271 at four radians. With contrast removed, every remote source difference contains zero and there are no coin changes. The archived six-radian conditional results include unresolved rare-word differences. The all-start contrasts do not select a particular coin history.

## What repeated readings do to the material

The native local generator has two coin-parity sectors. In one the connector couples the neighboring loads, and in the other it contributes a constant energy. A single sharp coin reading interrogates the complementary basis. The probability that a second reading changes depends on the relative phase built up by the two intervening matter evolutions. The calculation uses the same energy origin for both. A direct four-dimer check constructs the entire coin and sign matrix from native 99-state columns and resolves the sharp outcomes in succession. It independently reproduces the selective calculation. [Native propagation](measure.py) · [Independent matrix and reduced-channel checks](test_inflight.py)

For a predetermined schedule at this head, pooling complete coin records gives the composition of equal mixtures of these two unitary evolutions. In the limit of frequent readings the matter dynamics retain half the original interaction at that connector. For a fixed source, an event probability differs from its frequent-read limit by at most the smaller of one and κ² times the sum of squared inter-read durations. This follows for every terminal detector and every source density matrix in the stated sector, including sources with a retained reference. The finite experiment resolves two time-ordered readings. [Channel identity and finite bound](FREQUENT.md)

**OpenWave status.** OpenWave lists M10 as untested on OW-30. This experiment establishes remote wave transport subject to two scheduled local observations and shows a change in its measured time dependence on a prepared one-dimensional material. A dispersion for this read schedule, the source and decoder admission, and the prospective confirmation remain open. [Official cell](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L467) · [Filing decision](../../../submission/ROW_30_DECISION.md)

Reproduce the simulated acquisition from the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_30_collective.committed_heads.during_flight
python -B -m "$p.measure" --write
python -B -m pytest -q -W error --pyargs "$p.test_inflight"
~~~
