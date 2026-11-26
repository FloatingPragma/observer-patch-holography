# Does a localized vacancy tick at its rest energy?

A vacancy released at one sharp port has a persistent local occupation beat on the 80- and 320-cell spheres. A nearby observer can isolate its frequency, **√5** in exchange-clock units, by monitoring the outgoing gateways before release and counting the vacancy at either end of a diameter. A separate energy reading on those localized packets gives **34.62 ± 0.02** and **34.58 ± 0.02** in the filled-host reference. The stationary vacancy calibrated by transfer from six-seam storage gives approximately 35. Even the largest possible single-vacancy beat on these boards is at most 22. OW-03 remains open. [Protocol](PROTOCOL.md) · [Registered experiment](measure.py) · [Sharp free launch](sharp.py) · [Observer counts](results.json)

Two other preparations separate timing from transport. A stationary neutral seam has a free beat at its independently measured excitation gap. On an empty host a pair of unequal loads carries the phase set by a local pulse to each of three receivers. The remote phase drifts relative to the isolated clock on all three bridges. [Stationary gap clock](written_coin/README.md) · [Moving pair and energy measurement](moving_pair/README.md)

## What starts the motion?

Put one vacancy at port 4 of an otherwise filled cell and open every contact immediately. The origin detects the missing load using the signed comparison to a filled neighbour; distinct signs resolve either end of that contact. The sharp initial record singles out a position and prepares a coherent superposition of free energies. The resulting return contains a permanent component at √5. No coherent splitter, periodic force, timed gateway check or reference path is used in this launch. The return also contains other frequencies.

There is a direct lower bound on that component. Number the local ports from 0 to 11 and let the entries of an integer vector at ports 4 through 7 be (1,1,−1,−1), with zero elsewhere. On each tested antipodal cell, let B be 69I − 2H on the **full board** and v = Bu. Integer matrix operations give Bv = 5u, with inner products u·u = 4, u·v = 4 and v·v = 20. The vectors √5u ± v are eigenvectors whose energies differ by √5. The squared overlaps with port 4 are (1 ± 1/√5)/8. In a return to the launch port, every energy level enters the amplitude with its nonnegative spectral weight. Consequently any contributions at the same frequency add. The coefficient of cos(√5 t) in the return probability is at least **1/40** on both boards for arbitrarily long free evolution. After four native exchange quarter-steps, the origin saw 88,166 and 116,409 returns out of 262,144 respective launches. The finite return counts sample this multi-frequency motion; the exact lower bound follows from the invariant subspace.

## Reading a single frequency

An observer can filter the vacancy locally from a sharp input. Isolate its cell, allow the load to move, and check the gateways after each of 1,024 free pulses. Accept when no vacancy is observed on any gateway, then open every contact for the entire release. The success fraction is 5/8 on the two antipodal-gateway preparations and 1/3 on the 20-cell control. Gateway monitoring and the pulse schedule are supplied. A shorter preparation or a different untouched launch diameter is a separate test of frequency rigidity.

The isolated cell has two antipodally odd adjacency eigenspaces at +√5 and −√5. The antipodally even nonconstant eigenspace has eigenvalue −1. In the retained compact part of this source the corresponding free vacancy energies are (69 − √5)/2, (69 + √5)/2 and 35. At a port and its antipode the odd components change sign and the even components keep their sign. Summing the two terminal probabilities cancels both odd–even beat terms. The odd–odd beat survives at √5. The finite 1,024 checks leave a residual outside this compact space. Its norm is below 6 × 10⁻¹⁵ on the two nominal releases and 9 × 10⁻⁹ with the shorter preparation. The simulations evolve the actual accepted vector on the full board. Pulse length can change the beat phase.

On fresh enrollments at each of 129 appointments, read the two positions and a separate material exchange timer. Half of the timing records determine the sinusoidal fit; the later 64 check its forecast. The time estimate comes from timer swaps during one of eight equal subintervals. Each appointment uses 262,144 enrollments including failed vacancy sources. RMS errors are in probability per enrollment.

| Board and release | Frequency per timer gap | Fitted peak-to-trough probability | Held RMS error |
| --- | ---: | ---: | ---: |
| 80 cells | 2.2367 | 0.0998 | 0.00157 |
| 320 cells | 2.2321 | 0.0999 | 0.00145 |
| 80 cells, second diameter | 2.2365 | 0.0998 | 0.00176 |
| 320 cells, half-length preparation pulses | 2.2373 | 0.0994 | 0.00162 |
| 80 cells, half contrast at release | 2.2338 | 0.0995 | 0.00163 |
| 20 cells, three gateways with no antipodal pair | Unresolved | Below 0.001 | 0.00079 |
| 80 cells, exchange removed at release | Unresolved | Below 0.001 | 0.00106 |

## What fixes the frequency?

A separate energy instrument samples the same prepared state on fresh copies. An unsigned occupancy check of a randomly selected internal bond admits either endpoint. On that branch, controlling a full isolated contact exchange against idle and recombining the control reads the symmetric and antisymmetric load modes. Missing the bond and a failed source each have their own output. Multiplying the conditional antisymmetric fraction by 30 gives the mean exchange energy in the cell's thirty bonds. External bonds add zero for an exactly compact state, and the finite-preparation correction is bounded by the quoted residual norm and the sector bandwidth. The static energy of a vacancy on a five-contact port is 4κ above the filled host, from the native seams and degree reference. In the six-seam storage convention a separate stationary control of exchange energy 3 gave 34.9943 ± 0.0047 and 35.0106 ± 0.0060 in earlier material-clock transfer measurements on 20 and 320 cells. This checks the energy reference independently of the beat. The samples here superpose several energies; the new instrument estimates their mean. [Transfer instrument and physical reference](../row_04_compact/README.md#measuring-energy-by-moving-the-same-load)

Halving the contrast at release on the same 80-cell material gives **18.62 ± 0.02** in place of **34.62 ± 0.02**, while the measured clock frequencies agree within 0.003. Removing exchange at release stops the clock. This identifies the clock's physical origin as interference between two spatial exchange modes.

For any initial density matrix and fixed terminal effect in a closed sector, each term in the detection probability oscillates at an energy **difference**. The vacancy Hamiltonian is the graph Laplacian divided by two plus the native site potential, 32 on degree-five ports and 48 on degree-six ports. Its row-sum energy enclosures give [32,54] on each board. Thus **every** single-vacancy internal line lies at or below 22, for any source or fixed terminal detector on that sector. The positive excess has energies in [48,69] and internal frequencies at most 21. Neither blank-collar bare-defect sector can carry an internal line as fast as its filled-host formation energy on these boards. A stationary energy eigenstate gives constant probabilities. A more general excitation with dynamical memory or a physical inertial mass determined independently of the host and storage conventions needs its own test.

The [formal rest-phase result](../../../../../Lean/Geometry/InternalClockRestFrequency.lean) gives the derivative of a declared plane-wave phase. Linking this phase to an observable rest-energy oscillation requires a suitable pair of interfering physical components. The [upstream M10 clock entry](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L437) describes a pumped construction with a supplied phase resource. The sharp spatial localization here supplies energy coherence at preparation; an initially energy-diagonal mixture remains stationary under the free law. This study introduces no mass in laboratory units.

## Replay

From the repository root, using its scientific Python environment, run

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_03_clock
```

The replay regenerates the sources, both native measurements, the sharp-return panels and the sector bounds. The independent checks reconstruct the constant-frequency prediction from the stored counts and verify the integer-pattern identities on the local incidence matrix. [Checks](test_clock.py)

## A neutral creation reference

An exact native-operation check extends the reference question beyond one fixed terminal detector. Echoes, free contacts and sharp memory/load pointers preserve the filled, zero-contrast subspace, so finite adaptive programs in this grammar neither create a neutral pair there nor read its coherence with a filled reference. An extra coherent analyzer can resolve a 32-unit phase. Separately, a uniform charge-reference shift is invisible to complete fixed-total-load transcripts, including compensating references and classical feedback. The neutral creation energy fixes a sum, without selecting the two charged formation energies. [Native instruments, boundary controls and proof](neutral_reference/README.md)

## A free-end clock on a collision daughter

A three-dimer receiving wing can resolve all three of its one-quantum energy costs by comparing successive current signs at its free end. The same three frequencies appear on an ordered sham. A sharp one-quantum herald followed by silent waiting gives zero mean terminal current in the native material; inserting the initial current read removes about half the heralded excitations from that sector. This distinguishes local spectroscopy of the material from a freely arriving excitation carrying its own rest clock. [Counts, controls and proof](free_end/README.md)

With a supplied coherent route and a ground-prepared memory, a second observer can measure the three relative energies on already heralded arrivals while preserving the one-quantum flag. Tuning the memory to the middle mode shifts the pointer by that reference energy. [Measurement against idle storage](free_end/ROUTING.md)
