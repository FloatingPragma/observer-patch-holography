# A standing wave written by an occupied seam

A native echo on a seam containing one load produces a stationary direction record entangled with the mobile material. On both the 20- and 80-cell boards, reversing the phase on just one occupied endpoint moves the dominant **measured** clock bin from −13.83037 to −5.89049 for the even source. All sharp input position probabilities are identical. The odd source moves from −14.82439 to the same −5.89049 bin. These energies are relative to a uniform blank load plus the conserved energy of a separated memory in the same direction parity. [Frozen protocol](PROTOCOL.md) · [Full native probabilities and simulated counts](results.json.gz)

| Cells | Direction parity | Ordinary clock peak, counts | Reversed clock peak, counts | Reversed counts in held $n=2$, $n=3$ windows |
| --- | --- | ---: | ---: | ---: |
| 20 | Even | 3,627,485 | 3,042,988 | 36, 15 |
| 80 | Even | 3,628,568 | 3,045,228 | 52, 13 |
| 20 | Odd | 4,122,444 | 3,602,517 | 15, 7 |
| 80 | Odd | 4,122,559 | 3,604,113 | 18, 7 |

Each clock bank retains 4,194,304 simulated starts. The ordinary peak fixes the ground binding. Dividing its energy by four and nine predicts −3.45759 and −1.53671 for even parity and −3.70610 and −1.64715 for odd parity. All counts within 0.1 exchange units on either side enter the held windows, including the finite clock's sidelobes. On 20 cells the even-parity bound energies are −13.829117, −6.136163, −5.888926 and −5.656226. The ordinary source places 99.452% of its spectral weight on the deepest state. The endpoint reversal leaves 0.091% on it and puts 90.176% on the −5.888926 state. The −6.136163 and −5.656226 states are dark to the ordinary occupied source. The odd-parity ground energy is −14.824011. [Spectral blocks and detector program](measure.py)

## Source, clocks and receivers

Start with a sharp amount-one load and blank memories. A shuttered native exchange at the core for time $\pi/2$ splits the load between its two endpoints. A private degree-reference pulse of duration $3\pi/16$ or $\pi/16$ sets the standing dimer even or odd. One native echo writes its signed direction and entangles the coin with position. A separate $\pi/8$ reference pulse on one endpoint reverses the relative orbital phase after writing. Each registered seam has a degree-five endpoint with a private reference. The shuttered source operations, phase controls, prepared uniform-load reference, synchronized clock with coherent conditional waits and the terminal signed-record read are supplied apparatus.

On independent copies, two blank-seam swaps park the load off both core endpoints, retaining its coin correlation. A shuttered comparison of the resulting empty core with a bypass, in two phase quadratures, reads 8.00609 and 8.98832 in even and odd parity. The uniform-load clock reads 48.02127. Each has simultaneous 99% allowance ±0.12237 exchange units. The summed native core and free-load reference has a conservative allowance of ±0.24474. Even that whole shift leaves each held target inside the spectral gap. The spectral interferometer compares coherent evolutions of the source and physical reference copies directly. Its 4,096 forward waits are spaced by one quarter and its frequency period is $8\pi$. The measured sector spectra fit inside one period. All material positions remain included at waits 0, 1, 8 and 32. Each owner reads the signed load differences along eleven blank seams in a tree disjoint from the written seam, restores the material between reads, and recovers its port word from those differences. With at most one load on all twelve ports, the minimum relative level is the empty level. The controller combines the local records to report the terminal address. [Preparation, physical reference evolution and local address construction](core.py)

For the native even source after 32 exchange units on 20 cells, the two cases return 4,144,257 and 4,021,036 loads on the original endpoints. Switching off contrast during flight leaves both even-source endpoint laws identical, with probability 0.24600 on 20 cells and 0.33697 on 80. The even changed-law spectrum has no negative level. Odd parity retains the exchange energy of its memory, supporting one weak sector-relative bound level at −0.06869 on 20 cells and −0.06894 on 80. The blank and vacuum-written controls, both phase settings and all waits are retained. [Uniform, vacuum-written and changed-law controls](results.json.gz)

## Which atomic spectra are excluded?

Against the cheapest separated reference, a uniform blank load and a zero-energy core, the full native Hamiltonian for **any coherent state of one eleven-state seam plus one unit load** has exactly six negative energies counted with multiplicity, all in $[-8,-5.5]$. The complete spectral gap $(-5.5,0)$ excludes the second and third hydrogenic bindings of any lowest line in this family. In each occupied parity sector referred to its *own* higher separated memory energy, there are exactly four negative levels, all at or below −5.5. This accounts for the measured −13.83 and −14.82 ground bindings as well as the held empty windows. Both results hold on the 20-, 80- and 320-cell boards. Only the first two sizes have sampled flight and clock records. Evidence for discrete standing waves is Demonstrated in this model. The full OW-31 hydrogen-like criterion remains Partial on the local scorecard; the OpenWave M10 cell remains untested. [Finite-board proof](PROOF.md)

From the repository root, with the scientific Python environment active, replay using one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_31_orbitals.occupied_write
nice -n 15 python -B -O -W error -m "$p.measure"
nice -n 15 python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_31_orbitals/occupied_write
```
