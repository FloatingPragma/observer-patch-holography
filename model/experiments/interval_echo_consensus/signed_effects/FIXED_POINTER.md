# A fixed neutral pointer measures either material sign

A neutral material dimer prepared with loads (2,0) records the source address for either vacancy or excess. Its preparation and pulse word are independent of the incoming sign. A native phase correction also preserves coherent superpositions whose total charge is balanced by an untouched reference. This removes the matched-sign pointer-preparation assumption from the two-populated-rail experiment. Protected contacts, a unit-filled host and the final sharp occupation instrument are supplied. [Pulse word and records](fixed_pointer.json) · [Native implementation](fixed_pointer.py) · [Owning checks](test_fixed_pointer.py)

## Physical preparation and domain

The source has one vacancy or excess at port 0 or 1. Two fixed neutral pointers have park/cross ports (6,10) and (4,5), each initialized to (2,0). Port 7 is a host of load one. Every other contact is parked. The pulse settings do not branch on the sign, a source census or an earlier outcome.

A classical unknown sign can be described by any mixture of the two source sectors. The pointer and host preparation are the same in each. Coherence between those sectors requires a charge reference. For the coherent test, port 11 carries load two with a vacancy and zero with an excess. The source and reference together have fixed total charge. This reference is part of the input, and is untouched throughout record acquisition. It supplies no sign instruction to the pointer or controller. A separate native analyzer couples to it only in the phase-control experiment.

Propagation retains the complete conserved multiset sector of the eight ports: three zeros, two ones and three twos, dimension 560 and total load eight. This is a reducing sector inside the **1,107-dimensional total-load-eight sector**, and inside the full ternary space of dimension 6,561. Every primitive contact swaps endpoint occupations or phases them, and every reference is diagonal, so the multiset is preserved. Mixed endpoint words (0,2) and (2,0), temporary escape from the dimer code and all 560 placements are retained. No signed hard-core projection is made during the pulses.

An independent second implementation propagates all 6,561 ternary configurations from four coherent input columns using local 9-by-9 contact and 3-by-3 reference exponentials. It imports the declared word and constructs its own source and ideal record maps. Both two-record histories agree, including all off-code rows and one common phase across signs. [Independent full-space regression](test_fixed_pointer_full_space.py)

## Refocusing all ternary endpoint pairs

The native blank-memory contact is

\[
h=(I-\mathrm{SWAP})/2+8(n_a-n_b)^2,
\qquad V(a,b)=8(a-b)^2+\tfrac12\mathbf1_{a\ne b}.
\]

The three unequal swap blocks are (0,1)/(1,0), (1,2)/(2,1) and (0,2)/(2,0). The first two have diagonal cost 17/2; the last has 65/2. A parity reference \(Z=\exp(-i\pi n_a^2)\) flips exchange in the odd-difference blocks. For contact duration \(\tau\), the chronological contact, parity, contact, parity word gives an operator \(E_\tau\). It is scalar on the odd blocks and retains exchange on the extreme block.

The reference \(T=\exp(-i\pi n_a^2/4)\) flips exchange on (0,2)/(2,0). Conjugating one of two \(E_\tau\) words by \(T\) cancels that exchange as well. The complete chronological word is

\[
E_\tau,\ T,\ E_\tau,\ T^*,
\]

and implements exactly \(\exp(-4i\tau V)\) on all nine endpoint words. The parity waits are \(\pi/8\), the \(T\) wait is \(\pi/32\), and its positive inverse wait is \(63\pi/32\). No arbitrary controlled unitary replaces this native pulse sequence.

For the source/pointer cross contact choose total contact dwell \(4\tau=2\pi/65\). The neutral pointer's cross rail is zero or two. The difference of the two costs is

| Source cross-port load | \(V(n,2)-V(n,0)\) |
| --- | ---: |
| 0, vacancy | \(+65/2\) |
| 1, host | 0 |
| 2, excess | \(-65/2\) |

Both marks therefore give a Ramsey flip, while a host does not. On the neutral dimer, reference \(\pi/64\), contact \(\pi/2\), reference \(\pi/64\) implements Hadamard up to a common scalar. Applying this basis word, the diagonal interaction and the actual positive-wait inverse gives the pointer action

\[
n=0:\ X,\qquad n=1:\ gI,\qquad n=2:\ -X,
\qquad g=e^{-17\pi i/65}.
\]

This bare word measures a classical unknown sign, but its relative sign and host phase change coherent source states.

## Correcting the branch phase without reading the sign

A source reference wait \(\pi/32\) gives phases \((1,e^{-i\pi/4},-1)\) on loads (0,1,2). It cancels the excess minus sign. Apply the same fully refocused diagonal interaction between source port 0 and the fixed load-one host at port 7, with total dwell \(133\pi/2210\). It gives phase \(h=e^{-133\pi i/260}\) on either mark and one on a host. Since \(g e^{-i\pi/4}=h\), the corrected action is

\[
|n\rangle|0_{\rm pointer}\rangle\longmapsto
h\,|n\rangle|\mathbf1_{n\ne1}\rangle_{\rm pointer}.
\]

The same scalar holds for all three source loads. All encoded pointer inputs and both source signs are checked. The copy is nondemolition on the addressed observable and preserves the relative charge phase. It copies an orthogonal address record, without cloning an arbitrary quantum state.

The phase test has its own native preparation and analyzer. On seam (0,11), the neutral-dimer basis word prepares an equal superposition of vacancy/reference-two and excess/reference-zero at source port 0. After one pointer copy, the inverse basis word returns the reference to load zero with probability one for the corrected word, and to load two with probability one for the bare word. All three reference-load outcomes are retained. These are exact ideal probabilities; the receipt retains unnormalized floating propagation centers and their roundoff. Both words give the same sharp address record on the classical sign inputs. This control tests the coherent phase left by the implementation.

## Complete sequential instrument

For each pointer, run the specified three-pulse source basis word, the corrected material copy and the inverse source word. No source or pointer read occurs between the two acquisitions. One final sharp read resolves the occupations of four pointer ports and the host. With each pointer word \(q=n_{\rm park}+3n_{\rm cross}\), the record is \(q_1+9q_2+81n_{\rm host}\). All 243 values are retained, including malformed pointers and a disturbed host.

On the encoded four-dimensional source/reference space, the branch maps are exactly one common phase \(e^{-133\pi i/130}\) times the direct sequential maps. Both sign-diagonal and cross-sign matrix units are checked. Inserting a native source swap between copies preserves the first measurement effects and yields the same distinct second-read statistics as the projective control. Twenty banks retain 81,920 synthetic starts, sampled as a binary tree of first-read and conditional second-read probabilities. Each tree probability is rounded to an integer weight on the \(2^{-40}\) grid, inside the \(3\times10^{-12}\) propagation budget, and each draw compares a PCG64 integer below \(2^{40}\) with that weight, so the counts replay identically on every platform. They include definite-sign coherent source preparations and two states coherent across the balanced charge sectors.

Leaving out the extreme-block refocusing pulses produces escape probability 0.0011670913 for either sign at source port 0. Source-port-1 inputs expose only odd differences and retain zero escape under this control. This tests the mixed-sign collision omitted by a binary-only derivation. The final exact refocused word has no code leakage. All control inputs and outcomes are retained.

There are 66 native waits, or 67 with the source swap. Contact generators have norm at most 33 on the full ternary sector and references at most 32. A copy has ten contacts and seventeen references. Each source rotation/undo adds four contacts and two references. Preparation half-trace error \(\epsilon_p\) and absolute error \(\delta\) on every duration therefore give complete two-record classical/quantum output half-trace error at most

\[
\min(1,\epsilon_p+2173\delta),\qquad
2173=2[(10+4)33+(17+2)32]+33.
\]

This is a conditional all-outcome bound; no herald probability divides it. Physical duration errors, trusted preparation and final sharp-pointer fidelity are not certified here.

The result implements a fixed neutral ancillary preparation on the two populated source rails. A general twelve-position source instrument, arbitrary adaptive schedules, preparation of arbitrary physical inputs, the final projective read, physical fidelity and external custody retain their own obligations. It does not retire PR-03, PR-64 or PR-65.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.signed_effects.fixed_pointer
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects/test_fixed_pointer.py model/experiments/interval_echo_consensus/signed_effects/test_fixed_pointer_full_space.py
```

Replay is read-only. Explicit `--write` creates a missing receipt exclusively. Source hashes bind producer-side software and supply no external authentication or prospective experimental claim.

A neutral output dimer also admits a native basis turn and echo that deposit its bit in a protected parity of the seam memory. Free board evolution preserves the joint probabilities of four such records for arbitrary material motion. Decoding consumes a fresh or restored pair of unit loads and a final sharp read. [Free storage and disturbance controls](COLLAR_STORAGE.md)
