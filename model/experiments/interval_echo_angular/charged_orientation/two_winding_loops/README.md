# Quantum information arrives at one fixed material port

A source owner prepares two possible histories of a single vacancy circulating its three-port loop. Both histories have the vacancy at the same source gateway. A second owner has its own prepared loop. After free Interval Echo evolution, the receiving owner first checks for the vacancy at a **single gateway**. On arrival, two states of its three local seam memories retain quantum information about the source. All three source heads are discarded. At exposure 16 the arrival probability ranges from **24.04% to 43.72%** across input states. The conditional states carry information about a winding history in the prepared apparatus. [Actual six ports and pulses](measure.py) · [All simulated counts](results.json)

## The apparatus and the measurement

The source triangle is 1–0–5–1 and the receiver triangle is 82–74–78–82. The blank 1–82 contact joins them. On a magnitude-one directed seam, a permitted vacancy hop changes the direction coin. Sending the vacancy once around a triangle changes all three direction coins from zero to one and returns the vacancy to the same material port. A balanced blank-dimer exchange, a relative private phase and a native echo write each bond in its specified contrast-sign parity. Both sharp load heralds can be corrected to the same direction. The first source bond uses odd contrast-sign parity, and the other five use even parity. The all-even comparison has identical geometry and the same registered free law.

Initially the vacancy visits the four degree-five loop tips to erase their private contacts. Eight separate shutter writes isolate the compensating excess on the 12–13 dimer. The board then has a single vacancy among its six active laboratory ports. All other exterior ports outside the compensating dimer stay filled. Eighteen off-loop seams carry unit-pair silent marks during flight. Two of these provide addressed phase access. The emitter keeps the 0–7 contact blank during preparation, uses full swaps to lend its charged tip coordinate to the ordinary private reference at 7 and returns it. Both winding words have the vacancy on gateway 1 after preparation, allowing the emitter to close 0–7 by an echo on two unit loads.

After free flight a sharp vacancy read at port 82 heralds arrival and preserves coherence between the local head words. On success both ends of the receiving 74–75 phase contact are filled. Its inverse unit-pair echo restores the blank contact before any relay pulse. Two addressed full directed exchanges unwind the receiving loop down to two adjacent material coordinates. Native swap, private dwell and return give the relative phase. A partial exchange and sharp gateway/tip load reads resolve the chosen Pauli axis. All three axes are determined through receiver-local controls; the three source heads and the rest of the exterior are traced. The program verifies the shutter inverse on the 99-state native seam, and the full phase relay restores the external port for either tip occupation. Every experiment counts the two local outcomes and a missed herald on original starts. [Preparation, exterior separation and local readout](measure.py)

## Why the two preparations differ

For the uniform all-even bench, reversing the source triangle swaps the two initial winding words and leaves the Hamiltonian invariant. On the accepted receiver port that reversal acts only on the discarded source heads. The success map therefore has the same output for an input density matrix and its conjugate by the logical Pauli X. Averaging the two inputs is dephasing in the X basis. The success map factors through this two-outcome classical measurement and breaks entanglement at every time. At the gateway the two discarded source histories span a doublet with one dimension of each reflection parity. A sharp parity input consequently gives a rank-at-most-one receiver state on success, which explains the zero eigenvalues in the control witness.

Writing odd contrast-sign parity on the first source contact breaks this symmetry with a native preparation. For each final source winding, the column pairs of the physical propagator define one conditional amplitude matrix. Summing their contributions gives the success map after discarding the three source heads. At 16, applying this map to half of a normalized Bell state gives trace **0.33878329**. Its partial transpose has minimum eigenvalue **−0.00810956** per original start, or **−0.02393730** after normalizing on success. An integer Taylor evaluation bounds a rational witness below **−0.0081** before normalization. The bound includes the propagation approximation, projection, discarded heads and reference transposition. A local miss flag can be appended on its own orthogonal outcome. Its inclusion preserves the negative success block. [Rational certificate](certificate.py)

Observers reconstruct the same linear map from six native input states and three receiver settings, without measuring a remote reference. Each setting uses 16,777,216 simulated starts, counting both outcomes and every miss. The minimum partial-transpose eigenvalue from the records is **−0.00814908** with a simultaneous 99% error radius of **0.00207398**. The all-even bank gives **−0.00003753** with the same radius and its exact value is zero. These records distinguish transfer of quantum information from retaining only a classical record of input parity. They measure an arrival-conditioned channel with input-dependent success. [Counts for every basis](results.json)

## Equal local energy and free dwell

Let P retain the fixed receiving gateway and both winding words, with both source histories allowed. In exchange units the native integer matrices give
[
 PHP=\tfrac92P,\qquad PH(I-P)HP=\tfrac{65}{4}P .
]
Every normalized received code state thus has mean formation energy 4.5 over the matched filled heads and variance 16.25. Starting on that gateway, free evolution for another quarter unit retains between 52.70% and 88.05% of the code on the asymmetrically prepared bench. Its complete 16-word generator has simple spectrum, verified by the square-free characteristic polynomial. The received winding pair provides a local quantum readout; its energy spread and leakage prevent interpreting it as a stationary spin doublet. Assigning an electron would additionally require a physical spatial-rotation response on the same state, exchange statistics and a rest-energy measurement. [Detector frame under relabeling](FRAME.md)

From the repository root, run the measurement with one numerical thread and check the explicit trace and two phase shutters.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_angular.charged_orientation.two_winding_loops.measure
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/two_winding_loops/test_channel.py
```
