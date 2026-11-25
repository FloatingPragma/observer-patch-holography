# Interference and Bell correlations in local charge records

Two sharp loads start at ports 1 and 3 of cell zero on the registered 20-cell Interval Echo board. All 630 seams act for eight exchange times. No source pulse entangles the loads and no guide routes them. Observers in cells 6 and 4 then isolate their local port pairs and ask whether exactly one load arrived. Among 268,435,456 simulated launches, both answered yes 533,643 times. Their later signed records on these trials give a Clauser–Horne–Shimony–Holt (CHSH) score of **2.7005**. Its simultaneous 99% lower bound is **2.6604** under the simulated independent trials. A Bell-local model with independent settings has an upper bound of 2. [Protocol](PROTOCOL.md) · [Experiment](experiment.py) · [Counts](results.json)

## What each observer does

Each observer uses the local echo, reads one bit for absolute contrast one, and undoes the echo. This measurement keeps the two possible occupied ports coherent. Both flags are recorded before fresh analyzer settings. Each wing has a supplied independent setting apparatus equivalent to a sharp one-load dimer, a native half contact and a signed pointer. Its initial state factors from the host and from the other setting apparatus. The simulator samples these independent setting distributions in batches. The selected measurement uses a degree-five port as a phase reference, rotates the local seam, and reads positive, negative, or quiet. The two local operations commute. Neither observer receives the other setting. Later comparison of the two observers' records gives the four correlations. The receipt keeps every capture category, setting pair and local outcome.

The captured state was found in an exploratory calculation over 870 pairs of receiver seams at four times. Only the detector geometry and the free-flight appointment were selected this way. Nine independent simulated calibration tables supplied the measurement axes. A fresh acquisition scored those axes. The 99% intervals cover simulated counting noise per panel. They quantify neither uncertainty in the law nor the frequency of violations on other preparations. With the same source, rails and axes, the capture scores at times 7.5 and 8.5 were 0.6301 and 2.1934 respectively. At the first appointment the selected detector misses the violation. The second gives a lower bound of 2.1582.

| Experiment | Both flags positive | Score across all launches | Score on joint captures |
|---|---:|---:|---:|
| Free flight, time 8 | 0.1988% | 1.6311 | 2.7005 |
| Same flight, measure which rail and forget | 0.1982% | 1.6550 | −0.0536 |
| Exchange alone, same apparatus | 0.1160% | 1.7105 | 1.2368 |
| Addressed split and swap routing, five cell steps | 50.000% | 2.4143 | 2.8284 |
| Same routing, measure which rail and forget | 50.005% | 1.7067 | 1.4136 |

Quiet always scores +1 in the all-launch column. On the free source, varying the four local two-rail rotations with this scoring convention cannot close the all-launch gap. Write the two single-capture Bloch vectors with their unconditional weights as $u,v$, and the joint failure and success fractions as $p_{00},p_{11}$. The quiet–quiet sector contributes $2p_{00}$. On each single-capture sector one wing has fixed output, giving absolute contributions bounded by $2\|u\|$ and $2\|v\|$. The joint-capture sector contributes at most $2\sqrt2\,p_{11}$. The nine calibration tables give a simultaneous 99% upper bound of **1.925** for their sum at time 8. This bound concerns these receivers and this binary scoring rule. For the routed apparatus the coherent source occupies the one-load-per-wing sector with probability one half. The remaining branches give two quiet readings. This gives the exact all-launch prediction $1+\sqrt2$.

## What happens if the gateways settle a record during flight?

At time four both source gateways publish their signed contrast and restore the material through the inverse echo. All host contacts then act for another four units. Keep every pair of midpoint signs, including two quiet signs. The same receivers jointly capture 427,481 of 268,435,456 launches. Their score is 0.7467 at the previously successful settings, or 1.6594 across all launches. This measurement distinguishes histories that previously interfered.

Fresh local calibration and a separate acquisition give **0.9800** on the pooled captures after choosing new analyzer axes. For a captured pair of two-rail qubits with correlation matrix $T$, the best CHSH score attainable by rotating these projective meters is twice the square root of the sum of the largest two squared singular values of $T$. This follows by writing the sum and difference of the two unit Bob directions as orthogonal vectors and optimizing the two Alice directions and the relative angle. The nine measured correlations give a simultaneous 99% upper bound of **1.011** for this class. The bound adds twice the Frobenius norm of the nine Hoeffding error limits to the fitted maximum. It applies when all midpoint sign pairs are pooled. These records distinguish a coarse local capture that preserves the two rail alternatives from a specified intermediate publication that suppresses their Bell signal. The experiment schedules the publications. Spontaneous selection of a measurement regime by the clock and repair process remains open.

## One freely propagating load interferes

Start at port zero and leave every seam open for one unit. Ports 5 and 7 serve as the two apertures. After another unit of free flight the preselected return detector sees **5.684%** of all launches with both apertures and zero added phase. The two single-aperture returns sum to **2.899%**. Reading the aperture and forgetting the result gives **2.845%**. A half-turn of relative phase brings the return to **0.013%**. The 240 terminal positions and absorption are all retained. With the apertures removed altogether, the return is **11.666%** for uninterrupted flight and **13.712%** after an unread midpoint position measurement. This is a finite spatial interference experiment with a return detector.

## What the simulation assumes

The registered interaction evolves coherent amplitudes, and the local meters use Born probabilities and projective state updates. A public flag is classical once read. The premeasurement amplitudes can still interfere within the recorded sector. During these free flights all collars are blank and the inherited free Hamiltonian acts. The clock's signed first-return service is idle. The experiment demonstrates what observers of this coherent model can record. Deriving that coherent instrument and its outcome probabilities from a lower-level irreversible repair process remains a separate task. Graph separation and local isolation are implemented here. A physical spacelike identification of the measurements would require independent clock and range evidence. [Candidate specification](../../candidate_laws/interval_echo/SPEC.md) · [Finite event-algebra theorem surface](../../../Lean/EventAlgebra/QuantumAdequacySurface.lean) · [Diagonal-state bound](../../../Lean/EventAlgebra/ProductSplitSeparability.lean)

Run from the repository root with the model dependencies installed.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_quantum.experiment
```

[Two independent sources and a middle Bell measurement](swapping/) give all-launch CHSH and independent-source violations after the three stations join their records.
